#' @title methPLIER
#' @description methPLIER
#'
#' @author Ken Takasawa

#### Slightly modified multiPLIER function -----------------------------------------------------------
#' PLIERNewData
#' 
#' @description A wrapper function for applying PLIER to a data set.
#' We use the following genesets that come with PLIER: bloodCellMarkersIRISDMAP, svmMarkers,
#' and canonicalPathways.
#' We set the k parameter for the PLIER model by
#' identifying the number of "significant PCs" with PLIER::num.pc and then
#' using sig PCs * 0.3. This is consistent with recommendations from the PLIER authors.
#' @references [Mao, W., et al., Nat. Meth. 16(7):607-610, (2019)]
#' @references [Taroni, N. J., et al., Cell Syst. 8(5):380-394.e4 (2019)]
#'
#' @import PLIER
#' @param D a DNA methylation matrix, rows are genes, columns are samples
#' @param C knowledge matrix
#' @param seed an integer to be supplied to set.seed() for reproducibility purposes, default is 12345
#' @return plier.res: output from PLIER::PLIER()
#' @export

PLIERNewData <- function(D, C = 'all', seed = 12345) {
  set.seed(seed)

  if (!is.matrix(C)){
    # load PLIER pathway and cell type data
    data(bloodCellMarkersIRISDMAP, package = 'PLIER')
    data(svmMarkers, package = 'PLIER')
    data(canonicalPathways, package = 'PLIER')
    data(oncogenicPathways, package = 'PLIER')

    # combine the pathway data from PLIER
    C <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers,
                                     canonicalPathways, oncogenicPathways)
  }

  # transform to matrix
  d <- d %>% inner_join(tibble(Gene=rownames(C)))
  cm.genes <- d %>% pull(Gene)
  d <- d %>% dplyr::select(-1) %>% as.matrix()

  # row normalize
  mat.norm <- rowNorm(d)

  # what should we set the minimum k parameter to in PLIER? estimate the number
  # of PC for the SVD decomposition
  set.k <- PLIER::num.pc(mat.norm)

  # PLIER main function + return results
  plier.res <- PLIER::PLIER(mat.norm, all.paths[cm.genes, ],
                            k = round((set.k + set.k * 0.3), 0), trace = TRUE)

  return(plier.res)
}


#' GetOrderedRowNorm
#' 
#' @description Given a DNA methylation matrix of data that was not used to train a PLIER model
#' and the output of PLIER::PLIER, row-normalize, remove genes not in the
#' training data, set missing genes to zero (the mean), and reorder to match
#' plier.model$Z
#' This makes the input DNA methylation data suitable for projection into
#' the training latent variable space (see GetNewDataB) and for evaluating
#' reconstruction (MASE, correlation)
#'
#' @param D a DNA methylation matrix, rows are genes, columns are samples
#' @param plier.model PLIER results, output from PLIER::PLIER()
#' @return meth.cg: a row normalized (z-scored) DNA methylation matrix that matches the Z matrix of the PLIER model -- ordered to match, contains the same genes
#' @export

GetOrderedRowNorm <- function(D, plier.model) {
  # first, row normalize the DNA methylation data
  meth.norm <- rowNorm(D)

  # get Z matrix from PLIER model
  z.mat <- plier.model$Z
  # what genes were used in the model?
  genes.in.model <- rownames(z.mat)

  # get the genes that are common to the PLIER model and the new expression
  # data
  meth.cg <- meth.norm[which(rownames(meth.norm) %in% genes.in.model), ]

  gene.ord.chk <- all(rownames(meth.cg) == rownames(z.mat))

  if (gene.ord.chk) {
    return(meth.cg)
  } else {
    stop("Something went wrong -- Z matrix gene order doesn't match with
         the ordered DNA methylation matrix")
  }

}


#' getNewDataB
#' 
#' @description Apply a previously computed PLIER to a new dataset to get the LV x sample
#' matrix (B)
#' see main PLIER function:
#' https://github.com/wgmao/PLIER/blob/a2d4a2aa343f9ed4b9b945c04326bebd31533d4d/R/Allfuncs.R#L227
#'
#' @import PLIER
#' @param D a DNA methylation matrix, rows are genes, columns are samples
#' @param plier.model PLIER results, output from PLIER::PLIER()
#' @return meth.new.b: a matrix that contains the values of each latent variable for each sample from the new dataset (D),
#' @export

getNewDataB <- function(D, plier.model) {
  # the error handling for reordering failing is now in GetOrderedRowNorm
  ord.rownorm <- GetOrderedRowNorm(D, plier.model)

  # get Z matrix from PLIER model
  z.mat <- plier.model$Z

  # get LV by sample (B) matrix from the DNA methylation using PLIER model
  meth.new.b <-
    solve(t(z.mat) %*% z.mat + plier.model$L2 * diag(ncol(z.mat))) %*% t(z.mat) %*% ord.rownorm

  # add in the rownames from the PLIER model B matrix
  rownames(meth.new.b) <- rownames(plier.model$B)

  # return B matrix
  return(meth.new.b)
}

#### My difined function -----------------------------------------------------------
#' rowNorm
#' 
#' @description Row normalization to apply PLIERNewData function.
#' Subtracting mean value of the row from data, and dividing by standard deviation value of the row.
#'
#' @param x a DNA methylation matrix, rows are genes, columns are samples
#' @return x: Normalized data matrix
#' @export

rowNorm <- function(x){
  s = apply(x, 1, sd, na.rm = TRUE)
  m = apply(x, 1, mean, na.rm = TRUE)
  x = sweep(x, 1, m)
  x = sweep(x, 1, s, "/")
  x = ifelse(is.na(x), 0, x)
  return(x)
}

#' getDmatrix
#' 
#' @description Getting D matrix for applying to methPLIER
#'
#' @import tidyverse
#' @param data data table of  DNA methylation matrix.  Rownames are 'TargetID', Columns are 'Sample'
#' @return D: D matrix, row is gene, column is sample
#' @export

getDmatrix <- function(data){
  if (!exists('pca.w')){
    data('pca.w')
  }
  if (!exists('pca.ord')){
    data("pca.ord")
  }

  D <- data %>% as_tibble(rownames = 'TargetID') %>%
    pivot_longer(-TargetID, names_to = 'AccessionNo', values_to = 'methylation') %>%
    dplyr::inner_join(pca.w) %>%
    dplyr::group_by(id.gene, PC, Gene, TargetID) %>%
    dplyr::mutate(methylation = dplyr::if_else(is.na(methylation),
                                               median(methylation, na.rm = TRUE),
                                               methylation)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(methylation = dplyr::if_else(is.na(methylation),
                                               0, methylation)) %>%
    dplyr::group_by(AccessionNo, PC, id.gene, Gene) %>%
    dplyr::summarise(value = sum(weight * methylation)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = AccessionNo, values_from = value) %>%
    dplyr::full_join(pca.ord, .) %>%
    dplyr::ungroup()
  id <- D %>% dplyr::pull(Gene)
  D <- D %>% dplyr::select(-id.gene, -PC, -Gene) %>%
    dplyr::mutate_all(~ dplyr::if_else(is.na(.), 0, .)) %>%
    as.matrix()
  rownames(D) <- id
  D <- (D + 5) / 10
  return(D)
}


#' getGenesInLV
#' 
#' @description This function gets list of genes of LV
#'
#' @import tidyverse
#' @param methPLIER methPLIER: methPLIER model
#' @param LV LV number of to be obtained pathway
#' @param frac A number specifying what percentage of the total weight of the genes to be extracted (default is 0.8)
#' @return genes: vector of genes in LV
#' @export

getGenesInLV <- function(methPLIER, LV, frac = 0.8){
  genes <- methPLIER$Z[, LV] %>% tibble(Gene = names(.), value = .) %>%
    mutate(value = value /sum(value)) %>% arrange(desc(value)) %>%
    mutate(cum = cumsum(value)) %>% filter(cum <= frac) %>%
    distinct(Gene) %>% pull(Gene)
  return(genes)
}


#' getPathway
#' 
#' @description This function gets pathway of top genes of LV
#'
#' @import tidyverse
#' @importFrom  DOSE
#'   enrichDGN
#'   setReadable
#' @importFrom enrichplot
#'   dotplot
#'   pairwise_termsim
#'   treeplot
#' @importFrom graphics
#'   barplot
#' @param genes vector of genes in LV
#' @param showCategory number of plotting ontology in barplot
#' @return list: edo, edox, edox2, barplot, dotplot, treeplot
#' @export

getPathway <- function(genes, showCategory=15){
  # get entrez id table
  data(geneIdTable)

  entrez.id <- tibble(hgnc_symbol = genes) %>% inner_join(geneIdTable) %>% pull(entrezgene_id)

  library(DOSE)
  edo <- DOSE::enrichDGN(entrez.id, qvalueCutoff = 0.05)

  library(enrichplot)
  barplot <- graphics::barplot(edo, showCategory=showCategory)
  dotplot <- enrichplot::dotplot(edo, showCategory=showCategory)

  edox <- DOSE::setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
  edox2 <- enrichplot::pairwise_termsim(edox)
  treeplot <- enrichplot::treeplot(edox2, hclust_method = 'ward.D2')

  return(list(edo = edo, edox = edox, edox2 = edox2,
              barplot = barplot, dotplot = dotplot, treeplot=treeplot))
}


#' getDMGsInLV
#' 
#' @description Wrapper function of getDMGs for TopLVs
#'
#' @param methPLIER methPLIER object
#' @param D D matrix obtained by 'getDmatrix' function
#' @param LV integer of LV
#' @param annot sample annotation tibble dataframe. It have to contain following column,
#'   AccessionNo: the name of samples,
#'   cluster.sample: category of each samples.
#' @param frac A number specifying what percentage of the total weight of the genes to be extracted (default is 0.8)
#' @param threshold q.value threshold (default is 0.05)
#' @param method method of p.adjust (default is 'fdr')
#' @return result: result of getDMGs function
#' @export

getDMGsInLV <- function(methPLIER, D, LV, annot, frac=0.8, threshold=0.05, method='fdr'){
  genes <- getGenesInLV(methPLIER, LV, frac=frac)
  result <- getDMGs(D[genes, ], annot, threshold = threshold, method = method)
  return(result)
}


#' getDMGs
#' 
#' @description This function gets differentially methylated score genes
#'
#' @import tidyverse
#' @param D D matrix obtained by 'getDmatrix' function
#' @param annot tibble data frame as the sample annotation.
#    The attribute "cluster.sample" obtained by 'plotHeatmap' function can be used as arg.
#      Column_1: 'cluster.sample'
#      Column_2: 'AccessionNo'
#' @param threshold q.value threshold (default is 0.05)
#' @param method method of p.adjust (default is 'fdr')
#' @return result: tibble data frame
#' @export

getDMGs <- function(D, annot, threshold=0.05, method='fdr'){
  D[apply(D, 1, var) > 0, ] %>%
    as_tibble(rownames = 'Gene') %>% rowid_to_column(var = 'id.gene') %>%
    pivot_longer(cols = c(-Gene, -id.gene), names_to = 'AccessionNo') %>%
    inner_join(annot) %>%
    group_by(cluster.sample, Gene, id.gene) %>%
    summarise(value = list(value)) %>% ungroup() %>%
    pivot_wider(names_from = cluster.sample, values_from = value) %>%
    rowwise() -> D
  tx <- annot %>% distinct(cluster.sample) %>% pull(cluster.sample) %>% as.vector()
  script <- paste0('D <- D %>% mutate(p.value = t.test(`', tx[1], '`, `', tx[2], '`)$p.value)')
  eval(parse(text = script))
  D %>% ungroup() %>% mutate(q.value = p.adjust(p.value, method = method)) %>%
    filter(q.value < threshold) -> result
  return(result)
}



#' getDMPs
#' 
#' @description This function gets differentially methylated probes
#'
#' @import tidyverse
#' @param data data frame of DNA methylation data.
#'    The attribute "methylation" obtained by 'makeGvizObj' function can be used as arg.
#' @param annot tibble data frame as the sample annotation.
#'    The attribute "cluster.sample" obtained by 'plotHeatmap' function can be used as arg.
#'      Column_1: 'cluster.sample'
#'      Column_2: 'AccessionNo'
#' @param genes vector of genes
#' @param threshold q.value threshold (default is 0.05)
#' @param method method of p.adjust (default is fdr)
#' @return result: tibble data frame
#' @export

getDMPs <- function(data, annot, genes, threshold=0.05, method='fdr'){
  if (!exists('methProbeAnnot')){
    data("methProbeAnnot")
  }
  probe.tbl <- tibble(genesUniq = genes) %>% inner_join(methProbeAnnot, .) %>% unnest(data) %>% distinct()
  probes <- probe.tbl %>% distinct(TargetID) %>% pull(TargetID)

  data %>% as_tibble(rownames = 'TargetID') %>% inner_join(tibble(TargetID = probes)) %>%
    pivot_longer(-TargetID, names_to = 'AccessionNo') %>%
    inner_join(annot) %>% group_by(cluster.sample, TargetID) %>%
    summarise(value = list(value)) %>%
    pivot_wider(names_from = cluster.sample, values_from = value) %>%
    rowwise() -> d
  tx <- annot %>% distinct(cluster.sample) %>% pull(cluster.sample) %>% as.vector()
  script <- paste0('d <- d %>% mutate(p.value = stats::t.test(`', tx[1], '`, `', tx[2], '`)$p.value)')
  eval(parse(text = script))
  d %>% ungroup() %>% mutate(q.value = p.adjust(p.value, method = method)) %>%
    filter(q.value < threshold) %>% inner_join(probe.tbl) -> result
  return(result)
}



#' getTopLVs
#' 
#' @description This function get top-LVs
#'
#' @import tidyverse
#' @param cluster.LV tibble talbe of LV cluster obtained by 'plotHeatmap' function
#' @param col.1 column A
#' @param col.2 column B
#' @param threshold threshold for q.value filtering (default is 0.05)
#' @param method method for p.adjust (default is 'fdr')
#' @return res: tibble data frame
#' @export

getTopLVs <- function(cluster.LV, col.1, col.2, threshold = 0.05, method = 'fdr'){
  script <- paste0('res <- cluster.LV %>% rowwise() %>% mutate(p.value = t.test(`',
                   col.1, '`,`', col.2, '`)$p.value) %>% ungroup() %>% mutate(q.value = p.adjust(p.value, method = "',
                   method, '")) %>% filter(q.value < ', threshold, ') %>% arrange(q.value)')
  eval(parse(text = script))
  return(res)
}





### Wrapper function for plotting data with Gviz package --------
#' plotDetails
#' 
#' @description This function takes the output of gene name, methylation data, and sample group,
#' and plot methylation value on genomic position.
#'
#' @importFrom lattice
#'   bwplot
#' @param identifier TargetID of DNA methylation data.
#' @param methylation matrix of DNA methylation. Row as TargetID, Column as Sample. Only numeric (DNA methylation value)
#' @param sgroups vector of sample group
#' @return plot bwplot
#' @export

plotDetails <- function(identifier, methylation, sgroups, fill=NULL, ...) {
  dt <- data.frame(methylation = methylation[identifier, ], group = sgroups)
  print(lattice::bwplot(methylation~group, group = group, data = dt,
                        main = list(label = identifier, cex = 0.7),
                        scales = list(x = list(draw = TRUE, rot = 45)),
                        fill = fill), newpage = FALSE, prefix = "plot")
}


#' makeGvizObj
#' 
#' @description This function generate Gviz object
#'
#' @importFrom Gviz
#'   AnnotationTrack
#'   IdeogramTrack
#'   GenomeAxisTrack
#'   GeneRegionTrack
#'   DataTrack
#' @importFrom GenomicRanges
#'   GRanges
#' @importFrom IRanges
#'   IRanges
#' @importFrom GenomeInfoDb
#'   Seqinfo
#' @import tidyverse
#' @param gene vector of gene name, i.e. c('PGAP3', 'ERBB2')
#' @param data matrix of DNA methylation. Row as TargetID, Column as Sample. Only numeric (DNA methylation value)
#' @param genome the version of reference genome. (default: 'hg38')
#' @return list of Gviz object
#' @export

makeGvizObj <- function(gene, data, genome='hg38', ...) {
  if (!exists('methProbeAnnot')){
    data("methProbeAnnot")
  }
  if (!exists('geneAnnot')){
    data("geneAnnot")
  }
  probe <- tibble(genesUniq = gene) %>% inner_join(methProbeAnnot, .) %>% unnest(data) %>%
    distinct(CpG_chrm, CpG_beg, CpG_end, probe_strand, TargetID)
  gene <- geneAnnot %>% inner_join(tibble(symbol = gene)) %>% unnest(data) %>% distinct()
  chr <- probe %>% dplyr::distinct(CpG_chrm) %>% dplyr::pull(1)
  d <- probe %>% distinct(TargetID) %>% pull(TargetID) %>% base::intersect(rownames(data)) %>%
    base::unique() %>% data[., ]
  probe <- rownames(d) %>% tibble(TargetID = .) %>% inner_join(probe)
  probeViz <- Gviz::AnnotationTrack(range = GenomicRanges::GRanges(seqnames = probe$CpG_chrm,
                                                                   ranges = IRanges::IRanges(probe$CpG_beg,
                                                                                             probe$CpG_end,
                                                                                             names = probe$TargetID),
                                                                   seqinfo = GenomeInfoDb::Seqinfo(seqnames = chr,
                                                                                                   genome = genome)),
                                    name = 'Illumina Probes', stacking = 'dense', background.title = 'black')
  ideoViz <- Gviz::IdeogramTrack(chr, genome = genome)
  scaleViz <- Gviz::GenomeAxisTrack()
  geneViz <- Gviz::GeneRegionTrack(gene, genome = genome,
                                   chromosome = as.vector(gene$seqid[1]),
                                   name = 'Gene Model', transcriptAnnotation = 'symbol',
                                   background.title = 'brown')
  dTrack <- probe %>% distinct(CpG_beg, CpG_end, CpG_chrm) %$%
    Gviz::DataTrack(start = .$CpG_beg, end = .$CpG_end, chromosome = .$CpG_chrm, genome = genome, data = t(d))
  st <- dTrack@range %>% start() %>% min()
  ed <- dTrack@range %>% end() %>% max()
  wd <- (ed - st + 1) %/% 2
  return(list(AnnotationTrack = probeViz,
              IdeogramTrack = ideoViz,
              GenomeAxisTrack = scaleViz,
              GeneRegionTrack = geneViz,
              DataTrack = dTrack,
              methylation = as.matrix(d),
              from = st - wd, to = ed + wd + 1,
              genome = genome, chr = chr,
              gene = gene,
              probe = probe))
}



#' plotGvizObj
#' 
#' @description This function plots Gviz object
#'
#' @importFrom Gviz
#'   AnnotationTrack
#'   IdeogramTrack
#'   GenomeAxisTrack
#'   GeneRegionTrack
#'   DataTrack
#' @importFrom GenomicRanges
#'   GRanges
#' @importFrom IRanges
#'   IRanges
#' @importFrom GenomeInfoDb
#'   Seqinfo
#' @import tidyverse
#' @param gvizObj Gviz object
#' @param type vector for types of plotting by Gviz::plotTacks()
#' @param filename optional, filename for saving plot (default is NULL).
#'   If you want to save plot, also 'dir' must be passed this function.
#' @param dir output directory path (default is NULL).
#'   If you want to save plot, also 'filename' must be passed this function.
#' @param pfx optional, prefix for saving plot filenames (default is NULL)
#' @param ... arguments for Gviz::plotTracks()
#' @return plot Gviz track
#' @export

plotGvizObj <- function(gvizObj, type = c('a', 'p'), filename=NULL, dir=NULL, pfx = NULL, ...){
  if(min(!sapply(list(filename, dir), is.null))){
    fn <- paste0(dir, '/') %>% str_replace(., '\\/\\/$', '/') %>%
      paste0(., filename, '_', pfx, '.pdf') %>%
      str_replace(., '_\\.', '.')
    pdf(file = fn)
  }
  Gviz::plotTracks(list(gvizObj$IdeogramTrack,
                        gvizObj$GenomeAxisTrack,
                        gvizObj$AnnotationTrack,
                        gvizObj$GeneRegionTrack,
                        gvizObj$DataTrack),
                   type = type, ...)
  if(exists('fn')){
    dev.off()
  }
}


#' plotHeatmap
#' 
#' @description This function plots Heatmap of B matrix with or without sample annotation
#'
#' @import ComplexHeatmap
#' @import tidyverse
#' @param B B matrix obtained by 'getNewDataB' function
#' @param k number of cluster splitting samples (default is 2)
#' @param k.lv number of cluster splitting for LVs (default is 4)
#' @param annot.df optional. data frame of sample annotation. rowid is sample name.
#' @param col optional. list of annotaiton color. heatmap.col.RData as a test file.
#' @return plot Heatmap and cluster table
#' @export

plotHeatmap <- function(B, k = 2, k.lv = 4, annot.df=NULL, col=NULL){
  require(ComplexHeatmap)

  if (!is.null(annot.df)){
    if (!is.null(col)){
      ha <- HeatmapAnnotation(df = annot.df, col = col, annotation_name_side = 'left')
    } else {
      ha <- HeatmapAnnotation(df = annot.df, annotation_name_side = 'left')
    }
    B <- annot.df %>% rownames() %>% B[, .]
    h <- Heatmap(B, name ='value', top_annotation = ha,
                 column_names_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 5.5),
                 column_split = k, row_km = k.lv, row_km_repeats = 100,
                 clustering_method_columns = 'ward.D2',
                 clustering_distance_columns = 'euclidean',
                 show_column_names = FALSE, show_row_names = FALSE)
  } else {
    h <- Heatmap(B, name ='value',
                 column_names_gp = gpar(fontsize = 8),
                 row_names_gp = gpar(fontsize = 5.5),
                 column_split = k, row_km = k.lv, row_km_repeats = 100,
                 clustering_method_columns = 'ward.D2',
                 clustering_distance_columns = 'euclidean',
                 show_column_names = FALSE, show_row_names = FALSE)
  }

  h.p <- ComplexHeatmap::draw(h,
                              column_title = 'Samples', column_title_side = 'bottom',
                              row_title = 'LVs', row_title_side = 'right', merge_legend = TRUE)

  cluster.sample <- h.p %>% ComplexHeatmap::column_order() %>% imap(~ tibble(cluster.sample = .y, rowid.sample = .x)) %>%
    purrr::reduce(bind_rows) %>% arrange(rowid.sample) %>%
    mutate(AccessionNo = colnames(B))
  cluster.LV <- h.p %>% ComplexHeatmap::row_order() %>% imap(~ tibble(cluster.LV = .y, LV = .x)) %>%
    purrr::reduce(bind_rows) %>% arrange(LV) %>%
    bind_cols(B) %>% pivot_longer(cols = c(-1:-2), names_to = 'AccessionNo') %>%
    inner_join(cluster.sample) %>% group_by(LV, cluster.sample) %>%
    summarise(value = list(value)) %>% ungroup() %>%
    pivot_wider(names_from = cluster.sample, values_from = value)
  return(list(h = h.p, cluster.sample = cluster.sample, cluster.LV = cluster.LV))
}


#' plotSurvival
#' 
#' @description This function plots Survival plot
#'
#' @import survival
#' @import tidyverse
#' @param cl data for survival analysis. Data must following columns.
#' @param cluster numeric. Cluster number to which the sample belongs.
#' @param Time numeric. Time to event.
#' @param Event numeric. Event value for survival analysis.
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @return plot survival plot
#' @export

plotSurvival <- function(cl,
                         xlab='Time Since Surgery (year)',
                         ylab='Relapce-Free Survival (probability)'){
  library(survival)
  fit <- survfit(Surv(Time, Event) ~ cluster, data = cl)
  surv <- survminer::ggsurvplot(fit, risk.table = TRUE,
                                tables.theme = survminer::theme_cleantable(),
                                tables.col = 'strata', pval = TRUE, pval.method = TRUE)
  surv$plot <- surv$plot + ylab(ylab) + xlab(xlab)
  surv$table <- surv$table + theme(axis.ticks.y = element_blank(),
                                   axis.text.y = element_blank())
  surv
}



#' plotBoxplot
#' 
#' @description This function plots Boxplot
#'
#' @import tidyverse
#' @param B B matrix obtained by 'getNewDataB' function
#' @param LV LV number of to be drawn
#' @param cl tibble data table. category of sample. "AccessionNo", and "cluster.sample" column needed.
#' @return plot box plot
#' @export

plotBoxplot <- function(B, LV, cl){
  d <- B[LV, ] %>% as_tibble(rownames = 'AccessionNo') %>%
    inner_join(cl) %>% mutate(cluster.sample = as.factor(cluster.sample))
  d %>% ggplot(aes(y = value, fill = cluster.sample)) + geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          aspect.ratio = 2) +
    ylab('Score') + ggtitle(paste0('LV ', LV))
}



#' plotBoxplot.gviz
#' 
#' @description This function plots detail Track and GeneRegion Track
#'
#' @import tidyverse
#' @param gvizObj Gviz object
#' @param sgroups vector of sample group
#' @param target.id vector of TargetID what you want to plot
#' @param sizes vector of sizes for detail Track and GeneRegion Track (default is c(5, 1))
#' @param fill: vector of fill color for bwplot (default is NULL)
#' @param filename optional, filename for saving plot (default is NULL).
#'   If you want to save plot, also 'dir' must be passed this function.
#' @param dir output directory path (default is NULL).
#'   If you want to save plot, also 'filename' must be passed this function.
#' @param pfx optional, prefix for saving plot filenames (default is NULL)
#' @param width width of plot, numeric
#' @param height height of plot, numeric
#' @param ... arguments for lattice::bwplot()#'
#' @return plot box plot with Gviz
#' @export

plotBoxplot.gviz <- function(gvizObj, sgroups, target.id, sizes=c(5, 1), fill=NULL,
                             filename=NULL, dir=NULL, pfx = NULL,
                             width = 10, height = 12, ...){
  if(min(!sapply(list(filename, dir), is.null))){
    fn <- paste0(dir, '/') %>% str_replace(., '\\/\\/$', '/') %>%
      paste0(., filename, '_', pfx, '.pdf') %>%
      str_replace(., '_\\.', '.')
    pdf(file = fn, width = width, height = height)
  }
  if(is.null(fill)){
    colPal <- colorRampPalette(RColorBrewer::brewer.pal(name = 'Set1', n = 9))
    fill <- unique(sgroups) %>% length() %>% colPal(.)
    names(fill) <- unique(sgroups)
  }
  id <- gvizObj$probe %>% rowid_to_column() %>%
    inner_join(tibble(TargetID = target.id)) %>% pull(rowid)
  deTrack <- gvizObj %$%
    Gviz::AnnotationTrack(range = DataTrack@range[id],
                          genome = genome,
                          chromosome = chr,
                          id = rownames(methylation[id, ]), fun = plotDetails,
                          detailsFunArgs = list(methylation = methylation[id, ],
                                                sgroups = sgroups, fill = fill))
  Gviz:::plotTracks(list(deTrack,
                         gvizObj$GeneRegionTrack),
                    sizes = sizes, ...)
  if(exists('fn')){
    dev.off()
  }
}

#### Initialize Load -----------------------------------------------------------
data(pca.ord)
data(pca.w)
data(methPLIER)