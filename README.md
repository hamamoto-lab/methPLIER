# methPLIER
R package for DNA methylation data analysis using NMF with sparse regularization

## Library dependencies

## Installation

### Option 1: Install local
```R
# Install remotes from CRAN
install.packages("remotes")

# Pre-requiring package install
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("enrichplot", type = 'source')
BiocManager::install("DOSE", type = 'source')
BiocManager::install("ggtree", type = 'source')
BiocManager::install("Gviz", type = 'source')
remotes::install_github("wgmao/PLIER")

# Install methPLIER
remotes::install_github("hamamoto-lab/methPLIER")
```

### Option 2: Docker
```
# build docker image from Dockerfile
docker build -t methplier:0.9 .

# run docker image (i.e. password as "methPLIER")
docker run --rm -it -e PASSWORD=methPLIER -p 8787:8787 --name methplier methplier:0.9

# username: rstudio
# password: methPLIER (or your passward setting above)
```

## What the 'methPLIER' can do.
![schema](https://user-images.githubusercontent.com/7193590/172372421-db129640-486f-4f8c-a8f3-015fba7c58ab.png)


### Download figshare data
download dataset from figshare and put it in `data` directory

## Citation

## License
