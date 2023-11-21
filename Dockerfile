# Use rocker/tidyverse as the base image
FROM rocker/tidyverse

# Update the system and install required dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    unzip

# Install R packages
RUN R -e "install.packages(c('devtools', 'BiocManager'), repos='http://cran.rstudio.com/')"

# Pre-requiring package install
RUN R -e "BiocManager::install(c('enrichplot', 'DOSE', 'ggtree', 'Gviz'), type = 'source')"
RUN R -e "remotes::install_github('wgmao/PLIER')"

# Execute installation using devtools
RUN R -e "devtools::install_github('hamamoto-lab/methPLIER')"

# Bind current directory
VOLUME ["/app"]
WORKDIR /app