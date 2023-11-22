# Use rocker/tidyverse as the base image
FROM rocker/tidyverse

# Update the system and install required dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    unzip \
    libbz2-dev \
    libglpk-dev

# Install R packages
RUN R -e "install.packages(c('devtools', 'BiocManager'), repos='http://cran.rstudio.com/')"

# Execute installation using devtools
RUN R -e "remotes::install_github('hamamoto-lab/methPLIER')"

# Bind current directory
VOLUME ["/app"]
WORKDIR /app