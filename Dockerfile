# Use rocker/tidyverse as the base image
FROM rocker/tidyverse

# Set environment variables
ENV DOWNLOAD_URL=https://figshare.com/ndownloader/articles/21938528

# Update the system and install required dependencies
RUN apt-get update && apt-get install -y \
    curl \
    unzip

# Download files from figshare
RUN curl -LO $DOWNLOAD_URL/files.zip && \
    unzip files.zip -d /data && \
    rm files.zip

# Install R packages
RUN R -e "install.packages(c('devtools', 'BiocManager'), repos='http://cran.rstudio.com/')"

# Execute installation using devtools
RUN R -e "devtools::install_github('hamamoto-lab/methPLIER')"

# Bind current directory
VOLUME ["/app"]
WORKDIR /app

# Run R when the container launches
CMD ["R"]