# methPLIER
R package for DNA methylation data analysis using NMF with sparse regularization

## Library dependencies

## Installation

### Option 1: Install local
```R
# Install devtools from CRAN
install.packages("devtools")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("hamamoto-lab/methPLIER")
```

### Option 2: Docker
```
# build docker image from Dockerfile
docker image build -t methplier:0.9 .

# run docker image (i.e. password as "methPLIER")
docker run --rm -it -e PASSWORD=methPLIER -p 8787:8787 --name methplier:0.9
```

## What the 'methPLIER' can do.
![schema](https://user-images.githubusercontent.com/7193590/172372421-db129640-486f-4f8c-a8f3-015fba7c58ab.png)


### Download figshare data
download dataset from figshare and put it in `data` directory

## Citation

## License
