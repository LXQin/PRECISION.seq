# DepthNorm: An R Package for Comparsion of Depth Normalization Methods

We propose the well-prepared datasets and pack commonly used normalization methods/functions in this package. This package is based on R 4.0.2. The general workflow of this package is: 
* Input dataset in the format of data frame or matrix with with columns for samples and raws for genes;
* Indicate the group for each sample using a vector of characters;
* Ensure the normalization and DEA method;
* Obtain the output including the DE genes and p values for each gene.

Before installing this package, please ensure that the dependency packages installed. The R codes for installing these packages are:

```R
## from CRAN
CRAN.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}
CRAN.packages(c("PoissonSeq", "DescTools", "BiocManager", "readr"))

## from Bioconductor
Bioconductor.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg, dependencies = TRUE)
}
Bioconductor.packages(c("DESeq", "edgeR", "affy", "sva", "RUVSeq", "EDASeq", "limma", "preprocessCore", "ffpe", "Biobase", "vsn"))
```

The original R code and figures for the paper [*Statistical Assessment of Depth Normalization for Small RNA Sequencing*](https://pubmed.ncbi.nlm.nih.gov/32598180/) is inclueded in for reference.
