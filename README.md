# precision.seq: An R Package for Comparison of Depth Normalization Methods

We propose the well-prepared datasets and pack commonly used normalization methods/functions in this package. This package provides a pipeline for evaluting depth normalization methods, especially for small RNA sequencing. The general workflow of this package is: 
* Input dataset in the format of data frame or matrix with columns for samples and rows for genes;
* Indicate the group for each sample using a vector of characters;
* Select the normalization and DEA method for analysis;
* Obtain the output including the DE genes and p-values for each gene;
* Visualize the outputs for methods comparision.

The package can be installed in R.

```R
devtools::install_github("LXQin/precision.seq")
```

If the package cannot be installed successfully, please ensure that the dependency packages are installed. This package is based on R 4.0.2, and the R codes for installing the dependent packages are:

```R
## from CRAN
CRAN.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
}
CRAN.packages(c("PoissonSeq", "DescTools", "BiocManager", "readr", "magrittr", "ggplot2", "ggrepel", "ggdendro", "data.table", "tidyr", "dplyr", "ggplotify", "cluster"))

## from Bioconductor
Bioconductor.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg, dependencies = TRUE)
}
Bioconductor.packages(c("DESeq2", "edgeR", "affy", "sva", "RUVSeq", "EDASeq", "limma", "preprocessCore", "ffpe", "Biobase", "vsn"))
```

The original R code and figures for the paper [*Statistical Assessment of Depth Normalization for Small RNA Sequencing*](https://pubmed.ncbi.nlm.nih.gov/32598180/) are included in the article *Pipeline of the Paper* for reference.
