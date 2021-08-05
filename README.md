# precision.seq: An R Package for Performance Assessment of Depth Normalization Methods in MicroRNA Sequencing

We make available a new R package PRECISION.seq (PaiREd miCrorna analysIs of differential expresSION for sequencing) for assessing the performance of depth normalization methods based on differential expression status in microRNA sequencing. The package provides a pair of microRNA sequencing data sets for the same set of tumor samples, additional simulated pairs of data sets under various patterns of differential expression, and a collection of numerical and graphical tools for normalization assessment. Users can easily assess their own normalization method and compare its performance to nine popular methods already implemented in the package. PRECISION.seq enables an objective and systemic evaluation of depth normalization methods in microRNA sequencing using realistically distributed and robustly benchmarked data. 

The package can be installed in R. The full *package documentation* can be found [here](https://lxqin.github.io/PRECISION.seq/).

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

The original R code along with the figures for the paper [*Statistical Assessment of Depth Normalization for Small RNA Sequencing*](https://pubmed.ncbi.nlm.nih.gov/32598180/) is included in the article *Pipeline of the Paper*.
