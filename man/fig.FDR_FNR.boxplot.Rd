% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VisualFunc.R
\name{fig.FDR_FNR.boxplot}
\alias{fig.FDR_FNR.boxplot}
\title{Boxplot of FDR and FNR for Simulated data}
\usage{
fig.FDR_FNR.boxplot(DEA.list, truth.DEA.list, title, subset = NULL)
}
\arguments{
\item{DEA.list}{A list with each element as a normalization method including the differential expression analysis (DEA) results from multiple simulated data}

\item{truth.DEA.list}{A list of DEA results of benchmark data for multiple simulated dataset}

\item{title}{The title of the figure}

\item{subset}{optional A vector of a subset of genes/markers for this analysis}
}
\value{

}
\description{
Function for the boxplots of FDR and FNR based on golden truth (typically the simulated benchmark data set \code{\link{simulated.data()}})
and the differential expression analysis (typically the results from the simulated normalized test data set \code{\link{simulated.data()}}).
}
\examples{
simulated <- simulated.data(proportion = c(0.15, 0.25),  median = c(2, 4), numsets = 10)
norm.methods <- c("norm.TMM", "norm.TC", "norm.UQ", "norm.med")
test.DEA.list <- list()
for (i in norm.methods){
temp = lapply(1:10, function(x) pip.norm.DE(raw = simulated[[x]]$simulated_test,
                                            groups = simulated[[x]]$simulated_group, norm.method = i))
test.DEA.list <- append(test.DEA.list, list(temp))
}
names(test.DEA.list) = norm.methods
benchmark.DEA.list <- lapply(1:10, function(x) DE.voom(RC = simulated[[x]]$simulated_benchmark,
                                                       groups = simulated[[x]]$simulated_group, P = 0.01))
fig.FDR_FNR.boxplot(DEA.list = test.DEA.list, truth.DEA.list = benchmark.DEA.list, title = "Example of FDR and FNR boxplot")
}
