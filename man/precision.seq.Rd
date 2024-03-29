% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CombAnalysis.R
\name{precision.seq}
\alias{precision.seq}
\title{Full normalization assessment for given normalized test data}
\usage{
precision.seq(
  norm.counts,
  adjust.factors,
  method.name,
  DE.method = "DE.voom",
  Pval = 0.01
)
}
\arguments{
\item{norm.counts}{Normalized counts of the test data for a
user-provided normalization method.}

\item{adjust.factors}{optional Adjustment factors for the normalization of the
test data for a user-provided normalization method.}

\item{method.name}{optional Name of the normalization method. Used for
naming of data frames and plot legends.}

\item{DE.method}{Method for computing differential expression statuses.
Available: "DE.voom" and "DE.edgeR".}

\item{Pval}{P-value cutoff for differential expression.}
}
\value{
List of assessment statistics and plots (ggplot objects):
\describe{
  \item{data.norm}{List of Normalized test data containing the test data
  for each package-supplied normalization methods and the normalized test data
  provided by the user.}
  \item{data.DE}{List of results of the differential expression analysis
  using the differential expression method \code{DE.method} for each
  normalization method using the test data.}
  \item{benchmark.DE}{Results of the differential expression analysis using
  the method \code{DE.method} for the benchmark data.}
  \item{DE.stats}{List of statistics (true positive rate, false positive
  rate, false discovery rate, and false negative rate) based on the
  comparison of differential expression statuses between the normalized test
  data and the unnormalized gold standard (the benchmark data).}
  \item{p.volcano}{List of volcano plots from the differential expression
  analysis for each normalization method.}
  \item{p.RLE}{List of Relative Log Expression plots for the normalized
  and un-normalized test data and the un-normalized benchmark data.}
  \item{p.CAT}{Concordance At The Top plot containing all normalization
  methods.}
  \item{p.FNR_FDR}{Plot of the False Discovery Rate over the False Negative
  Rate of the concordance of the differential expression statuses of the
  normalized test data with the ones of the gold standard (the benchmark)
  for each normalization method under study.}
  \item{p.dendrogram}{Dendrogram from clustering p-values from the
  differential expression analysis for each normalization method.}
  \item{p.venn}{List of venn diagrams of the correspondance of differential
  expression statuses between the normalized test data and the gold standard
  (the benchmark).}
}
}
\description{
Full normalization assessment for given normalized test data
}
