#' Full pipeline function
#'
#' This internal function performs the full precision.seq pipeline.
#' It is used by the functions \code{\link{precision.seq}} and
#' \code{\link{pip.simulated.data}}.
#' @keywords internal
full.pipeline <- function(bench.data, test.data, groups.data,
                          norm.counts, adjust.factors,
                          method.name, DE.method, Pval, generate.plots) {
  ### Check Inputs

  if(any(dim(bench.data) != dim(test.data))) {
    stop("Dimensions of benchmark and test data must agree.")
  }

  if(!(DE.method %in% c("DE.voom", "DE.edgeR"))) {
    stop("Unknown method for differential expression (DE). Options are \"DE.voom\" and \"DE.edgeR\".")
  }

  if(missing(method.name)) {
    method.name <- "new"
  }

  ### Data Normalization
  ## Convert normalized data to internal format

  if(missing(adjust.factors)){
    has.adjust <- F
  } else {
    has.adjust <- T
  }

  if(has.adjust) {
    norm.new <- list(dat.normed=norm.counts, adjust.factor=adjust.factors)
  } else {
    norm.new <- list(dat.normed=norm.counts)
  }

  ## Normalize test data using internal normalization functions
  test.norm <- pip.norm(raw=test.data, groups=groups.data, norm.method = "all")
  # Append new (externally) normalized data
  test.norm <- append(list(norm.new), test.norm)
  names(test.norm)[1] <- method.name
  # Append un-normalized data
  test.norm <- append(test.norm, list(noNorm=list(dat.normed=test.data)))
  # names of all tested normalization methods
  norm.names <- names(test.norm)


  ### Differential Expression Analysis

  # Wrapper for DEA function
  DEA <- function(RC, normalized = TRUE, adjust = NULL){
    FUN <- match.fun(DE.method)
    data.DE <- FUN(RC = RC, groups = groups.data, Pval = Pval, normalized = normalized, adjust = adjust)
    return(data.DE)
  }

  ## DE for new external normalization
  if(has.adjust) {
    new.DE <- DEA(RC=test.data, normalized=FALSE, adjust=norm.new$adjust.factor)
  } else {
    new.DE <- DEA(RC=norm.new$dat.normed)
  }

  ## DE for internal normalization
  test.DE <- list(
    TMM = DEA(RC=test.norm$TMM$dat.normed),
    TC = DEA(RC=test.norm$TC$dat.normed),
    UQ = DEA(RC=test.norm$UQ$dat.normed),
    med = DEA(RC=test.norm$med$dat.normed),
    DESeq = DEA(RC=test.norm$DESeq$dat.normed),
    PoissonSeq = DEA(RC=test.norm$PoissonSeq$dat.normed),
    QN = DEA(RC=test.norm$QN$dat.normed),
    RUVg = DEA(RC=test.data, normalized=FALSE, adjust=test.norm$RUVg$adjust.factor),
    RUVs = DEA(RC=test.data, normalized=FALSE, adjust=test.norm$RUVs$adjust.factor),
    RUVr = DEA(RC=test.data, normalized=FALSE, adjust=test.norm$RUVr$adjust.factor),
    SVA = DEA(RC=test.data, normalized=FALSE, adjust=test.norm$SVA$adjust.factor),
    noNorm = DEA(RC=test.data)
  )
  # Add new method to data.DE list
  test.DE <- append(list(new.DE), test.DE)
  names(test.DE)[1] <- method.name

  ## DE for benchmark data
  bench.DE <- DEA(RC=bench.data)

  ## Compute DE statistics
  test.DE.stats <- list()
  for(DE.res in test.DE) {
    test.DE.stats <-
      append(test.DE.stats,
             list(DE.statistics(markers=rownames(test.data),
                                id.list=DE.res$id.list,
                                truth=bench.DE$id.list,
                                selected.marker=NULL)))
    # TODO: Add option for selecting a subset of markers for DEA analysis.
  }
  names(test.DE.stats) <- names(test.DE)


  ### Generate Figures

  if(generate.plots) {  # Case: Generate additional plots

    # Volcano plots
    p.volcano <- vector("list", length(test.DE.stats))
    for(i in 1:length(test.DE.stats)) {
      p.volcano[[i]] <- fig.volcano(test.DE[[i]], title=names(test.DE)[i])
    }
    names(p.volcano) <- names(test.DE.stats)
    p.volcano <- append(p.volcano, list(benchmark=fig.volcano(bench.DE, title="benchmark")))

    # RLE plots
    p.RLE <- vector("list", length(test.norm))
    for(i in 1:length(test.norm)) {
      p.RLE[[i]] <- fig.RLE(test.norm[[i]]$dat.normed, groups.data,
                            paste0("RLE for ", names(test.DE)[i]))
    }
    names(p.RLE) <- names(test.norm)
    p.RLE <- append(p.RLE, list(benchmark=fig.RLE(bench.data, groups.data,
                                                  "RLE for benchmark")))

    # CAT plot
    p.CAT <- fig.CAT(test.DE, bench.DE, title="Corncordance at the top")

    # Scatterplot (FNR-FDR)
    p.FNR_FDR <- fig.FDR_FNR(test.DE, bench.DE, title="FNR vs. FDR")

    # Dendrogram
    p.dendro <- fig.dendrogram(lapply(test.DE, function(x) x$p.val), title = "")

    # Venn Diagrams
    p.venn <- vector("list", length(test.DE))
    for(i in 1:length(test.DE)) {
      p.venn[[i]] <- fig.venn(bench.DE$p.val, test.DE[[i]]$p.val, Pvalue = Pval, title = paste0("Venn diagram for ", names(test.DE)[i]))
    }
    names(p.venn) <- names(test.norm)

    # List of return values
    res <- list(data.norm=test.norm,
                data.DE=test.DE,
                benchmark.DE=bench.DE,
                DE.stats=test.DE.stats,
                p.volcano=p.volcano,
                p.RLE=p.RLE,
                p.CAT=p.CAT,
                p.FNR_FDR=p.FNR_FDR,
                p.dendrogram=p.dendro,
                p.venn=p.venn)
  } else {
    # Case: Dont generate plots
    # List of return values
    res <- list(data.norm=test.norm,
                data.DE=test.DE,
                benchmark.DE=bench.DE,
                DE.stats=test.DE.stats)
  }


  return(res)
}


#' Full normalization assessment for given normalized test data
#'
#' @param norm.counts Normalized counts of the test data for a
#' user-provided normalization method.
#' @param adjust.factors optional Adjustment factors for the normalization of the
#' test data for a user-provided normalization method.
#' @param method.name optional Name of the normalization method. Used for
#' naming of data frames and plot legends.
#' @param DE.method Method for computing differential expression statuses.
#' Available: "DE.voom" and "DE.edgeR".
#' @param Pval P-value cutoff for differential expression.
#'
#' @return List of assessment statistics and plots (ggplot objects):
#' \describe{
#'   \item{data.norm}{List of Normalized test data containing the test data
#'   for each package-supplied normalization methods and the normalized test data
#'   provided by the user.}
#'   \item{data.DE}{List of results of the differential expression analysis
#'   using the differential expression method \code{DE.method} for each
#'   normalization method using the test data.}
#'   \item{benchmark.DE}{Results of the differential expression analysis using
#'   the method \code{DE.method} for the benchmark data.}
#'   \item{DE.stats}{List of statistics (true positive rate, false positive
#'   rate, false discovery rate, and false negative rate) based on the
#'   comparison of differential expression statuses between the normalized test
#'   data and the unnormalized gold standard (the benchmark data).}
#'   \item{p.volcano}{List of volcano plots from the differential expression
#'   analysis for each normalization method.}
#'   \item{p.RLE}{List of Relative Log Expression plots for the normalized
#'   and un-normalized test data and the un-normalized benchmark data.}
#'   \item{p.CAT}{Concordance At The Top plot containing all normalization
#'   methods.}
#'   \item{p.FNR_FDR}{Plot of the False Discovery Rate over the False Negative
#'   Rate of the concordance of the differential expression statuses of the
#'   normalized test data with the ones of the gold standard (the benchmark)
#'   for each normalization method under study.}
#'   \item{p.dendrogram}{Dendrogram from clustering p-values from the
#'   differential expression analysis for each normalization method.}
#'   \item{p.venn}{List of venn diagrams of the correspondance of differential
#'   expression statuses between the normalized test data and the gold standard
#'   (the benchmark).}
#' }
#' @export
precision.seq <- function(norm.counts, adjust.factors, method.name,
                          DE.method="DE.voom", Pval=0.01) {

  ### Check Inputs

  if(missing(method.name)) {
    method.name <- "new"
  }

  # Compute full precision seq pipeline using the full.pipeline function
  if(missing(adjust.factors)) {
    res <- full.pipeline(bench.data=data.benchmark,
                         test.data=data.test,
                         groups.data = data.group,
                         norm.counts=norm.counts,
                         method.name=method.name,
                         DE.method=DE.method,
                         Pval=Pval,
                         generate.plots=TRUE)
  } else {
    res <- full.pipeline(bench.data=data.benchmark,
                         test.data=data.test,
                         groups.data = data.group,
                         norm.counts=norm.counts,
                         adjust.factors=adjust.factors,
                         method.name=method.name,
                         DE.method=DE.method,
                         Pval=Pval,
                         generate.plots=TRUE)
  }

  return(res)
}




#' Normalization for RNASeq Data
#'
#' @param raw raw count data in the format of data frame or matrix, with columns
#' for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample
#' (only 2 groups allowed).
#' @param norm.method the method for normalization selected from
#'  \code{all}, \code{norm.none}
#' \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med},
#' \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA},
#' \code{norm.RUVr}, \code{norm.RUVg}, and \code{norm.RUVs}.
#' If \code{all} is selected, the function applies all provided normalization
#' methods to the data and returns a list normalized counts (and scaling factors
#' or adjusting factors depending on the method).
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#'
#' @return list, containing normalized dataset, and scaling factors or adjusting factors.
#' If method="all" the list contains the normalized data set and
#' scaling/adjusting factors for all provided methods.
#'
#' @export
#'
#' @examples
#' res <- pip.norm(data.test, data.group, "norm.TMM")
pip.norm <- function(raw, groups, norm.method="all",
                     QN_filter = FALSE) {
  if(norm.method == "all") {
    data.norm <- list(
      TMM = norm.TMM(raw, groups),
      TC = norm.TC(raw, groups),
      UQ = norm.UQ(raw, groups),
      med = norm.med(raw, groups),
      DESeq = norm.DESeq(raw, groups),
      PoissonSeq = norm.PoissonSeq(raw),
      RUVg = norm.RUVg(raw, groups),
      RUVs = norm.RUVs(raw, groups),
      RUVr = norm.RUVr(raw, groups),
      SVA = norm.SVA(raw, groups),
      QN = norm.QN(raw, filter=QN_filter)
    )
  } else if (norm.method %in% c("norm.TMM", "norm.TC", "norm.UQ", "norm.med",
                                "norm.DESeq", "norm.RUVg",
                                "norm.RUVs", "norm.RUVr", "norm.SVA")) {
    FUN <- match.fun(norm.method)
    data.norm <- FUN(raw, groups)
  } else if (norm.method == "norm.PoissonSeq") {
    data.norm <- norm.PoissonSeq(raw)
  } else if (norm.method == "norm.QN") {
    data.norm <- norm.QN(raw, filter = QN_filter)
  } else {
    stop("Error: Unknown value for parameter norm.method in function pip.norm.")
  }

  return(data.norm)
}


#' Pipeline of Differential Expression Analysis for RNASeq Data
#'
#' @param raw raw count data in the format of data frame or matrix, with columns
#' for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample
#' (only 2 groups allowed).
#' @param norm.method the method for normalization selected from \code{all}
#' \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ},
#' \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN},
#' \code{norm.SVA}, \code{norm.RUVr}, \code{norm.RUVg}, and \code{norm.RUVs}.
#' The function applies all available normalization methods if \code{all}
#' is selected.
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#' @param DE.method the method for differential expression analysis from
#' \code{DE.voom} and \code{DE.edgeR}, default to be \code{DE.voom}.
#' @param Pval p-value for identifying DE genes, default to be 0.01
#'
#' @return list, containing \code{id.list} (names of DE genes), \code{p.val},
#' and \code{log2.FC} for a single normalization. If method="all", a list of
#' of DEAs is returned for the raw data normalized with each supported
#' normaliztion methods vontaining containing \code{id.list} (names of DE genes),
#' \code{p.val}, and \code{log2.FC} each.
#'
#' @export
#'
#' @examples
#' res <- pip.norm.DE(data.test, data.group, "norm.TMM")
pip.norm.DE <- function(raw, groups, norm.method, QN_filter = FALSE, DE.method = "DE.voom", Pval = 0.01) {
  DEA <- function(RC, groups, Pval, DE.method, normalized = TRUE, adjust = NULL){
    FUN <- match.fun(DE.method)
    data.DE <- FUN(RC = RC, groups = groups, Pval = Pval, normalized = normalized, adjust = adjust)
    return(data.DE)
  }

  if (norm.method == "all") {
    data.DE <- list(
      noNorm = DEA(RC=raw, groups=groups, Pval=Pval, DE.method=DE.method),
      TMM = DEA(RC=norm.TMM(raw, groups)$dat.normed, groups=groups, Pval=Pval, DE.method=DE.method),
      TC = DEA(RC=norm.TC(raw, groups)$dat.normed, groups=groups, Pval=Pval, DE.method=DE.method),
      UQ = DEA(RC=norm.UQ(raw, groups)$dat.normed, groups=groups, Pval=Pval, DE.method=DE.method),
      med = DEA(RC=norm.med(raw, groups)$dat.normed, groups=groups, Pval=Pval, DE.method=DE.method),
      DESeq = DEA(RC=norm.DESeq(raw, groups)$dat.normed, groups=groups, Pval=Pval, DE.method=DE.method),
      PoissonSeq = DEA(RC=norm.PoissonSeq(raw)$dat.normed, groups=groups, Pval=Pval, DE.method=DE.method),
      QN = DEA(RC=norm.QN(raw, QN_filter)$dat.normed, groups=groups, Pval=Pval, DE.method=DE.method),
      RUVg = DEA(RC=raw, groups=groups, Pval=Pval, DE.method = DE.method, normalized=FALSE, adjust=norm.RUVg(raw, group = groups)$adjust.factor),
      RUVs = DEA(RC=raw, groups=groups, Pval=Pval, DE.method = DE.method, normalized=FALSE, adjust=norm.RUVs(raw, group = groups)$adjust.factor),
      RUVr = DEA(RC=raw, groups=groups, Pval=Pval, DE.method = DE.method, normalized=FALSE, adjust=norm.RUVr(raw, group = groups)$adjust.factor),
      SVA = DEA(RC=raw, groups=groups, Pval=Pval, DE.method = DE.method, normalized=FALSE, adjust=norm.SVA(raw, group = groups)$adjust.factor)
    )
  }
  else if (norm.method == "norm.none") {
    data.DE = DEA(RC=raw, groups=groups, Pval, DE.method = DE.method)
  }
  else if (norm.method %in% c("norm.TMM", "norm.TC", "norm.UQ",
                              "norm.med", "norm.DESeq")) {
    FUN <- match.fun(norm.method)
    test.norm <- FUN(raw, groups)$dat.normed
    data.DE <- DEA(RC = test.norm, groups=groups, Pval, DE.method = DE.method)
  }
  else if (norm.method == "norm.PoissonSeq") {
    FUN <- match.fun(norm.method)
    test.norm <- norm.PoissonSeq(raw)$dat.normed
    data.DE <- DEA(RC = test.norm, groups=groups, Pval, DE.method = DE.method)
  }
  else if (norm.method == "norm.QN") {
    test.norm <- norm.QN(raw, filter=QN_filter)$dat.normed
    data.DE <- DEA(RC = test.norm, groups=groups, Pval, DE.method = DE.method)
  }
  else if (norm.method %in% c("norm.SVA", "norm.RUVg", "norm.RUVs", "norm.RUVr")) {
    FUN <- match.fun(norm.method)
    test.norm <- FUN(raw, group = groups)$adjust.factor
    data.DE <- DEA(RC=raw, groups=groups, Pval=Pval, DE.method = DE.method, normalized = FALSE, adjust = test.norm)
  } else {
    stop("Error: Unknown value for parameter norm.method in function pip.norm.DE")
  }

  return(data.DE = data.DE)
}

#' Full normalization assessment for simulated data
#'
#' Simulated data can be sampled using the function \code{\link{simulated.data}}.
#'
#' @param data list of paired data sets. Each pair must consist of
#' a benchmark data set "simulated_benchmark" and a test data set
#' "simulated_test".
#' @param groups Sample groups in the simulated test and benchmark data sets.
#' Must be the same across all pairs.
#' @param norm.counts Normalized counts of each simulated test data set for a
#' user-provided normalization method.
#' @param adjust.factors optional Adjustment factors for the normalization of the
#' test data for a user-provided normalization method.
#' @param method.name optional Name of the normalization method. Used for
#' naming of data frames and plot legends.
#' @param DE.method Method for computing differential expression statuses.
#' Available: "DE.voom" and "DE.edgeR".
#' @param Pval P-value cutoff for differential expression.
#'
#'
#' @return List of assessment statistics and plots (ggplot objects):
#' \describe{
#'   \item{DE.stats}{List of statistics (true positive rate, false positive
#'   rate, false discovery rate, and false negative rate) based on the
#'   comparison of differential expression statuses between the normalized test
#'   data and the unnormalized gold standard (the benchmark data).
#'   Contains results for each pair of simulated data sets.}
#'   \item{p.boxplot}{Boxplot for False Negative Rate (FNR) and False
#'   Discovery Rate (FDR) of the agreement of DE statuses between the simulated
#'   test and benchmark data sets.}
#' }
#'
#' @export
pip.simulated.data <- function(data, groups, norm.counts, adjust.factors, method.name,
                               DE.method="DE.voom", Pval=0.01) {
  ### Check Inputs
  if(!is.list(data)) {
    stop("Unknown data. Must be a list of pairs of data sets containing simulated_benchmark and simulated_test.")
  }
  if(names(data[[1]]) != c("simulated_benchmark", "simulated_test")) {
    stop("Unknown data. Must be a list of pairs of data sets containing simulated_benchmark and simulated_test.")
  }
  num.datasets <- length(data)

  if(num.datasets != length(norm.counts)) {
    stop("Number of normalized data sets (", length(norm.counts), ") does not correspond to the number of un-nomrlized data sets (", num.datasets, ")")
  }

  if(!(DE.method %in% c("DE.voom", "DE.edgeR"))) {
    stop("Unknown method for differential expression (DE). Options are \"DE.voom\" and \"DE.edgeR\".")
  }

  if(missing(method.name)) {
    method.name <- "new"
  }

  ### Data Normalization
  ## Convert normalized data to internal format

  if(missing(adjust.factors)){
    has.adjust <- F
  } else {
    has.adjust <- T
    if(num.datasets != length(adjust.factors)) {
      stop("Number of sets of adjustment factors (", length(adjust.factors), ") does not correspond to the number of un-nomrlized data sets (", num.datasets, ")")
    }
  }

  DE.stats <- list()
  for(i in 1:num.datasets) {
    if(has.adjust) {
      # Case: using adjustment factors
      pipeline.res <- full.pipeline(bench.data=data[[i]]$simulated_benchmark,
                                    test.data=data[[i]]$simulated_test,
                                    groups.data = groups,
                                    norm.counts=norm.counts[[i]],
                                    adjust.factors=adjust.factors[[i]],
                                    method.name=method.name,
                                    DE.method=DE.method,
                                    Pval=Pval,
                                    generate.plots=FALSE)
    } else {
      # Case: normalized counts only (no adjustment factors)
      pipeline.res <- full.pipeline(bench.data=data[[i]]$simulated_benchmark,
                                    test.data=data[[i]]$simulated_test,
                                    groups.data = groups,
                                    norm.counts=norm.counts[[i]],
                                    method.name=method.name,
                                    DE.method=DE.method,
                                    Pval=Pval,
                                    generate.plots=FALSE)
    }
    DE.stats <- append(DE.stats, list(pipeline.res$DE.stats))
  }

  FNR <- FDR <- data.frame()
  for(element in DE.stats){
    FNR <- rbind(FNR, data.frame(lapply(element, function(x) x$FNR)))
    FDR <- rbind(FDR, data.frame(lapply(element, function(x) x$FDR)))
  }

  FNR_FDR_tab <- data.frame(
    value=c(unname(unlist(FDR)), unname(unlist(FNR))),
    method=rep(as.vector(sapply(colnames(FNR), function(x) rep(x,num.datasets))),2),
    stat=c(rep("FDR", length(unlist(FDR))), rep("FNR", length(unlist(FDR)))))


  p <- ggplot(data = FNR_FDR_tab , aes(x=method, y=value)) +
    geom_boxplot(aes(fill=stat))  +
    theme_bw() +
    scale_fill_manual(values = c("white", "grey")) +
    xlab("") + ylab("")



  return(list(DE.stats=DE.stats,
              p.boxplot=p))
}

