
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
#'
#' @examples TODO
precision.seq <- function(norm.counts, adjust.factors, method.name,
                          DE.method="DE.voom", Pval=0.01) {

  ### Check Inputs

  if(!(DE.method %in% c("DE.voom", "DE.edgeR"))) {
    stop("Unknown method for differential expression (DE). Options are \"DE.voom\" and \"DE.edgeR\".")
  }

  if(missing(method.name)) {
    method.name <- "new"
  }

  ### Data Normalization
  ## Convert normalized data to internal format
  cat("Data normalization\n")

  if(missing(adjust.factors)){
    has.adjust <- F
    cat("No adjustment factors\n")
  } else {
    has.adjust <- T
    cat("Use adjustment factors\n")
  }

  if(has.adjust) {
    norm.new <- list(dat.normed=norm.counts, adjust.factor=adjust.factors)
    cat("Use adjustment factors\n")
  } else {
    norm.new <- list(dat.normed=norm.counts)
    cat("No adjustment factors\n")
  }

  ## Normalize test data using internal normalization functions
  test.norm <- pip.norm(raw=data.test, groups=data.group, norm.method = "all")
  # Append new (externally) normalized data
  test.norm <- append(list(norm.new), test.norm)
  names(test.norm)[1] <- method.name
  # Append un-normalized data
  test.norm <- append(test.norm, list(noNorm=list(dat.normed=data.test)))
  # names of all tested normalization methods
  norm.names <- names(test.norm)


  ### Differential Expression Analysis

  cat("\n\nDEA\n")

  # Wrapper for DEA function
  DEA <- function(RC, normalized = TRUE, adjust = NULL){
    FUN <- match.fun(DE.method)
    data.DE <- FUN(RC = RC, groups = data.group, Pval = Pval, normalized = normalized, adjust = adjust)
    return(data.DE)
  }

  ## DE for new external normalization
  if(has.adjust) {
    cat("DEA of new norm using adjust factors\n")
    new.DE <- DEA(RC=data.test, normalized=FALSE, adjust=norm.new$adjust.factor)
  } else {
    cat("DEA of new norm (no adjust)\n")
    new.DE <- DEA(RC=norm.new$dat.normed)
  }

  ## DE for internal normalization
  cat("DEA for provided norms\n")
  test.DE <- list(
    TMM = DEA(RC=test.norm$TMM$dat.normed),
    TC = DEA(RC=test.norm$TC$dat.normed),
    UQ = DEA(RC=test.norm$UQ$dat.normed),
    med = DEA(RC=test.norm$med$dat.normed),
    DESeq = DEA(RC=test.norm$DESeq$dat.normed),
    PoissonSeq = DEA(RC=test.norm$PoissonSeq$dat.normed),
    QN = DEA(RC=test.norm$QN$dat.normed),
    RUVg = DEA(RC=data.test, normalized=FALSE, adjust=test.norm$RUVg$adjust.factor),
    RUVs = DEA(RC=data.test, normalized=FALSE, adjust=test.norm$RUVs$adjust.factor),
    RUVr = DEA(RC=data.test, normalized=FALSE, adjust=test.norm$RUVr$adjust.factor),
    SVA = DEA(RC=data.test, normalized=FALSE, adjust=test.norm$SVA$adjust.factor),
    noNorm = DEA(RC=data.test)
  )
  # Add new method to data.DE list
  test.DE <- append(list(new.DE), test.DE)
  names(test.DE)[1] <- method.name

  ## DE for benchmark data
  cat("DE for benchmark data\n")
  bench.DE <- DEA(RC=data.benchmark)

  ## Compute DE statistics
  cat("Compute DE statistics\n")
  test.DE.stats <- list()
  for(DE.res in test.DE) {
    test.DE.stats <-
      append(test.DE.stats,
             list(DE.statistics(markers=rownames(data.test),
                                id.list=DE.res$id.list,
                                truth=bench.DE$id.list,
                                selected.marker=NULL)))
    # TODO: Add option for selecting a subset of markers for DEA analysis.
  }
  names(test.DE.stats) <- names(test.DE)


  ### Generate Figures

  # Volcano plots
  cat("Volcano plots\n")
  p.volcano <- vector("list", length(test.DE.stats))
  for(i in 1:length(test.DE.stats)) {
    p.volcano[[i]] <- fig.volcano(test.DE[[i]], names(test.DE)[i])
  }
  names(p.volcano) <- names(test.DE.stats)
  p.volcano <- append(p.volcano, list(benchmark=fig.volcano(bench.DE, "benchmark")))

  # RLE plots
  cat("RLE plots\n")
  p.RLE <- vector("list", length(test.norm))
  for(i in 1:length(test.norm)) {
    p.RLE[[i]] <- fig.RLE(test.norm[[i]]$dat.normed, data.group,
                          paste0("RLE for ", names(test.DE)[i]))
  }
  names(p.RLE) <- names(test.norm)
  p.RLE <- append(p.RLE, list(benchmark=fig.RLE(data.benchmark, data.group,
                                                "RLE for benchmark")))

  # CAT plot
  cat("CAT plot\n")
  ## CAT plot function copied from DANA project
  # Comcordance at the top plot
  #
  # Function for generating concordance at the top plot, which compares the
  # concordance of the p-values for given differential expression analyzes
  # to an assumed truth (a benchmark).
  #
  # @param DEA result for the data under study
  # @param truth.DEA DEA result for the assumed truth (the gold standard)
  # @param title Plot title
  # @param maxrank Optionally specify the maximum size of top-ranked items
  # that you want to plot.
  # @param subset vector of a subset of genes/markers for this analysis
  #
  # @return a list of values about TPR, FPR, FDR, FNR

  plot.CAT <- function(DEA, truth.DEA, title="", maxrank=100, subset=NULL){
    # Subset DEAs to a given set of miRNAs
    if (!is.null(subset)) {
      DEA <- DEA[subset, ]
      truth.DEA <- truth.DEA[subset, ]
    }
    # Reduce Data to named p.values
    truth <- truth.DEA$p.val
    names(truth) <- truth.DEA$id.list
    compare <- list()
    for(i in 1:length(DEA)) {
      compare <- append(compare, list(DEA[[i]]$p.val))
      names(compare[[i]]) <- DEA[[i]]$id.list
    }
    names(compare) <- names(DEA)

    catplots <- list()
    for(i in 1:length(compare)) {
      catplots <- append(catplots,
                         list(CATplot(compare[[i]], truth,
                                      maxrank = maxrank, make.plot=FALSE)))
    }
    names(catplots) <- names(compare)

    data_cat <- data.frame(x = numeric(), y = numeric(), curve = character())
    for(i in 1:length(catplots)){
      y = catplots[[i]][,2]
      x = 1:length(y)
      data_cat = rbind(data_cat, data.frame(x, y, curve = names(catplots)[i]))
    }
    data_cat$curve = factor(data_cat$curve)

    return(ggplot(data_cat, aes(x, y, color = curve, linetype = curve)) +
             theme(legend.title = element_blank()) +
             geom_line(size=.75) +
             ylab("Rate of Agreement with Benchmark") +
             xlab("Significance Rank") +
             theme(legend.title=element_blank()) +
             ggtitle(title) +
             ylim(c(0,1)) +
             theme_bw())
  }
  p.CAT <- plot.CAT(test.DE, bench.DE, title="Corncordance at the top")

  # Scatterplot (FNR-FDR)
  test.DEA.res <- data.frame(method=names(test.DE.stats),
                             FDR=sapply(test.DE.stats, function(x) x$FDR),
                             FNR=sapply(test.DE.stats, function(x) x$FNR))
  p.FNR_FDR <- ggplot(test.DEA.res, aes(x=FNR, y=FDR, label=method)) +
    theme_bw() +
    geom_point(alpha=1, color = "red") +
    xlab("FNR") + # True positive rate(1-FNR)
    ylab("FDR") + # Positive predictive value (1-FDR)
    geom_text_repel(aes(label = method)) +
    scale_x_continuous(labels = scales::percent, limits=c(0,1),
                       breaks = scales::pretty_breaks(n = 4)) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1))

  # Dendrogram
  cat("Dendrogram plot\n")
  genes <- rownames(data.test)
  pval_frame <- data.frame(sapply(test.DE, function(x) x$p.val[genes]))
  hc <- hclust(dist(t(-log10(pval_frame))), method = "ward.D")
  p.dendro <- ggdendrogram(hc, rotate = FALSE, size = 2)
  rm(genes, pval_frame, hc)

  # Venn Diagrams
  cat("Venn plots\n")
  p.venn <- vector("list", length(test.DE))
  bench.sig <- (bench.DE$p.val < Pval)
  test.sig <- lapply(test.DE, function (x) (x$p.val[names(bench.sig)] < Pval))
  for(i in 1:length(test.sig)) {
    p.venn[[i]] <- as.ggplot(function() vennDiagram(
      vennCounts(cbind(bench.sig, test.sig[[i]])),
      names = c("Benchmark", names(test.sig)[i]),
      cex = 1.5, counts.col = rainbow(1)))
  }
  names(p.venn) <- names(test.sig)
  # TODO -> Shows empty white plot in output. Why?



  return(list(data.norm=test.norm,
              data.DE=test.DE,
              benchmark.DE=bench.DE,
              DE.stats=test.DE.stats,
              p.volcano=p.volcano,
              p.RLE=p.RLE,
              p.CAT=p.CAT,
              p.FNR_FDR=p.FNR_FDR,
              p.dendrogram=p.dendro,
              p.venn=p.venn))
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
