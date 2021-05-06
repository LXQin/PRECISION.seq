#' Full normalization assessment for given normalized test data
#'
#' @param normalized.test
#'
#' @return
#' @export
#'
#' @examples
precision.seq <- function(norm.counts, adjust.factors, method.name,
                          DE.method="DE.voom", Pval=0.01) {
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



  return(list(data.norm=test.norm,
              data.DE=test.DE,
              data.DE.stats=test.DE.stats,
              p.volcano=p.volcano))
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
pip.norm.DE <- function(raw, groups, norm.method,
                        QN_filter = FALSE,
                        DE.method = "DE.voom", Pval = 0.01) {
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


#' Statistics for DEA Results Based Golden Standards
#'
#' Computing the true positive rate, false positive rate, false discovery rate, and false negative rate based on given golden standards.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample (only 2 groups allowed).
#' @param truth vector of genes that are truly differential expressed
#' @param DE.method the method for differential expression analysis from \code{DE.voom} and \code{DE.edgeR}, default to be \code{DE.voom}.
#' @param norm.method the method for normalization selected from \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#' @param Pval p-value for identifying DE genes, default to be 0.01
#' @param selected.marker if given, provide vector of a subset of genes/markers for
#' this analysis. Leave \code{NULL} if all markers are considered for the analysis.
#' @param new.norm.list a vector of DE genes selected from the data normalized by the new method (the method which is not in the provided in normalization methods).
#'
#' @return a list of values about TPR, FPR, FDR, FNR
#' @export
#'
#' @examples
#' truthgene <- DE.voom(data.benchmark, data.group)$id.list
#' stat <- pip.statistics(data.test, data.group, truth=truthgene,
#'                        DE.method="DE.voom", norm.method="norm.TMM")
#'
pip.statistics <- function(raw, groups, truth, DE.method = "DE.voom",
                           norm.method, QN_filter = FALSE, Pval = 0.01,
                           selected.marker = NULL, new.norm.list = NULL) {
  norm.method.list <- c("norm.TMM", "norm.TC", "norm.UQ", "norm.med",
                        "norm.DESeq", "norm.RUVg", "norm.RUVs", "norm.RUVr",
                        "norm.SVA", "norm.PoissonSeq", "norm.QN", "norm.none")

  if(norm.method %in% norm.method.list) {
    voom.test = pip.norm.DE(raw = raw, groups = groups,
                          norm.method = norm.method, QN_filter = QN_filter,
                          DE.method = DE.method, Pval = Pval)
  }

  DEgenelist = if(norm.method %in% norm.method.list) {voom.test$id.list} else {new.norm.list}
  statistics <- DE.statistics(markers=rownames(raw),
                              id.list=DEgenelist,
                              truth=truth,
                              selected.marker=selected.marker)
  return(statistics)
}


#' Pipeline for visualization Using
#'
#' A pipeline for visualization to compare the normalization method with selected existing methods
#' @param norm.method the normalization method selected for comparsion from \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param QN_filter whether the filtering is performed if \code{norm.QN} selected
#' @param new.method.name the name of new normalization method
#' @param normalized.test the normalized test data by the new normalization method
#' @param adjust.factor the adjusting factors for the test data by the new normalization method
#' @param DE.method the method for differential expression analysis from \code{DE.voom} and \code{DE.edgeR}, default to be \code{DE.voom}
#' @param Pval_cutoff  p-value for identifying DE genes, default to be 0.01
#'
#' @return
#' @export
#'
#' @examples
#' vsn.norm <- vsn::justvsn(data.test)
#' vsn.norm <- ifelse(vsn.norm<0, 1, vsn.norm)
#' pip.fig(c("norm.none", "norm.TMM", "norm.SVA", "norm.TC", "norm.RUVr"),
#'  new.method.name = "vsn", normalized.test = vsn.norm)
pip.fig <- function(norm.method, QN_filter = FALSE,
                    new.method.name, normalized.test = NULL, adjust.factor = NULL,
                    DE.method = "DE.voom", Pval_cutoff = 0.01) {

  DEA <- function(RC, groups, Pval, DEmethod = DE.method, normalized = TRUE, adjust = NULL){
    FUN <- match.fun(DEmethod)
    data.DE <- FUN(RC = RC, groups = groups, Pval = Pval, normalized = normalized, adjust = adjust)
    return(data.DE)
  }

  if (!is.null(normalized.test) & !is.null(adjust.factor)) {
    stop("Only one of normalized.test and adjust.factor can be input")
  }

  if (is.null(normalized.test)){
    new.method.dea = DEA(RC = data.test, groups = data.group, Pval = Pval_cutoff, normalized = FALSE, adjust = adjust.factor)
  } else{
    new.method.dea = DEA(RC = normalized.test, groups = data.group, Pval = Pval_cutoff, normalized = TRUE, adjust = NULL)
  }

  ## fig.volcano
  volcano <- fig.volcano(new.method.dea, title = "Volcano Plot")

  ## fig.RLE
  if (is.null(normalized.test)) {
    warning("No relative log expression generated beacuse of no normalized test provided")
  } else {
    RLE <- fig.RLE(normalized.test, data.group, paste0("RLE for ", new.method.name))
  }

  ## fig.CAT
  CAT <- fig.CAT(MethodsCompare = norm.method, benchmark = data.benchmark, test = data.test,
                 group = data.group, DE.method = DE.method, QN_filter = QN_filter,
                 Pval = Pval_cutoff, MethodNew = new.method.name,
                 pvalues = new.method.dea$p.val)

  ## fig.FNR_FDR
  DEFUN <- match.fun(DE.method)
  truthgene <- DEFUN(data.benchmark, data.group)
  FDR_FNR <- fig.FDR_FNR(raw = data.test, groups = data.group, MethodsCompare = norm.method,
                         truth = truthgene$id.list, QN_filter = QN_filter,
                         DE.method = DE.method, Pval = Pval_cutoff,
                         MethodNew = new.method.name, new.norm.list = new.method.dea$id.list, selected.marker = NULL)

  ## fig.venn
  venn <- fig.venn(truthgene$p.val, new.method.dea$p.val, Pvalue = Pval_cutoff)

  ## fig.dendrogram
  dendrogram <- fig.dendrogram(MethodsCompare = norm.method, DE.method = DE.method, QN_filter = QN_filter,
                               Pval = Pval_cutoff, MethodNew = new.method.name, pvalues = new.method.dea$p.val,
                               Methods_visual = c(MethodsCompare, MethodNew),
                               test = data.test, group = data.group)

  return(list(volcano = volcano,
              RLE = RLE,
              CAT = CAT,
              FDR_FNR = FDR_FNR,
              venn = venn,
              dendrogram = dendrogram))
}

