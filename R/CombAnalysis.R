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
pip.norm <- function(raw, groups, norm.method,
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
