#' Normalization for RNASeq Data
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample (only 2 groups allowed).
#' @param norm.method the method for normalization selected from \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, \code{norm.RUVr}, \code{norm.RUVg}, and \code{norm.RUVs}.
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#'
#' @return list, containing normalized dataset, and scaling factors or adjusting factors.
#'
#' @export
#'
#' @examples
#' res <- pip.norm(data.test, data.group, "norm.TMM")
pip.norm <- function(raw, groups, norm.method,
                     QN_filter = FALSE) {

  if (norm.method %in% c("norm.TMM", "norm.TC", "norm.UQ", "norm.med",
                         "norm.DESeq", "norm.PoissonSeq", "norm.RUVg",
                         "norm.RUVs", "norm.RUVr", "norm.SVA")) {
    FUN = match.fun(norm.method)
    test.norm = FUN(raw, groups)
  }

  if (norm.method == "norm.QN") {
    test.norm = norm.QN(raw, groups, filter = filter)
  }

  return(normalized = test.norm)
}


#' Pipeline of Differential Expression Analysis for RNASeq Data
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample (only 2 groups allowed).
#' @param norm.method the method for normalization selected from \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, \code{norm.RUVr}, \code{norm.RUVg}, and \code{norm.RUVs}.
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#' @param DE.method the method for differential expression analysis from \code{DE.voom} and \code{DE.edgeR}, default to be \code{DE.voom}.
#' @param Pval p-value for identifying DE genes, default to be 0.01
#'
#' @return list, containing \code{id.list} (names of DE genes), \code{p.val}, and \code{log2.FC}.
#'
#' @export
#'
#' @examples
#' res <- pip.norm.DE(data.test, data.group, "norm.TMM")
pip.norm.DE <- function(raw, groups, norm.method,
                        QN_filter = FALSE,
                        DE.method = "DE.voom", Pval = 0.01) {
  DEA <- function(raw, DE.method, normalized = TRUE, adjust = NULL){
    FUN = match.fun(DE.method)
    test.DE = FUN(RC = raw, groups = groups, Pval = Pval, normalized = normalized, adjust = adjust)
    return(test.DE = test.DE)
  }

  if (norm.method == "norm.none") {
    test.DE = DEA(raw = raw, DE.method = DE.method)
  }

  if (norm.method %in% c("norm.TMM", "norm.TC", "norm.UQ", "norm.med",
                         "norm.DESeq", "norm.PoissonSeq")) {
    FUN = match.fun(norm.method)
    test.norm = FUN(raw, groups)$dat.normed
    test.DE = DEA(raw = test.norm, DE.method = DE.method)
  }

  if (norm.method == "norm.QN") {
    test.norm = norm.QN(raw, groups, filter = filter)$dat.normed
    test.DE = DEA(raw = test.norm, DE.method = DE.method)
  }

  if (norm.method %in% c("norm.SVA", "norm.RUVg", "norm.RUVs", "norm.RUVr")) {
    FUN = match.fun(norm.method)
    test.norm = FUN(raw, group = groups)$adjust.factor
    test.DE = DEA(raw = raw, DE.method = DE.method, normalized = FALSE, adjust = test.norm)
  }
  return(test.DE = test.DE)
}


#' Statistics for DEA Results Based Golden Standards
#'
#' Computing the true positive rate, false positive rate, false discovery rate, and false negative rate based on given golden standards.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample (only 2 groups allowed).
#' @param norm.method the method for normalization selected from \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#' @param DE.method the method for differential expression analysis from \code{DE.voom} and \code{DE.edgeR}, default to be \code{DE.voom}.
#' @param Pval p-value for identifying DE genes, default to be 0.01
#' @param truth vector of genes that are truly differential expressed
#' @param marker_selection whether selecting a subset of genes/markers for this analysis
#' @param selected_marker vector of a subset of genes/markers for this analysis
#'
#' @return a list of values about TPR, FPR, FDR, FNR
#' @export
#'
#' @examples
#' truthgene <- DE.voom(data.benchmark, data.group)$id.list
#' stat <- pip.statistics(data.test, data.group, "norm.TMM", truth = truthgene)
#'
pip.statistics <- function(raw, groups, norm.method,
                          QN_filter = FALSE,
                          DE.method = "DE.voom", Pval = 0.01, truth,
                          marker_selection = FALSE, selected_marker = NULL) {
  voom.test = pip.norm.DE(raw = raw, groups = groups,
                          norm.method = norm.method, QN_filter = QN_filter,
                          DE.method = DE.method, Pval = Pval)

  stat.DE = matrix(, nrow = nrow(raw), ncol = 2)
  rownames(stat.DE) = rownames(raw)
  colnames(stat.DE) = c("Benchmark", "NormTest")
  stat.DE[truth, 1] = "DE"
  stat.DE[voom.test$id.list, 2] = "DE"
  stat.DE[is.na(stat.DE)] = "non-DE"
  if (marker_selection == TRUE) {stat.DE = stat.DE[selected_marker,]}
  t = table(prediction = stat.DE[,2], truth = stat.DE[,1])
  return(list(table = t,
              TPR = t[1,1]/sum(t[,1]),
              FPR = t[1,2]/sum(t[,2]),
              FDR = t[1,2]/sum(t[1,]),
              FNR = t[2,1]/sum(t[,1])))
}
