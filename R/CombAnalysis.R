#' Pipeline of Differential Expression Analysis for RNASeq Data
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample (only 2 groups allowed).
#' @param norm.method the method for normalization selected from \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param RUV_method the exact RUV method used from \code{RUVg}, \code{RUVr} and \code{RUVs} if \code{method = norm.RUV}
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#' @param DE.method the method for differential expression analysis from \code{DE.voom} and \code{DE.edgeR}, default to be \code{DE.voom}.
#' @param Pval p-value for identifying DE genes, default to be 0.01
#'
#' @return list, containing \code{id.list} (names of DE genes) and \code{p.val} (corresponding p-values and log2-fold-change for all genes).
#'
#' @export
#'
#' @examples res <- pip.norm.DE(data.test, data.group, "norm.TMM")
pip.norm.DE <- function(raw, groups, norm.method,
                        RUV_method = NULL, QN_filter = FALSE,
                        DE.method = "DE.voom", Pval = 0.01) {
  DEA <- function(raw, DE.method, normalized = TRUE, adjust = NULL){
    FUN = match.fun(DE.method)
    test.DE = FUN(RC = raw, groups = groups, P = Pval, normalized = normalized, adjust = adjust)
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

  if (norm.method %in% c("norm.SVA", "norm.RUV")) {
    if (norm.method == "norm.SVA") {
      test.norm = norm.SVA(raw, group = groups)$adjust.factor
    } else if (norm.method == "norm.RUV") {
      test.norm = norm.RUV(raw, group = groups, method = RUV_method)$adjust.factor
    }
    test.DE = DEA(raw = raw, DE.method = DE.method, normalized = FALSE, adjust = test.norm)
  }
  return(test.DE = test.DE)
}
