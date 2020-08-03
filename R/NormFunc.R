#' Normalization By Trimmed Mean of M-values (TMM)
#'
#' Normalize the dataset using TMM, and return the normalized dataset with scaling factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#'
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) for each sample.
#'
#' @import edgeR
#' @export
#'
#' @references \href{https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf}{edgeR User Guide}
#'
#' @examples
#' group = c(rep(c(rep('MXF',9),rep('PMFH',9)),3))
#' test.TMM = norm.TMM(test, group)
norm.TMM = function(raw, groups) {
  dat.DGE = DGEList(counts = matrix(raw, ncol = length(groups)),
                    group = factor(groups),
                    genes = rownames(raw))
  d = calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor = d$samples$norm.factors * d$samples$lib.size / 1e6
  dat.normed = t(t(raw)/scaling.factor)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


#' Normalization By Total Count (TC)
#'
#' Normalize the dataset using TC, and return the normalized dataset with scaling factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#'
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) for each sample.
#'
#' @import edgeR
#' @export
#'
#' @references Dillies, Marie-Agnès, et al. "A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis." \emph{Briefings in bioinformatics} 14.6 (2013): 671-683.
#'
#' @examples
#' group = c(rep(c(rep('MXF',9),rep('PMFH',9)),3))
#' test.TC = norm.TC(test, group)
norm.TC = function(raw, groups) {
  dat.DGE = DGEList(counts = matrix(raw, ncol = length(groups)),
                    group = factor(groups),
                    genes = rownames(raw))
  scaling.factor = dat.DGE$samples$lib.size/1e6
  dat.normed = t(t(raw)/scaling.factor)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


#' Normalization By Upper Quantile (UQ)
#'
#' Normalize the dataset using UQ, and return the normalized dataset with scaling factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#'
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) for each sample.
#'
#' @import edgeR
#' @export
#'
#' @references Dillies, Marie-Agnès, et al. "A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis." \emph{Briefings in bioinformatics} 14.6 (2013): 671-683.
#'
#' @examples
#' group = c(rep(c(rep('MXF',9),rep('PMFH',9)),3))
#' test.UQ = norm.UQ(test, group)
norm.UQ = function(raw, groups) {
  dat.DGE = DGEList(counts = matrix(raw, ncol = length(groups)),
                    group = factor(groups),
                    genes = rownames(raw))
  q.factor = apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor = q.factor/1e6
  dat.normed = t(t(raw)/scaling.factor)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


#' Normalization By Median (Med)
#'
#' Normalize the dataset using Med, and return the normalized dataset with scaling factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#'
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) for each sample.
#'
#' @import edgeR
#' @export
#'
#' @references Dillies, Marie-Agnès, et al. "A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis." \emph{Briefings in bioinformatics} 14.6 (2013): 671-683.
#'
#' @examples
#' group = c(rep(c(rep('MXF',9),rep('PMFH',9)),3))
#' test.med = norm.med(test, group)
norm.med = function(raw, groups) {
  dat.DGE = DGEList(counts = matrix(raw, ncol = length(groups)),
                    group = factor(groups),
                    genes = rownames(raw))
  m.factor = apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor = m.factor/1e6
  dat.normed = t(t(raw)/scaling.factor)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


#' Normalization By DESeq (DESeq)
#'
#' Normalize the dataset using DESeq, and return the normalized dataset with scaling factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#'
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) for each sample.
#'
#' @import DESeq
#' @export
#'
#' @references \href{https://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf}{DESeq Package}
#'
#' @examples
#' group = c(rep(c(rep('MXF',9),rep('PMFH',9)),3))
#' test.DESeq = norm.DESeq(test, group)
norm.DESeq = function(raw, groups) {
  condition = factor(groups)
  dat.DGE = estimateSizeFactors(newCountDataSet(raw, condition))
  scaling.factor = sizeFactors(dat.DGE)
  dat.normed = counts(dat.DGE, normalized = T)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


#' Normalization By PoissonSeq (PoissonSeq)
#'
#' Normalize the dataset using PoissonSeq, and return the normalized dataset with scaling factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#'
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) for each sample.
#'
#' @import PoissonSeq
#' @export
#'
#' @references \href{https://cran.r-project.org/web/packages/PoissonSeq/PoissonSeq.pdf}{PossionSeq Package}
#'
#' @examples
#' test.PoissonSeq = norm.PoissonSeq(test)
norm.PoissonSeq = function(raw) {
  scaling.factor = PS.Est.Depth(raw)
  dat.normed = t(t(raw)/scaling.factor)
  return(list(dat.normed = dat.normed,
              scaling.factor = scaling.factor))
}


#' Normalization By Quantile Normalization (QN)
#'
#' Normalize the dataset using QN, and return the normalized dataset with scaling factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param filter whether filter before normalization
#'
#' @return list, containing \code{dat.normed} (normalized dataset).
#'
#' @import preprocessCore
#' @export
#'
#' @references \href{http://jtleek.com/genstats/inst/doc/02_05_normalization.html}{Quantile Normalization Tutorial}
#'
#' @examples
#' test.QN = norm.QN(test)
norm.QN = function(raw, filter = FALSE) {
  if (filter == TRUE) {
    raw = log2(raw + 1)
    raw = raw[rowMeans(raw) > 2, ]
  } else {
    raw = log2(raw + 1)
  }
  dat.log.normed = normalize.quantiles(as.matrix(raw))
  dat.normed = 2^dat.log.normed - 1
  colnames(dat.normed) = colnames(raw)
  rownames(dat.normed) = rownames(raw)
  return(list(dat.normed = dat.normed))
}


#' Normalization By Surrogate Variable Analysis for Sequencing Data (SVA)
#'
#' Normalize the dataset using SVA, and return the normalized dataset with adjusting factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#' @param calibrator whether the raw dataset includes the known negative control genes which should contain "cali" in the gene names.
#'
#' @return list, containing \code{dat.normed} (normalized dataset) without the calibrators, and the \code{adjust.factor} (adjusting factors) for the design matrix. The normalized dataset could only used for exploration.
#'
#' @import sva
#' @export
#'
#' @references \href{https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf}{SVA Package}
#'
#' @examples
#' group = c(rep(c(rep('MXF',9),rep('PMFH',9)),3))
#' test.SVA = norm.SVA(test, group)
norm.SVA = function(raw, groups) {
  filter = apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.sva = raw[filter, ]
  genes = rownames(dat.sva)
  mod1 = model.matrix(~ groups)
  mod0 = cbind(mod1[,1])
  dat0 = as.matrix(dat.sva)
  svseq = svaseq(dat0, mod1, mod0, n.sv = 1)$sv
  adjust = cbind(mod1, svseq)
  hat = solve(t(adjust) %*% adjust) %*% t(adjust)
  beta = (hat %*% t(raw))
  P = ncol(mod1)
  dat.normed = raw - t(as.matrix(adjust[,-c(1:P)]) %*% beta[-c(1:P),])
  return(list(dat.normed = dat.normed,
              adjust.factor = svseq))
}


#' Normalization By Remove Unwanted Variation (RUV)
#'
#' Normalize the dataset using RUV (including RUVg, RUVr and RUVs), and return the normalized dataset with adjusting factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#' @param method indicate the exact RUV method used from \code{RUVg}, \code{RUVr} and \code{RUVs}.

#' @return list, containing \code{dat.normed} (normalized dataset), and the \code{adjust.factor} (adjusting factors) for the design matrix. The normalized dataset could only used for exploration, and adjusting factors are recommended as a covariate in the downstream analysis.
#'
#' @import EDASeq
#' @import DESeq
#' @import edgeR
#' @import RUVSeq
#' @export
#'
#' @references \href{http://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf}{RUVSeq Tutorial}
#'
#' @examples
#' group = c(rep(c(rep('MXF',9),rep('PMFH',9)),3))
#' test.RUVr = norm.RUV(test, group, method = "RUVr")
norm.RUV = function(raw, groups, method = c("RUVg", "RUVs", "RUVr")) {
  filter = apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv = raw[filter, ]
  genes = rownames(dat.ruv)
  condition = factor(groups)
  set = newSeqExpressionSet(as.matrix(dat.ruv),
                            phenoData = data.frame(condition,
                                                   row.names = colnames(dat.ruv)))
  design = model.matrix(~ condition, data = data.frame(condition,
                                                       row.names = colnames(dat.ruv)))
  y = DGEList(counts = counts(set), group = condition)
  y = calcNormFactors(y, method = "upperquartile")
  y = estimateGLMCommonDisp(y, design)
  y = estimateGLMTagwiseDisp(y, design)
  fit = glmFit(y, design)
  lrt = glmLRT(fit, coef = 2)
  top = topTags(lrt, n = nrow(set))$table
  spikes = rownames(set)[which(!(rownames(set) %in% rownames(top)[1:0.15*nrow(raw)]))]

  if (method == "RUVg") {
    t = RUVg(set, spikes, k = 1)
    dat.normed = normCounts(t)
    return(list(dat.normed = dat.normed,
                adjust.factor = t$W))
  }else if (method == "RUVs") {
    differences = makeGroups(condition)
    controls = rownames(dat.ruv)
    t = RUVs(set, controls, k = 1, differences)
    dat.normed = normCounts(t)
    return(list(dat.normed = dat.normed,
                adjust.factor = t$W))
  }else if (method == "RUVr") {
    design = model.matrix(~ condition, data = pData(set))
    y = DGEList(counts = counts(set), group = condition)
    y = calcNormFactors(y, method = "upperquartile")
    y = estimateGLMCommonDisp(y, design)
    y = estimateGLMTagwiseDisp(y, design)
    fit = glmFit(y, design)
    res = residuals(fit, type = "deviance")
    setUQ = betweenLaneNormalization(set, which = "upper")
    controls = rownames(dat.ruv)
    t = RUVr(setUQ, controls, k = 1, res)
    dat.normed = normCounts(t)

    return(list(dat.normed = dat.normed,
                adjust.factor = t$W))
  }
}
