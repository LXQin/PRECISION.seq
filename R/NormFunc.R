#' Normalization By Trimmed Mean of M-values (TMM)
#'
#' Normalize the dataset using TMM, and return the normalized dataset with scaling factor scaled by 1e6.
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
#' test.TMM <-  norm.TMM(data.test, data.group)
norm.TMM <- function(raw, groups) {
  dat.DGE <- DGEList(counts = matrix(raw, ncol = length(groups)),
                     group = factor(groups),
                     genes = rownames(raw))
  d <- calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor <- d$samples$norm.factors * d$samples$lib.size / 1e6
  dat.normed <- t(t(raw)/scaling.factor)
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
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) scaled by 1e6 for each sample.
#'
#' @import edgeR
#' @export
#'
#' @references Dillies, Marie-Agnès, et al. "A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis." \emph{Briefings in bioinformatics} 14.6 (2013): 671-683.
#'
#' @examples
#' test.TC <-  norm.TC(data.test, data.group)
norm.TC <- function(raw, groups) {
  dat.DGE <- DGEList(counts = matrix(raw, ncol = length(groups)),
                     group = factor(groups),
                     genes = rownames(raw))
  scaling.factor <- dat.DGE$samples$lib.size/1e6
  dat.normed <- t(t(raw)/scaling.factor)
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
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) scaled by 1e6 for each sample.
#'
#' @import edgeR
#' @export
#'
#' @references Dillies, Marie-Agnès, et al. "A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis." \emph{Briefings in bioinformatics} 14.6 (2013): 671-683.
#'
#' @examples
#' test.UQ <- norm.UQ(data.test, data.group)
norm.UQ <- function(raw, groups) {
  dat.DGE <- DGEList(counts = matrix(raw, ncol = length(groups)),
                     group = factor(groups),
                     genes = rownames(raw))
  q.factor <- apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor <- q.factor/1e6
  dat.normed <- t(t(raw)/scaling.factor)
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
#' @return list, containing \code{dat.normed} (normalized dataset) and \code{scaling.factor} (scaling factor) scaled by 1e6 for each sample.
#'
#' @import edgeR
#' @export
#'
#' @references Dillies, Marie-Agnès, et al. "A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis." \emph{Briefings in bioinformatics} 14.6 (2013): 671-683.
#'
#' @examples
#' test.med <- norm.med(data.test, data.group)
norm.med <- function(raw, groups) {
  dat.DGE <- DGEList(counts = matrix(raw, ncol = length(groups)),
                     group = factor(groups),
                     genes = rownames(raw))
  m.factor <- apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor/1e6
  dat.normed <- t(t(raw)/scaling.factor)
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
#' @import DESeq2
#' @export
#'
#' @references \href{https://bioconductor.org/packages/release/bioc/vignettes/DESeq/inst/doc/DESeq.pdf}{DESeq Package}
#'
#' @examples
#' test.DESeq <- norm.DESeq(data.test, data.group)
norm.DESeq <- function(raw, groups) {
  condition <- data.frame(SampleName = colnames(raw), Condition = factor(groups))
  rownames(condition) = colnames(raw)
  dat.DGE <- DESeqDataSetFromMatrix(countData = raw, colData = condition, design = ~ Condition)
  dat.DGE <- estimateSizeFactors(dat.DGE)
  scaling.factor <- sizeFactors(dat.DGE)
  dat.normed <- DESeq2::counts(dat.DGE, normalized = T)
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
#' test.PoissonSeq <- norm.PoissonSeq(data.test)
norm.PoissonSeq <- function(raw) {
  invisible(capture.output(scaling.factor <- PS.Est.Depth(raw)))
  dat.normed <- t(t(raw)/scaling.factor)
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
#' test.QN <- norm.QN(data.test)
norm.QN <- function(raw, filter = FALSE) {
  if (filter == TRUE) {
    raw <- log2(raw + 1)
    raw <- raw[rowMeans(raw) > 2, ]
  } else {
    raw <- log2(raw + 1)
  }
  dat.log.normed <- normalize.quantiles(as.matrix(raw))
  dat.normed <- 2^dat.log.normed - 1
  colnames(dat.normed) <- colnames(raw)
  rownames(dat.normed) <- rownames(raw)
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
#' test.SVA <- norm.SVA(data.test, data.group)
norm.SVA <- function(raw, groups) {
  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.sva <- raw[filter, ]
  genes <- rownames(dat.sva)
  mod1 <- model.matrix(~ groups)
  mod0 <- cbind(mod1[,1])
  dat0 <- as.matrix(dat.sva)
#  svseq <- svaseq(dat0, mod1, mod0, n.sv = 1)$sv
  invisible(capture.output(svseq <- svaseq(dat0, mod1, mod0, n.sv = 1)$sv))
  adjust <- cbind(mod1, svseq)
  hat <- solve(t(adjust) %*% adjust) %*% t(adjust)
  beta <- (hat %*% t(raw))
  P <- ncol(mod1)
  dat.normed <- raw - t(as.matrix(adjust[,-c(1:P)]) %*% beta[-c(1:P),])
  return(list(dat.normed = dat.normed,
              adjust.factor = svseq))
}


#' Normalization By Remove Unwanted Variation Using Control Genes (RUVg)
#'
#' Normalize the dataset using RUV Using Control Genes, and return the normalized dataset with adjusting factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.

#' @return list, containing \code{dat.normed} (normalized dataset), and the \code{adjust.factor} (adjusting factors) for the design matrix. The normalized dataset could only used for exploration, and adjusting factors are recommended as a covariate in the downstream analysis.
#'
#' @import EDASeq
#' @import DESeq2
#' @import edgeR
#' @import RUVSeq
#' @import Biobase
#' @import BiocGenerics
#' @export
#'
#' @references \href{http://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf}{RUVSeq Tutorial}
#'
#' @examples
#' test.RUVg <- norm.RUV(data.test, data.group)
norm.RUVg <- function(raw, groups) {
  if (!require("Biobase")) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }

  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  genes <- rownames(dat.ruv)
  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv),
                             phenoData = data.frame(condition,
                                                    row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition,
                                                        row.names = colnames(dat.ruv)))
  y <- DGEList(counts = DESeq2::counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15*nrow(raw))]))]


  t <- RUVg(set, spikes, k = 1)
  dat.normed <- normCounts(t)
  return(list(dat.normed = dat.normed,
              adjust.factor = t$W))
}


#' Normalization By Remove Unwanted Variation Using Replicate Samples (RUVs)
#'
#' Normalize the dataset using RUV using replicate samples, and return the normalized dataset with adjusting factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.

#' @return list, containing \code{dat.normed} (normalized dataset), and the \code{adjust.factor} (adjusting factors) for the design matrix. The normalized dataset could only used for exploration, and adjusting factors are recommended as a covariate in the downstream analysis.
#'
#' @import EDASeq
#' @import DESeq2
#' @import edgeR
#' @import RUVSeq
#' @import Biobase
#' @import BiocGenerics
#' @export
#'
#' @references \href{http://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf}{RUVSeq Tutorial}
#'
#' @examples
#' test.RUVs <- norm.RUVs(data.test, data.group)
norm.RUVs <- function(raw, groups) {
  if (!require("Biobase")) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }

  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  genes <- rownames(dat.ruv)
  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv),
                             phenoData = data.frame(condition,
                                                    row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition,
                                                        row.names = colnames(dat.ruv)))
  y <- DGEList(counts = counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15*nrow(raw))]))]

  differences <- makeGroups(condition)
  controls <- rownames(dat.ruv)
  t <- RUVs(set, controls, k = 1, differences)
  dat.normed <- normCounts(t)
  return(list(dat.normed = dat.normed,
              adjust.factor = t$W))
}


#' Normalization By Remove Unwanted Variation Using Residuals (RUVr)
#'
#' Normalize the dataset using RUV using residuals, and return the normalized dataset with adjusting factor.
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.

#' @return list, containing \code{dat.normed} (normalized dataset), and the \code{adjust.factor} (adjusting factors) for the design matrix. The normalized dataset could only used for exploration, and adjusting factors are recommended as a covariate in the downstream analysis.
#'
#' @import EDASeq
#' @import DESeq2
#' @import edgeR
#' @import RUVSeq
#' @import Biobase
#' @import BiocGenerics
#' @export
#'
#' @references \href{http://www.bioconductor.org/packages/devel/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf}{RUVSeq Tutorial}
#'
#' @examples
#' test.RUVr <- norm.RUVr(data.test, data.group)
norm.RUVr <- function(raw, groups) {
  if (!require("Biobase")) {
    stop("Package \"Biobase\" needed for this function to work. Please install it.",
         call. = FALSE
    )
  }

  filter <- apply(raw, 1, function(x) length(x[x > 5]) >= 2)
  dat.ruv <- raw[filter, ]
  genes <- rownames(dat.ruv)
  condition <- factor(groups)
  set <- newSeqExpressionSet(as.matrix(dat.ruv),
                             phenoData = data.frame(condition,
                                                    row.names = colnames(dat.ruv)))
  design <- model.matrix(~ condition, data = data.frame(condition,
                                                        row.names = colnames(dat.ruv)))
  y <- DGEList(counts = DESeq2::counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  top <- topTags(lrt, n = nrow(set))$table
  spikes <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:(0.15*nrow(raw))]))]

  design <- model.matrix(~ condition, data = pData(set))
  y <- DGEList(counts = DESeq2::counts(set), group = condition)
  y <- calcNormFactors(y, method = "upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  res <- residuals(fit, type = "deviance")
  setUQ <- betweenLaneNormalization(set, which = "upper")
  controls <- rownames(dat.ruv)
  t <- RUVr(setUQ, controls, k = 1, res)
  dat.normed <- normCounts(t)

  return(list(dat.normed = dat.normed,
              adjust.factor = t$W))
}


