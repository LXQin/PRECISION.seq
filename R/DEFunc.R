#' Differential Expression Analysis Using Voom-limma Pipeline
#'
#' Perform DEA by voom-limma pipeline on normalized dataset or dataset with adjusting factor.
#'
#' @param RC data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#' @param Pval cut-off point for p-values to identify DE genes.
#' @param normalized whether the dataset is normalized (by scaling method)
#' @param adjust the adjusting factor if the dataset itself is not normalized
#'
#' @return list, containing \code{id.list} (names of DE genes) and \code{p.val} (corresponding p-values and log2-fold-change for all genes).
#'
#' @import edgeR
#' @import limma
#' @export
#'
#' @references \href{https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html}{Differential Expression with Limma-Voom}
#'
#' @examples
#' voom.benchmark <- DE.voom(data.benchmark, data.group)
DE.voom <- function(RC, groups, Pval = 0.01, normalized = TRUE, adjust = NULL) {
  event <- factor(groups)
  if (normalized == TRUE) {
    design <- model.matrix(~ 0 + event)
    colnames(design) <- levels(event)
    } else {
    design <- model.matrix(~ 0 + event + adjust)
    colnames(design)[1:2] <- levels(event)
  }
  contr <- paste(levels(event)[2], "-", levels(event)[1])
  contrast.mx <- limma::makeContrasts(contrasts = contr, levels = design)
  d <- DGEList(RC, genes = rownames(RC))
  d.voom <- voom(d, design)
  fit <- lmFit(d.voom, design)
  fit.contr <- contrasts.fit(fit, contrast.mx)
  fit.eb <- eBayes(fit.contr)
  P.value <- fit.eb$p.value
  fc.log2 <- fit.eb$coef

  out <- cbind(P.value,fc.log2)
  out <- out[order(out[,1]),]
  colnames(out) <- c("Pvalue","log2.FC")
  return(list(id.list = rownames(out[out[,1] < Pval,]),
              p.val = out))
}

#' Differential Expression Analysis Using EdgeR
#'
#' Perform DEA by EdgeR on normalized dataset or dataset with adjusting factor.
#'
#' @param RC data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#' @param Pval cut-off point for p-values to identify DE genes.
#' @param normalized whether the dataset is normalized (by scaling method)
#' @param adjust the adjusting factor if the dataset itself is not normalized
#'
#' @return list, containing \code{id.list} (names of DE genes) and \code{p.val} (corresponding p-values and log2-fold-change for all genes).
#'
#' @import edgeR
#' @import limma
#' @export
#'
#' @references
#' \url{https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html}
#'
#' @examples
#' edgeR.benchmark <- DE.edgeR(data.benchmark, data.group)
DE.edgeR <- function(RC, groups, Pval = 0.01, normalized = TRUE, adjust = NULL) {
  event <- factor(groups);
  if (normalized == TRUE) {
    design <- model.matrix(~ 0 + event)
    colnames(design) <- levels(event)
  } else {
    design <- model.matrix(~ 0 + event + adjust)
    colnames(design)[1:2] <- levels(event)
  }
  d <- DGEList(RC, group = event, genes = rownames(RC))
  d <- estimateDisp(d, design)
  DE <- exactTest(d, pair = levels(d$samples$group)[1:2])
  out <- DE$table[,c(3,1)]
  out <- out[order(out[,1]),]
  colnames(out) <- c("Pvalue","log2.FC")
  return(list(id.list = rownames(out[out[,1] < Pval,]),
              p.val = out))
}
