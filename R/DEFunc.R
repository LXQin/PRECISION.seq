#' Differential Expression Analysis Using Voom-limma Pipeline
#'
#' Perform DEA using the voom-limma pipeline on a normalized dataset.
#' The normalized data can be provided as normalized counts or by adjusting
#' factor for the original count data.
#'
#' @param RC Data in the format of a data frame or matrix, with columns for samples and rows for genes.
#' @param groups Vector of characters indicating the group label for each sample.
#' @param Pval Cut-off point for p-values to identify differentially expressed genes.
#' @param normalized Logical, whether the data is provided as normalized counts.
#' If set to FALSE, adjustment factors must be provided using the \code{adjust}
#' parameter.
#' @param adjust Adjusting factors for normalizing the count data. Must be
#' provided if the data is not normalized beforehand as indicated by the
#' \code{normalized} parameter.
#'
#' @return list, containing \code{id.list} (names of DE genes), \code{p.val}, and \code{log2.FC}.
#'
#' @import edgeR
#' @import limma
#' @export
#'
#' @references \href{https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html}{Differential Expression with Limma-Voom}
#'
#' @examples
#' voom.benchmark <- DE.voom(data.benchmark, data.group)
#' str(voom.benchmark)
DE.voom <- function(RC, groups, Pval = 0.01, normalized = TRUE, adjust = NULL) {
  event <- factor(groups)
  if (normalized == TRUE) {
    design <- model.matrix(~ 0 + event)
    colnames(design) <- levels(event)
  } else {
    if(is.null(adjust)) {
      stop("Error in DE.voom: If the data is not normalized (normalized=FALSE), adjustment factors must be provided using the adjust parameter.")
    }
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
              p.val = out[,1],
              log2.FC = out[,2]))
}


#' Differential Expression Analysis Using EdgeR
#'
#' Perform DEA by EdgeR on normalized dataset or dataset with adjusting factor.
#' The normalized data can be provided as normalized counts or by adjusting
#' factor for the original count data.
#'
#' @param RC data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample.
#' @param Pval cut-off point for p-values to identify DE genes.
#' @param normalized Logical, whether the data is provided as normalized counts.
#' If set to FALSE, adjustment factors must be provided using the \code{adjust}
#' parameter.
#' @param adjust Adjusting factors for normalizing the count data. Must be
#' provided if the data is not normalized beforehand as indicated by the
#' \code{normalized} parameter.
#'
#' @return list, containing \code{id.list} (names of DE genes), \code{p.val}, and \code{log2.FC}.
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
#' str(edgeR.benchmark)
DE.edgeR <- function(RC, groups, Pval = 0.01, normalized = TRUE, adjust = NULL) {
  event <- factor(groups);
  if (normalized == TRUE) {
    design <- model.matrix(~ 0 + event)
    colnames(design) <- levels(event)
  } else {
    if(is.null(adjust)) {
      stop("Error in DE.voom: If the data is not normalized (normalized=FALSE), adjustment factors must be provided using the adjust parameter.")
    }
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
              p.val = out[,1],
              log2.FC = out[,2]))
}



#' Statistics for DEA Results Based Golden Standards
#'
#' Computes the agreement of differential expression (DE) wihtin one data set
#' with an assumed true DE status usually given based on a gold standard.
#' The following statistics are computed based of the differentially expressed
#' markers \code{id.list} compared to the list truly differentially expressed
#' markers \code{truth}:
#' \itemize{
#'   \item False Discovery Rate (FDR)
#'   \item False Negative Rate (FNR)
#'   \item True Positive Rate (TPR)
#'   \item False Positive Rate (FPR)
#' }
#'
#' Both, \code{id.list} and \code{truth} are a subset of markers from
#' \code{markers}.
#'
#' @param markers Vector of all markers considered in the analysis.
#' @param id.list Vector of markers that are differentially expressed in a
#' data set --- typically in the data set \code{\link{data.test}}.
#' Differential expression can be computed via \code{\link{DE.edgeR}}
#' or \code{\link{DE.voom}}.
#' @param truth Vector of genes that are (assumedly) truly differential
#' expressed, e.g. basen on a gold standard --- typically the data set
#' \code{\link{data.benchmark}}.
#' @param selected.marker optional Vector of a subset of \code{markers}.
#' If given, the analysis will be limited to the given subset.
#' Leave \code{NULL} if all markers are considered for the analysis.
#'
#' @return A list of:
#' \describe{
#'   \item{TPR}{True positive Rate}
#'   \item{FPR}{False Positive Rate}
#'   \item{FDR}{False Discovery Rate}
#'   \item{FNR}{False Negative Rate}
#' }
#'
#' @export
#' @examples
#' DE.bench <- DE.voom(data.benchmark, data.group)
#' DE.test <- DE.voom(data.test, data.group)
#' stats <- DE.statistics(rownames(data.benchmark), DE.test$id.list, DE.bench$id.list)
#' print(stats)
DE.statistics <- function(markers, id.list, truth, selected.marker=NULL) {

  if (!is.null(selected.marker)) {
    id.list = id.list[id.list %in% selected.marker]
    truth = truth[truth %in% selected.marker]
    markers = selected.marker
  }

  tp = sum(id.list %in% truth)
  fp = length(id.list) - tp
  fn = length(truth) - tp
  tn = length(markers) - length(truth) - fp

  return(list(TPR = tp/(tp + fn),
              FPR = fp/(fp + tn),
              FDR = fp/(tp + fp),
              FNR = fn/(tp + fn)))
}




