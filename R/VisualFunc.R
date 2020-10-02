#' Volcano Figure
#'
#' Function for generating volcano figures with log fold change as x-axis and -log10(p-value) as y-axis.
#'
#' @param DEA.res the list of differential expression analysis results obtained from DE.voom, DE.edge, or any method from the users by storing the results as same as DE methods in the package (including DE genes, p-values and log2 fold changes).
#' @param title the title of the figure
#'
#' @return volcano plot
#' @export
#'
#' @examples
#' voom.benchmark <- DE.voom(data.benchmark, data.group)
#' fig.volcano(voom.benchmark, title = "Volcano Plot")
fig.volcano <- function(DEA.res, title){
  DE.list <- DEA.res$id.list
  dat.DE.frame <- data.frame(dm = DEA.res$log2.FC,
                             p.value = DEA.res$p.val)
  mask <- with(dat.DE.frame, p.value < .01)
  cols <- ifelse(mask, "red", "black")

  xlim = round(1.15 * max(abs(dat.DE.frame$dm)), 1)
  ylim = 1.15 * round(max(-log10(dat.DE.frame$p.value)), 1)
  with(dat.DE.frame, plot(dm, -log10(p.value), cex = .5, pch = 16,
                          col = cols, xlim = c(-xlim, xlim),
                          ylim = c(0, ylim),
                          xlab = "Mean Difference",
                          main = title))
  abline(h = 2, lty = 2)
}


#' Relative Log Expression Plot
#'
#' Function for generating relative log expression plot using the read count data as the input.
#'
#' @param data normalized testing using the method which the researchers desire to compare with the other methods
#' @param groups vector of characters indicating the group for each sample (only 2 groups are allowed).
#' @param title the title of the figure
#'
#' @return boxplot for relative log expression
#' @export
#'
#' @examples
#' fig.RLE(data.test, data.group, "test without normalization")
fig.RLE = function(data, groups, title) {
  raw.log <- log2(data + 1)
  rle <- t(apply(raw.log, 1, function(x) x - median(x)))
  color <- groups
  color[color == levels(factor(color))[1]] <- rainbow(2)[1]
  color[color == levels(factor(color))[2]] <- rainbow(2)[2]
  ylim = round(max(1.5*(apply(rle, 2, IQR)) + max(abs(apply(rle, 2, function(x) quantile(x, c(0.25, 0.75)))))), 1)
  boxplot(rle, col = color, ylab = "RLE", ylim = c(-ylim, ylim),
          outline = FALSE, xaxt = "n", main = title)
  legend("topright",levels(factor(groups)), bty = "n",
         pch = "x", cex = 1, col = c(rainbow(2)[1], rainbow(2)[2]))
}


#' Concordance At The Top Plot
#'
#' Function for generating concordance at the top plot, which compares concordance of the p-values obtained
#' from benchmark data without normalization and normalized test data.
#'
#'
#' @param MethodsCompare the vector of methods that researchers would like to compare with, and selected from  \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param benchmark benchmark data with the same format of data.benchmark, which are used for benchmark p-values, and is set to use data.benchmark by default
#' @param test test data with the same format of data.test, which are used for test p-values, and is set to use data.test by default
#' @param groups vector of characters indicating the group for each sample.
#' @param DE.method the DEA method selected from \code{DE.voom} and \code{DE.edgeR}, and default to be \code{DE.voom}.
#' @param QN_filter whether the filtering is performed if \code{norm.QN} is used.
#' @param Pval p-value cut-off point for identifying DE genes, default to be 0.01.
#' @param MethodNew the name of the new method which will be presented in the legend.
#' @param pvalues a vector of p values computed from DE analysis using the new method normalized test data.
#' @param Methods_visual a vector of methods' names shown in the figure.
#'
#' @return figure of concordance for comparison
#'
#' @import ffpe
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' t <-  runif(1033)
#' names(t) <-  rownames(data.test)
#' fig.CAT(MethodsCompare = c("norm.none", "norm.TMM", "norm.SVA", "norm.TC", "norm.RUVr"), MethodNew = "Example", pvalues = t)
fig.CAT <- function(MethodsCompare, benchmark = data.benchmark, test = data.test,
                    group = data.group, DE.method = "DE.voom", QN_filter = FALSE,
                    Pval = 0.01, MethodNew, pvalues, Methods_visual = c(MethodsCompare, MethodNew)){
  DEFUN = match.fun(DE.method)
  benchmark.p <-  DEFUN(benchmark, group)$p.val
  numMethods <- length(MethodsCompare) + 1
  catplots <- vector("list", numMethods)
  for (i in 1:(numMethods - 1)) {
    catplots[[i]] <- CATplot(pip.norm.DE(test, group, MethodsCompare[i],
                                         QN_filter = QN_filter,
                                         DE.method = DE.method, Pval = Pval)$p.val,
                             benchmark.p, maxrank = 100, make.plot = F)
  }
  catplots[[numMethods]] <- CATplot(pvalues, benchmark.p, maxrank = 100, make.plot = F)

  data_cat = data.frame(x = numeric(), y = numeric(), curve = character())
  for(i in 1:numMethods){
    y = catplots[[i]][,2]
    x = 1:length(y)
    data_cat = rbind(data_cat, data.frame(x, y, curve = Methods_visual[i]))
  }
  data_cat$curve = factor(data_cat$curve)

  return(ggplot(data_cat, aes(x, y, color = curve)) +
           geom_line() +
           ylab("Rate of Agreement with Benchmark") +
           xlab("Significance Rank") +
           ggtitle("CATplot") +
           theme_bw())
}


#' Selection of normalization methods based on golden standards
#'
#' @param raw raw count data in the format of data frame or matrix, with columns for samples and raws for genes.
#' @param groups vector of characters indicating the group for each sample (only 2 groups allowed).
#' @param truth vector of genes that are truly differential expressed
#' @param MethodsCompare the vector of methods that researchers would like to compare with, and selected from  \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#' @param DE.method the method for differential expression analysis from \code{DE.voom} and \code{DE.edgeR}, default to be \code{DE.voom}.
#' @param Pval p-value for identifying DE genes, default to be 0.01
#' @param selected.marker if given, vector of a subset of genes/markers for
#' this analysis. Leave \code{NULL} if all markers are considered for the analysis.
#'
#' @return figure for selection of normalization methods
#'
#' @import ggrepel
#'
#' @export
#'
#' @examples
#' truthgene <- DE.voom(data.benchmark, data.group)$id.list
#' fig.FDR_FNR(data.test, data.group, MethodsCompare = c("norm.none", "norm.TMM", "norm.SVA", "norm.TC"), truth = truthgene)
fig.FDR_FNR <- function(raw, groups, MethodsCompare, truth,
                        QN_filter = FALSE,
                        DE.method = "DE.voom", Pval = 0.01,
                        selected.marker = NULL) {
  numMethods <- length(MethodsCompare)
  FNR <- FDR <- c()
  for (i in 1:numMethods){
    temp <- pip.statistics(raw = raw, groups = groups, truth = truth,
                           DE.method = DE.method,
                           norm.method = MethodsCompare[i],
                           QN_filter = QN_filter, Pval = Pval,
                           selected.marker = selected.marker)
    FNR <- c(FNR, temp$FNR)
    FDR <- c(FDR, temp$FDR)
  }
  stat <- data.frame(FNR = FNR, FDR = FDR, names = MethodsCompare)
  p <- ggplot(stat, aes(FNR, FDR, label = names)) + geom_point(color = "red") +
    geom_text_repel() + theme_bw()

  return(p)
}


#' Venn diagram for p-values
#'
#' Venn diagram is used to identify the performance of different normalization methods based on intersection of differential expressed genes.
#' @param benchmark.pval p-values obtained from benchmark data
#' @param test.pval p-values obtained from normalized test data
#' @param Pvalue Cut-off point for p-values for identifying significance
#'
#' @return A Venn diagram
#' @import limma
#' @export
#'
#' @examples
#' benchmark.voom <- DE.voom(RC = data.benchmark, groups = data.group, P = 0.01)
#' test.voom <- DE.voom(RC = data.test, groups = data.group, P = 0.01)
#' fig.venn(benchmark.voom$p.val, test.voom$p.val, 0.01)
fig.venn <- function(benchmark.pval, test.pval, Pvalue){
  bench.sig <- (benchmark.pval < Pvalue)
  test.sig <- (test.pval[names(benchmark.pval)] < Pvalue)
  venn2 <- cbind(bench.sig, test.sig)
  p <- vennDiagram(vennCounts(venn2),
                   names = c("Benchmark", "Test"),
                   cex = 1.5, counts.col = rainbow(1))
  return(p)
}


#' Dendrogram for clustering p-values
#'
#' Function for clustering normalization methods based on the p-values pattern calculated from the same dataset.
#'
#' @param MethodsCompare the vector of methods that researchers would like to compare with, and selected from  \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param DE.method the DEA method selected from \code{DE.voom} and \code{DE.edgeR}, and default to be \code{DE.voom}.
#' @param QN_filter whether the filtering is performed if \code{norm.QN} is used.
#' @param MethodNew the name of the new method which will be presented in the legend.
#' @param pvalues a vector of p values computed from DE analysis using the new method normalized test data.
#' @param test test data with the same format of data.test, which are used for test p-values, and is set to use data.test by default
#' @param group vector of characters indicating the group for each sample (only 2 groups allowed).
#'
#' @return figure of dendrogram
#'
#' @import ggdendro
#' @import Biobase
#' @export
#'
#' @examples
#' t <-  runif(1033)
#' names(t) <-  rownames(data.test)
#' fig.dendrogram(MethodsCompare = c("norm.none", "norm.TMM", "norm.SVA", "norm.TC", "norm.RUVr"), MethodNew = "Example", pvalues = t)
fig.dendrogram <- function(MethodsCompare, DE.method = "DE.voom", QN_filter = FALSE,
                           Pval = 0.01, MethodNew, pvalues, Methods_visual = c(MethodsCompare, MethodNew),
                           test = data.test, group = data.group){
  pval.index = rownames(test)
  numMethods <- length(MethodsCompare) + 1
  pval_list <- vector("list", numMethods)
  for (i in 1:(numMethods - 1)) {
    pval_list[[i]] <- pip.norm.DE(test, group, MethodsCompare[i],
                                  QN_filter = QN_filter,
                                  DE.method = DE.method, Pval = Pval)$p.val[pval.index]
  }
  pval_list[[numMethods]] <- pvalues[pval.index]
  pval_frame <- data.frame(matrix(unlist(pval_list), nrow = length(pval.index), byrow=T))
  colnames(pval_frame) <- c(MethodsCompare, MethodNew)
  dendrogram.p <- t(-log10(pval_frame))
  ds <- dist(dendrogram.p)
  hc <- hclust(ds, method = "ward.D")
  p <- ggdendrogram(hc, rotate = FALSE, size = 2)
  return(p)
}
