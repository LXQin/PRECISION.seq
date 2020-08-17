#' Volcano Figure
#'
#' Function for generating volcano figures with log fold change as x-axis and -log10(p-value) as y-axis.
#'
#' @param normtest normalized data using the method which the researchers desire to compare with the other methods, or it could be the raw data if the adjusting parameters available. The rows are for genes/markers, and the columns are for samples.
#' @param groups vector of characters indicating the group for each sample (only 2 groups are allowed).
#' @param adjust vector of adjusting parameters for each sample.
#' @param main the title of the figure
#'
#' @return volcano plot
#' @export
#'
#' @examples
#' fig.volcano(data.benchmark, data.group, main = "Volcano Plot")
fig.volcano <- function(normtest, groups, adjust = NULL, main){
  if (is.null(adjust)) {
    dat.voom = DE.voom(RC = normtest, groups = groups, P = 0.01)
  } else {
      dat.voom = DE.voom(RC = normtest, groups = groups, P = 0.01, adjust = adjust)
  }

  DE.list <- dat.voom$id.list
  dat.voom.frame <- data.frame(dm = dat.voom$p.val[,2],
                               p.value = dat.voom$p.val[,1])
  mask <- with(dat.voom.frame, p.value < .01)
  cols <- ifelse(mask, "red", "black")

  xlim = round(1.15 * max(abs(dat.voom.frame$dm)), 1)
  ylim = 1.15 * round(max(-log10(dat.voom.frame$p.value)), 1)
  with(dat.voom.frame, plot(dm, -log10(p.value), cex = .5, pch = 16,
                            col = cols, xlim = c(-xlim, xlim),
                            ylim = c(0, ylim),
                            xlab = "Mean Difference",
                            main = main))
  abline(h = 2, lty = 2)
}


#' Relative Log Expression Plot
#'
#' Function for generating relative log expression plot using the normalized test data as the input.
#'
#' @param normtest normalized testing using the method which the researchers desire to compare with the other methods
#' @param groups vector of characters indicating the group for each sample (only 2 groups are allowed).
#' @param main the title of the figure
#'
#' @return boxplot for relative log expression
#' @export
#'
#' @examples
#' fig.RLE(data.test, data.group, "test without normalization")
fig.RLE = function(normtest, groups, main) {
  raw.log <- log2(normtest + 1)
  rle <- t(apply(raw.log, 1, function(x) x - median(x)))
  color <- groups
  color[color == levels(factor(color))[1]] <- rainbow(2)[1]
  color[color == levels(factor(color))[2]] <- rainbow(2)[2]
  ylim = round(max(1.5*(apply(rle, 2, IQR)) + max(abs(apply(rle, 2, function(x) quantile(x, c(0.25, 0.75)))))), 1)
  boxplot(rle, col = color, ylab = "RLE", ylim = c(-ylim, ylim),
          outline = FALSE, xaxt = "n", main = main)
  legend("topright",levels(factor(groups)), bty = "n",
         pch = "x", cex = 1, col = c(rainbow(2)[1], rainbow(2)[2]))
}


#' Concordance At The Top Plot
#'
#' Function for generating concordance at the top plot, which compares concordance of the p-values obtained
#' from benchmark data without normalization and normalized test data.
#'
#' @param MethodsCompare the vector of methods that researchers would like to compare with, and selected from  \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param DE.method the DEA method selected from \code{DE.voom} and \code{DE.edgeR}, and default to be \code{DE.voom}.
#' @param RUV_method the exact RUV method used from \code{RUVg}, \code{RUVr} and \code{RUVs} if \code{norm.RUV} is used.
#' @param QN_filter whether the filtering is performed if \code{norm.QN} is used.
#' @param Pval p-value for identifying DE genes, default to be 0.01.
#' @param MethodNew the name of the new method which will be presented in the legend.
#' @param pvalues a vector of p values computed from DE analysis using the new method normalized data.
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
#' fig.CAT(MethodsCompare = c("norm.none", "norm.TMM", "norm.SVA", "norm.TC"), MethodNew = "Example", pvalues = t)
fig.CAT <- function(MethodsCompare, DE.method = "DE.voom", RUV_method = NULL, QN_filter = FALSE,
                    Pval = 0.01, MethodNew, pvalues, Methods_visual = c(MethodsCompare, MethodNew)){
  DEFUN = match.fun(DE.method)
  benchmark.p <-  DEFUN(data.benchmark, data.group)$p.val[,1]
  numMethods <- length(MethodsCompare) + 1
  catplots <- vector("list", numMethods)
  for (i in 1:(numMethods - 1)) {
    catplots[[i]] <- CATplot(pip.norm.DE(data.test, data.group, MethodsCompare[i],
                                         RUV_method = RUV_method, QN_filter = QN_filter,
                                         DE.method = DE.method, Pval = Pval)$p.val[,1],
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

  ggplot(data_cat, aes(x, y, color = curve)) +
    geom_line() +
    ylab("Rate of Agreement with Benchmark") +
    xlab("Significance Rank") +
    ggtitle("CATplot") +
    theme_bw()
}
