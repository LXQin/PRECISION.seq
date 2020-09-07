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
  dat.DE.frame <- data.frame(dm = DEA.res$p.val[,2],
                             p.value = DEA.res$p.val[,1])
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
#' @param MethodsCompare the vector of methods that researchers would like to compare with, and selected from  \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param DE.method the DEA method selected from \code{DE.voom} and \code{DE.edgeR}, and default to be \code{DE.voom}.
#' @param RUV_method the exact RUV method used from \code{RUVg}, \code{RUVr} and \code{RUVs} if \code{norm.RUV} is used.
#' @param QN_filter whether the filtering is performed if \code{norm.QN} is used.
#' @param Pval p-value cut-off point for identifying DE genes, default to be 0.01.
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
#' @param MethodsCompare the vector of methods that researchers would like to compare with, and selected from  \code{norm.none}, \code{norm.TMM}, \code{norm.TC}, \code{norm.UQ}, \code{norm.med}, \code{norm.DESeq}, \code{norm.PoissonSeq}, \code{norm.QN}, \code{norm.SVA}, and \code{norm.RUV}.
#' @param RUV_method the exact RUV method used from \code{RUVg}, \code{RUVr} and \code{RUVs} if \code{method = norm.RUV}
#' @param QN_filter whether the filtering is performed if \code{method = norm.QN}.
#' @param DE.method the method for differential expression analysis from \code{DE.voom} and \code{DE.edgeR}, default to be \code{DE.voom}.
#' @param Pval p-value for identifying DE genes, default to be 0.01
#' @param truth vector of genes that are truly differential expressed
#' @param marker_selection whether selecting a subset of genes/markers for this analysis
#' @param selected_marker vector of a subset of genes/markers for this analysis
#'
#'
#' @return plotly figure for selection of normalization methods
#'
#' @import plotly
#'
#' @export
#'
#' @examples
#' truthgene <- DE.voom(data.benchmark, data.group)$id.list
#' fig.FDR_FNR(data.test, data.group, MethodsCompare = c("norm.none", "norm.TMM", "norm.SVA", "norm.TC"), truth = truthgene)
fig.FDR_FNR <- function(raw, groups, MethodsCompare,
                        RUV_method = NULL, QN_filter = FALSE,
                        DE.method = "DE.voom", Pval = 0.01, truth,
                        marker_selection = FALSE, selected_marker = NULL) {
  numMethods <- length(MethodsCompare)
  FNR <- FDR <- c()
  for (i in 1:numMethods){
    temp <- pip.statistics(raw = raw, groups = groups, norm.method = MethodsCompare[i],
                           RUV_method = RUV_method, QN_filter = QN_filter,
                           DE.method = DE.method, Pval = Pval, truth = truth,
                           marker_selection = marker_selection, selected_marker = selected_marker)
    FNR <- c(FNR, temp$FNR)
    FDR <- c(FDR, temp$FDR)
  }
  stat <- data.frame(FNR = FNR, FDR = FDR, names = MethodsCompare)
  return(plot_ly(data = stat, x = ~FNR, y = ~FDR,
                 text = ~paste('Method: ', names), type = "scatter",  mode = "markers"))
}
