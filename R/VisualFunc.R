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

  p <- ggplot(aes(x = dm, y = -log10(p.value)), data = dat.DE.frame) +
    geom_point(color = cols) +
    geom_hline(yintercept = 2, lty = 2) +
    xlim(-xlim, xlim) +
    ylim(0, ylim) +
    xlab("Mean Difference") +
    ggtitle(title) +
    theme_bw()

  return(p)
}


#' Relative Log Expression Plot
#'
#' Function for generating relative log expression plot using the read count data as the input.
#'
#' @param data normalized testing using the method which the researchers desire to compare with the other methods
#' @param groups vector of characters indicating the group for each sample (only 2 groups are allowed).
#' @param title the title of the figure
#'
#' @return ggplot boxplot for relative log expression
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' fig.RLE(data.test, data.group, "test without normalization")
fig.RLE = function(data, groups, title) {
  data <- ifelse(data<0, 0, data)  # Catch negative read counts and set them to 0
  raw.log <- log2(data + 1)
  rle <- t(apply(raw.log, 1, function(x) x - median(x)))

  ylim = round(max(1.5*(apply(rle, 2, IQR)) + max(abs(apply(rle, 2, function(x) quantile(x, probs=c(0.25, 0.75), na.rm=TRUE))))), 1)
  df <- tidyr::gather(as.data.frame(rle), Sample, RLE)
  df$group <- rep(groups, each = nrow(rle))
  df$Sample <- factor(df$Sample, levels = colnames(rle))
  p <- ggplot(data = df, aes(x = Sample, y = RLE, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    ylim(-ylim, ylim) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="bottom") +
    ggtitle(title)
  return(p)
 }


#' Concordance At The Top Plot
#'
#' Function for generating concordance at the top plot, which compares concordance of the p-values obtained
#' from benchmark data without normalization and normalized test data.
#'
#'
#' @param DEA a list of differential expression analysis results with the element names to be the normalization methods
#' @param truth.DEA differential expression analysis results from the benchmark (gold standard) obtained from DE.voom, DE.edge, or any method from the users by storing the results as same as DE methods in the package (including DE genes, p-values and log2 fold changes)
#' @param title figure title
#' @param maxrank optionally specify the maximum size of top-ranked items that you want to plot.
#' @param subset vector of a subset of genes/markers for this analysis
#
#' @return figure of concordance for comparison
#'
#' @import ggplot2
#' @import ffpe
#'
#' @export
#'
#' @examples
fig.CAT <- function(DEA, truth.DEA, title, maxrank=100, subset=NULL){
  # Subset DEAs to a given set of miRNAs
  if (!is.null(subset)) {
    DEA <- DEA[subset, ]
    truth.DEA <- truth.DEA[subset, ]
  }
  # Reduce Data to named p.values
  truth <- truth.DEA$p.val
  names(truth) <- truth.DEA$id.list
  compare <- list()
  for(i in 1:length(DEA)) {
    compare <- append(compare, list(DEA[[i]]$p.val))
    names(compare[[i]]) <- DEA[[i]]$id.list
  }
  names(compare) <- names(DEA)

  catplots <- list()
  for(i in 1:length(compare)) {
    catplots <- append(catplots,
                       list(CATplot(compare[[i]], truth,
                                    maxrank = maxrank, make.plot=FALSE)))
  }
  names(catplots) <- names(compare)

  data_cat <- data.frame(x = numeric(), y = numeric(), curve = character())
  for(i in 1:length(catplots)){
    y = catplots[[i]][,2]
    x = 1:length(y)
    data_cat = rbind(data_cat, data.frame(x, y, curve = names(catplots)[i]))
  }
  data_cat$curve = factor(data_cat$curve)

  p <- ggplot(data_cat, aes(x, y, color = curve, linetype = curve)) +
    theme(legend.title = element_blank()) +
    geom_line(size=.75) +
    ylab("Rate of Agreement with Benchmark") +
    xlab("Significance Rank") +
    theme(legend.title=element_blank()) +
    ggtitle(title) +
    ylim(c(0,1)) +
    theme_bw()

  return(p)
}


#' Selection of normalization methods based on golden standards
#'
#' @param DEA a list of differential expression analysis results with the element names to be the normalization methods
#' @param truth.DEA differential expression analysis results from the benchmark (gold standard) obtained from DE.voom, DE.edge, or any method from the users by storing the results as same as DE methods in the package (including DE genes, p-values and log2 fold changes)
#' @param title figure title
#' @param subset vector of a subset of genes/markers for this analysis
#'
#' @return figure for selection of normalization methods
#'
#' @import ggrepel
#' @import ggplot2
#'
#' @export
#'
#' @examples
fig.FDR_FNR <- function(DEA, truth.DEA, title, subset=NULL) {
  DE.stats <- list()
  for(DE.res in DEA) {
    DE.stats <-
      append(DE.stats,
             list(DE.statistics(markers=names(truth.DEA$p.val),
                                id.list=DE.res$id.list,
                                truth=truth.DEA$id.list,
                                selected.marker=subset)))
  }
  names(DE.stats) <- names(DEA)
  FNR <- sapply(1:length(names(DE.stats)), function(x)DE.stats[[x]]$FNR)
  FDR <- sapply(1:length(names(DE.stats)), function(x)DE.stats[[x]]$FDR)
  stat <- data.frame(FNR = FNR, FDR = FDR, names = names(DE.stats))
  p <- ggplot(stat, aes(FNR, FDR, label = names)) + geom_point(color = "red") +
    geom_text_repel() + theme_bw() +
    ggtitle(title)

  return(p)
}


#' Venn diagram for p-values
#'
#' Venn diagram is used to identify the performance of different normalization methods based on intersection of differential expressed genes.
#' @param truth.DEA differential expression analysis results from the benchmark (gold standard) obtained from DE.voom, DE.edge, or any method from the users by storing the results as same as DE methods in the package (including DE genes, p-values and log2 fold changes)
#' @param DEA.res p-values as a result from prior differential expression
#' analysis, e.g. using \code{\link{DE.voom}} or \code{\link{DE.edgeR}}.
#' @param Pvalue Cut-off point for p-values for identifying significant
#' differential expression.
#'
#' @return A Venn diagram
#' @import limma
#' @import ggplotify
#'
#' @export
#'
#' @examples
#' benchmark.voom <- DE.voom(RC = data.benchmark, groups = data.group, P = 0.01)
#' test.voom <- DE.voom(RC = data.test, groups = data.group, P = 0.01)
#' fig.venn(benchmark.voom$p.val, test.voom$p.val, 0.01)
fig.venn <- function(truth.DEA, DEA.res, Pvalue){
  bench.sig <- (truth.DEA$p.val < Pvalue)
  test.sig <- (DEA.res$p.val[names(truth.DEA$p.val)] < Pvalue)
  venn <- vennCounts(cbind(bench.sig, test.sig))
  p <- as.ggplot(function() vennDiagram(venn,
                                        names = c("Benchmark", "Test"),
                                        cex = 1.5, counts.col = rainbow(1)))

  return(p)
}


#' Dendrogram for clustering p-values
#'
#' Function for clustering normalization methods based on the p-values pattern calculated from the same dataset.
#'
#' @param DEA a list of differential expression analysis results with the element names to be the normalization methods
#' @param title figure title
#' @param subset vector of a subset of genes/markers for this analysis
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
fig.dendrogram <- function(DEA, title, subset=NULL){
  if(is.null(subset)){
    genes <- names(DEA[[1]]$p.val)
  } else{
    genes <- subset
  }
  pval_frame <- data.frame(sapply(DEA, function(x) x$p.val[genes]))
  hc <- hclust(dist(t(-log10(pval_frame))), method = "ward.D")
  p.dendro <- ggdendrogram(hc, rotate = FALSE, size = 2) + ggtitle(title)
  rm(genes, pval_frame, hc)
  return(p.dendro)
}


#' Boxplot of FDR and FNR for Simulated data
#'
#' @param DEA.list a list with each element as a normalization method including the differential expression analysis (DEA) results from multiple simulated data
#' @param truth.DEA.list a list of DEA results of benchmark data for multiple simulated dataset
#' @param title the title of the figure
#' @param subset a vector of a subset of genes/markers for this analysis
#'
#' @import ggplot2
#' @return
#' @export
#'
#' @examples
fig.FDR_FNR.boxplot <- function(DEA.list, truth.DEA.list, title, subset=NULL){
  DE.stats <- list()
  num.simulated.sets <- length(truth.DEA.list)
  for(DE.res in DEA.list) {
    DE.stats <- append(DE.stats,
                       list(lapply(1:num.simulated.sets,
                                   function(x) DE.statistics(markers=names(truth.DEA.list[[x]]$p.val),
                                                             id.list=DE.res[[x]]$id.list,
                                                             truth=truth.DEA.list[[x]]$id.list,
                                                             selected.marker=subset))))
  }
  names(DE.stats) <- names(DEA.list)

  FNR <- FDR <- list()
  for (element in DE.stats){
    FNR <- c(FNR, list(c(sapply(1:num.simulated.sets, function(x) element[[x]]$FNR))))
    FDR <- c(FDR, list(c(sapply(1:num.simulated.sets, function(x) element[[x]]$FDR))))
  }
  names(FDR) <- names(FNR) <- names(DEA.list)

  FNR_FDR_tab <- data.frame(rbind(reshape2::melt(FDR), reshape2::melt(FNR)))
  FNR_FDR_tab$stat <- c(rep("FDR", nrow(reshape2::melt(FDR))), rep("FNR", nrow(reshape2::melt(FNR))))

  p <- ggplot(data = FNR_FDR_tab , aes(x=L1, y=value)) +
    geom_boxplot(aes(fill=stat))  +
    theme_bw() +
    ggtitle(title) + scale_fill_manual(values = c("white", "grey")) +
    xlab("") + ylab("")

  return(p)
}
