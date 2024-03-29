#' Volcano Figure
#'
#' Function for generating volcano plot with log fold change as x-axis and
#' -log10(p-value) as y-axis. Volcano plots are typically used for evaluation
#' of differential expression statuses in a data set.
#'
#' @param DEA.res Results of a differential expression analysis using the
#' included methods \code{\link{DE.voom}} or \code{\link{DE.edgeR}}.
#' Must be a list including
#' \describe{
#'   \item{p.val}{p-values for differential expression}
#'   \item{log2.FC}{Fold change on a log2 scale}
#' }
#' @param Pval optional p-value cutoff for differential expression
#' @param title optional Figure title
#'
#' @return volcano plot (ggplot2 object)
#' @export
#'
#' @examples
#' voom.benchmark <- DE.voom(data.benchmark, data.group)
#' fig.volcano(voom.benchmark, title = "Benchmark data (differential expression)")
fig.volcano <- function(DEA.res, Pval=0.01, title){
  if(missing(title)) {
    title <- ""
  }
  dat.DE.frame <- data.frame(dm = DEA.res$log2.FC,
                             p.value = DEA.res$p.val)
  mask <- with(dat.DE.frame, p.value < Pval)
  cols <- ifelse(mask, "red", "black")

  xlim = round(1.15 * max(abs(dat.DE.frame$dm)), 1)
  ylim = 1.15 * round(max(-log10(dat.DE.frame$p.value)), 1)

  p <- ggplot(aes(x = dm, y = -log10(p.value)), data = dat.DE.frame) +
    geom_point(color = cols) +
    geom_hline(yintercept = -log10(Pval), lty = 2) +
    xlim(-xlim, xlim) +
    ylim(0, ylim) +
    xlab("Mean Difference") +
    ggtitle(title) +
    theme_bw()

  return(p)
}


#' Relative Log Expression Plot
#'
#' Function for generating relative log expression plots based on raw or
#' normalized count data as the input.
#'
#' @param data Raw or normalized count data
#' @param groups vector of characters indicating the group for each sample
#' (only 2 groups are allowed).
#' @param title optional Figure title
#'
#' @return ggplot boxplot for relative log expression
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' fig.RLE(data.test, data.group, title="test data (no normalization)")
fig.RLE = function(data, groups, title) {
  if(missing(title)) {
    title <- ""
  }
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
#' Function for generating a concordance at the top plot, which compares
#' concordance of the p-values of differential expression --- typically from the
#' normalized or raw data set \code{\link{data.test}} --- to the p-values of
#' differential expression of a gold standard --- typically
#' the data set \code{\link{data.benchmark}}.
#'
#'
#' @param DEA Results of a differential expression analysis using the
#' included methods \code{\link{DE.voom}} or \code{\link{DE.edgeR}}.
#' Must be a list including
#' \describe{
#'   \item{id.list}{List of differentially expressed markers}
#'   \item{p.val}{p-values for differential expression}
#' }
#' @param truth.DEA Gold standard (assumed truth) for differential expression.
#' Must be in the same format as \code{DEA}.
#' @param title optional Figure title
#' @param maxrank optional specify the maximum size of top-ranked items that you want to plot.
#' @param subset optional vector of a subset of genes/markers for this analysis
#
#' @return figure of concordance for comparison
#'
#' @import ggplot2
#' @import ffpe
#'
#' @export
#'
#' @examples
#' voom.benchmark <- DE.voom(data.benchmark, data.group)
#' test.norm <- pip.norm(raw=data.test, groups=data.group, norm.method = "all")
#' test.DE <- list(
#' TMM = DE.voom(RC=test.norm$TMM$dat.normed, groups = data.group),
#' TC = DE.voom(RC=test.norm$TC$dat.normed, groups = data.group),
#' UQ = DE.voom(RC=test.norm$UQ$dat.normed, groups = data.group),
#' med = DE.voom(RC=test.norm$med$dat.normed, groups = data.group),
#' DESeq = DE.voom(RC=test.norm$DESeq$dat.normed, groups = data.group),
#' PoissonSeq = DE.voom(RC=test.norm$PoissonSeq$dat.normed, groups = data.group),
#' QN = DE.voom(RC=test.norm$QN$dat.normed, groups = data.group),
#' RUVg = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVg$adjust.factor),
#' RUVs = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVs$adjust.factor),
#' RUVr = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVr$adjust.factor),
#' SVA = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$SVA$adjust.factor),
#' noNorm = DE.voom(RC=data.test, groups = data.group))
#'
#' fig.CAT(DEA = test.DE, truth.DEA = voom.benchmark, title = "Example of CAT plot")
fig.CAT <- function(DEA, truth.DEA, title, maxrank=100, subset=NULL){
  if(missing(title)) {
    title <- ""
  }
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


#' Selection of normalization methods based on golden standards (FDR and FNR)
#'
#' Function for boxplot of FDR and FNR based on the DEA from the golden benchmark (typically
#' the data set \code{\link{data.benchmark}}) and the normalized test dataset (typically the
#' normalized or raw data set \code{\link{data.test}})
#'
#' @param DEA a list of differential expression analysis results with the element names to be the normalization methods
#' @param truth.DEA differential expression analysis results from the benchmark (gold standard) obtained from DE.voom, DE.edge, or any method from the users by storing the results as same as DE methods in the package (including DE genes, p-values and log2 fold changes)
#' @param title Figure title
#' @param subset optional Vector of a subset of markers.
#' If given, the plot will be limited to the given subset
#' of markers. Leave \code{NULL} if all markers should be considered.
#'
#' @return figure for selection of normalization methods
#'
#' @import ggrepel
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' voom.benchmark <- DE.voom(data.benchmark, data.group)
#' test.norm <- pip.norm(raw=data.test, groups=data.group, norm.method = "all")
#' test.DE <- list(
#' TMM = DE.voom(RC=test.norm$TMM$dat.normed, groups = data.group),
#' TC = DE.voom(RC=test.norm$TC$dat.normed, groups = data.group),
#' UQ = DE.voom(RC=test.norm$UQ$dat.normed, groups = data.group),
#' med = DE.voom(RC=test.norm$med$dat.normed, groups = data.group),
#' DESeq = DE.voom(RC=test.norm$DESeq$dat.normed, groups = data.group),
#' PoissonSeq = DE.voom(RC=test.norm$PoissonSeq$dat.normed, groups = data.group),
#' QN = DE.voom(RC=test.norm$QN$dat.normed, groups = data.group),
#' RUVg = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVg$adjust.factor),
#' RUVs = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVs$adjust.factor),
#' RUVr = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVr$adjust.factor),
#' SVA = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$SVA$adjust.factor),
#' noNorm = DE.voom(RC=data.test, groups = data.group))
#'
#' fig.FDR_FNR(DEA = test.DE, truth.DEA = voom.benchmark, title = "Example of FDR FNR scatterplot")
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
    ggtitle(title)+
    scale_x_continuous(labels = scales::percent, limits=c(0,1),
                       breaks = scales::pretty_breaks(n = 4)) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1))


  return(p)
}


#' Venn diagram for p-values
#'
#' Venn diagram is used to identify the performance of different normalization methods based on intersection of differential expressed genes.
#' @param truth.DEA.pval p-values of differential expression analysis results from the benchmark (gold standard) obtained from DE.voom, DE.edge, or any method from the users by storing the results as same as DE methods in the package (including DE genes, p-values and log2 fold changes)
#' @param DEA.res.pval p-values as a result from prior differential expression
#' analysis, e.g. using \code{\link{DE.voom}} or \code{\link{DE.edgeR}}.
#' @param Pvalue Cut-off point for p-values for identifying significant
#' differential expression.
#' @param title optional, figure title
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
#' fig.venn(truth.DEA = benchmark.voom$p.val, DEA.res = test.voom$p.val, Pvalue = 0.01)
fig.venn <- function(truth.DEA.pval, DEA.res.pval, Pvalue, title){
  if(missing(title)) {
    title <- ""
  }
  bench.sig <- (truth.DEA.pval < Pvalue)
  test.sig <- (DEA.res.pval[names(truth.DEA.pval)] < Pvalue)
  venn <- vennCounts(cbind(bench.sig, test.sig))
  p <- as.ggplot(function() vennDiagram(venn,
                                        names = c("Benchmark", "Test"),
                                        cex = 1.5, counts.col = rainbow(1))) +
    ggtitle(title)

  return(p)
}


#' Dendrogram for clustering p-values
#'
#' Function for clustering normalization methods based on the p-values pattern calculated from the same dataset.
#'
#' @param DEA.pval.list A list of p-values from differential expression analysis results with the element names to be the normalization methods
#' @param title optional Figure title
#' @param subset optional Vector of a subset of markers.
#' If given, the dendrogram analysis will be limited to the given subset
#' of markers. Leave \code{NULL} if all markers should be considered.
#' @return Figure of dendrogram
#'
#' @import ggdendro
#' @import Biobase
#' @export
#'
#' @examples
#'
#' test.norm <- pip.norm(raw=data.test, groups=data.group, norm.method = "all")
#' test.DE <- list(
#' TMM = DE.voom(RC=test.norm$TMM$dat.normed, groups = data.group),
#' TC = DE.voom(RC=test.norm$TC$dat.normed, groups = data.group),
#' UQ = DE.voom(RC=test.norm$UQ$dat.normed, groups = data.group),
#' med = DE.voom(RC=test.norm$med$dat.normed, groups = data.group),
#' DESeq = DE.voom(RC=test.norm$DESeq$dat.normed, groups = data.group),
#' PoissonSeq = DE.voom(RC=test.norm$PoissonSeq$dat.normed, groups = data.group),
#' QN = DE.voom(RC=test.norm$QN$dat.normed, groups = data.group),
#' RUVg = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVg$adjust.factor),
#' RUVs = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVs$adjust.factor),
#' RUVr = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$RUVr$adjust.factor),
#' SVA = DE.voom(RC=data.test, groups = data.group, normalized=FALSE, adjust=test.norm$SVA$adjust.factor),
#' noNorm = DE.voom(RC=data.test, groups = data.group))
#' test.DE.pval <- lapply(1:12, function(x) test.DE[[x]]$p.val)
#' names(test.DE.pval) <- names(test.DE)
#'
#' fig.dendrogram(DEA.pval.list = test.DE.pval, title = "Example of dendrogram")
fig.dendrogram <- function(DEA.pval.list, title, subset=NULL){
  if(missing(title)) {
    title <- ""
  }
  if(is.null(subset)){
    genes <- names(DEA.pval.list[[1]])
  } else{
    genes <- subset
  }
  pval_frame <- data.frame(sapply(DEA.pval.list, function(x) x[genes]))
  hc <- hclust(dist(t(-log10(pval_frame))), method = "ward.D")
  p.dendro <- ggdendrogram(hc, rotate = FALSE, size = 2) + ggtitle(title)
  rm(genes, pval_frame, hc)
  return(p.dendro)
}


#' Boxplot of FDR and FNR for Simulated data
#'
#' Function for the boxplots of FDR and FNR based on golden truth (typically the simulated benchmark data set \code{\link{simulated.data()}})
#' and the differential expression analysis (typically the results from the simulated normalized test data set \code{\link{simulated.data()}}).
#'
#' @param DEA.list A list with each element as a normalization method including the differential expression analysis (DEA) results from multiple simulated data
#' @param truth.DEA.list A list of DEA results of benchmark data for multiple simulated dataset
#' @param title The title of the figure
#' @param subset optional A vector of a subset of genes/markers for this analysis
#'
#' @import ggplot2
#' @return
#' @export
#'
#' @examples
#' simulated <- simulated.data(proportion = c(0.15, 0.25),  median = c(2, 4), numsets = 10)
#' norm.methods <- c("norm.TMM", "norm.TC", "norm.UQ", "norm.med")
#' test.DEA.list <- list()
#' for (i in norm.methods){
#' temp = lapply(1:10, function(x) pip.norm.DE(raw = simulated[[x]]$simulated_test,
#'                                             groups = simulated[[x]]$simulated_group, norm.method = i))
#' test.DEA.list <- append(test.DEA.list, list(temp))
#' }
#' names(test.DEA.list) = norm.methods
#' benchmark.DEA.list <- lapply(1:10, function(x) DE.voom(RC = simulated[[x]]$simulated_benchmark,
#'                                                        groups = simulated[[x]]$simulated_group, P = 0.01))
#' fig.FDR_FNR.boxplot(DEA.list = test.DEA.list, truth.DEA.list = benchmark.DEA.list, title = "Example of FDR and FNR boxplot")
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
