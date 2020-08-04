#' Volcano Figure
#'
#' Function for generating volcano figures with log fold change as x-axis and -log10(p-value) as y-axis.
#'
#' @param normtest normalized data using the method which the researchers desire to compare with the other methods, or it could be the raw data if the adjusting parameters available.
#' @param groups vector of characters indicating the group for each sample.
#' @param adjust adjusting parameters for each sample.
#'
#' @return
#' @export
#'
#' @examples
#' fig.volcano(benchmark, group)
fig.volcano <- function(normtest, groups, adjust = NULL){
  dat.voom <- ifelse(is.null(adjust),
                      DE.voom(RC = normtest, groups = groups, P = 0.01),
                      DE.voom(RC = normtest, groups = groups, P = 0.01, adjust = adjust))
  DE.list <- dat.voom$id.list
  dat.voom.frame <- data.frame(dm = dat.voom$p.val[,2],
                               p.value = dat.voom$p.val[,1])
  mask <- with(dat.voom.frame, p.value < .01)
  cols <- ifelse(mask, "red", "black")

  with(dat.voom.frame, plot(dm, -log10(p.value), cex = .5, pch = 16,
                            col = cols, xlim = c(-3.6, 3.6),
                            ylim = c(0, 6),
                            xlab = "Mean Difference: PMFH - MXF",
                            main = "Volcano Plot"))
  abline(h = 2, lty = 2)
}


#' Relative Log Expression Plot
#'
#' Function for generating relative log expression plot using the normalized test data as the input.
#' The order of the columns in the normalized test data is transformed to be as same as the sorted benchmark data, which makes the figure comparable to the other methods.
#'
#' @param normtest normalized testing using the method which the researchers desire to compare with the other methods
#' @param main the title of the figure
#'
#' @return boxplot for relative log expression
#' @export
#'
#' @examples
#' fig.RLE(data.test, "test without normalization")
fig.RLE = function(normtest, main) {
  benchmark_sort <- c("MXF874", "MXF623", "MXF663", "MXF662", "MXF852", "MXF2516", "MXF776", "MXF832", "MXF771", "MXF847",
                      "MXF632", "MXF911", "MXF912", "MXF908", "MXF848", "MXF861", "MXF863", "MXF856", "MXF906", "MXF866",
                      "MXF836", "MXF879", "MXF875", "MXF772", "MXF742", "MXF7671", "MXF871", "PMFH814", "PMFH837", "PMFH877",
                      "PMFH760", "PMFH898", "PMFH30", "PMFH833", "PMFH6", "PMFH859", "PMFH766", "PMFH817", "PMFH2", "PMFH759",
                      "PMFH816", "PMFH762", "PMFH860", "PMFH41", "PMFH765", "PMFH774", "PMFH857", "PMFH865", "PMFH870", "PMFH878",
                      "PMFH899", "PMFH763", "PMFH904", "PMFH910")
  raw.log <- log2(normtest[, benchmark_sort] + 1)
  rle <- t(apply(raw.log, 1, function(x) x - median(x)))
  boxplot(rle, col = c(rep(rainbow(2)[1], 27), rep(rainbow(2)[2], 27)), ylab = "RLE", ylim = c(-6, 6),
          outline = FALSE, xaxt = "n", main = main)
  legend("topright",c("MXF", "PMFH"), bty = "n",
         pch = "x", cex = 1, col = c(rainbow(2)[1], rainbow(2)[2]))
}


#' Concordance At The Top Plot
#'
#' Function for generating concordance at the top plot, which compares concordance of the p-values obtained
#' from benchmark data without normalization and normalized test data.
#'
#' @param MethodType the type of method that researchers would like to compare with, either "Scale" or "Regression"
#' @param MethodName the name of the new method which will be presented in the legend
#' @param pvalues a vector of p values computed from DE analysis using the normalized data
#'
#' @return
#'
#' @import ffpe
#' @export
#'
#' @examples
#' t <-  runif(1033)
#' names(t) <-  rownames(data.test)
#' fig.CAT("Scale", "Example", t)
fig.CAT <- function(MethodType, MethodName, pvalues){
  benchmark.p <-  DE.voom(data.benchmark, c(rep(c(rep('MXF',9),rep('PMFH',9)),3)))$p.val[,1]
  if(MethodType == "Scale"){
     catplots <-  vector("list", 8)
     catplots[[1]] <- CATplot(subp.test$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[2]] <- CATplot(subp.TMM$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[3]] <- CATplot(subp.TC$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[4]] <- CATplot(subp.UQ$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[5]] <- CATplot(subp.med$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[6]] <- CATplot(subp.DESeq$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[7]] <- CATplot(subp.PoissonSeq$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[8]] <- CATplot(pvalues, benchmark.p, maxrank = 100, make.plot = F)

  color1 <-  rainbow(8)
  par(pty = "s", mar = c(2, 5, 2, 10))
  plot(catplots[[1]], ylim = c(0, 0.7), col = color1[1],
       lwd = 2, type = "l", ylab = "Rate of Agreement with Benchmark",
       xlab = "Significance Rank of MiRNAs",
       main = "CATplot: Scaling-based Normalization Methods")
  lines(catplots[[2]], col = color1[2], lwd = 2)
  lines(catplots[[3]], col = color1[3], lwd = 2)
  lines(catplots[[4]], col = color1[4], lwd = 2)
  lines(catplots[[5]], col = color1[5], lwd = 2)
  lines(catplots[[6]], col = color1[6], lwd = 2)
  lines(catplots[[7]], col = color1[7], lwd = 2)
  lines(catplots[[8]], col = color1[8], lwd = 2)
  legend(x = 110, y = 0.72, xpd = TRUE,
         legend = c("No Norm", "TMM", "TC", "UQ", "Med", "DESeq", "PoissonSeq", MethodName),
         col = color1, lwd = 2)
  }else if (MethodType == "Regression"){
     catplots <-  vector("list", 8)
     catplots[[1]] <- CATplot(subp.test$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[2]] <- CATplot(subp.QN$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[3]] <- CATplot(subp.QN.filter$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[4]] <- CATplot(subp.SVA$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[5]] <- CATplot(subp.RUVg$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[6]] <- CATplot(subp.RUVr$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[7]] <- CATplot(subp.RUVs$comp.test, benchmark.p, maxrank = 100, make.plot = F)
     catplots[[8]] <- CATplot(pvalues, benchmark.p, maxrank = 100, make.plot = F)

  color2 = rainbow(8)
  par(pty = "s", mar = c(2, 5, 2, 10))
  plot(catplots[[1]], ylim = c(0, 0.7), col = color1[1],
       lwd = 2, type = "l",ylab = "Rate of Agreement with Benchmark",
       xlab = "Significance Rank of MiRNAs",
       main = "CATplot: Regression-based Normalization Methods")
  lines(catplots[[2]], col = color2[2], lwd = 2, lty = 3)
  lines(catplots[[3]], col = color2[3], lwd = 2, lty = 3)
  lines(catplots[[4]], col = color2[4], lwd = 2, lty = 3)
  lines(catplots[[5]], col = color2[5], lwd = 2, lty = 3)
  lines(catplots[[6]], col = color2[6], lwd = 2, lty = 3)
  lines(catplots[[7]], col = color2[7], lwd = 2, lty = 3)
  lines(catplots[[8]], col = color1[8], lwd = 2, lty = 3)
  legend(x = 110, y = 0.72, xpd = TRUE,
         legend = c("No Norm", "QN", "QN with Filtering", "SVA", "RUVg", "RUVr", "RUVs", MethodName),
         col = c(color1[1], color2[2:8]), lwd = 2, lty=c(1, rep(3,6)))
  }
}


