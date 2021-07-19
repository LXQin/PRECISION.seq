#' MiRNA Sequencing Benchmark Data
#'
#' Myxofibrosarcoma (MXF) and pleomorphic malignant fibrous histiocytoma (PMFH) are the two most common and aggressive subtypes of genetically complex soft tissue sarcoma.
#' This dataset includes three libraries used for sequencing 54 individual tumor samples. Library preparation and read capture were each processed by a single experienced technician in one run.
#'
#'
#' @format A data frame with 1033 rows and 54 columns. Here are the examples of the column and row naming rule:
#' \describe{
#'   \item{MXF2516_D13}{One sample belonging to MXF and library D with sample ID MXF2516 and library ID D13.}
#'   \item{hsa-let-7a-2*}{Gene ID.}
#' }
"data.benchmark"

#' MiRNA Sequencing Test Data
#'
#' Myxofibrosarcoma (MXF) and pleomorphic malignant fibrous histiocytoma (PMFH) are the two most common and aggressive subtypes of genetically complex soft tissue sarcoma.
#' MiRNAs for the same 54 tumors used for the benchmark data were re-sequenced using neither uniform handling nor balanced library assignment.
#' In this study these samples were sequenced in the order of sample collection and processed in multiple runs.
#'
#'
#' @format A dataframe with 1033 rows and 54 columns. Here are the examples of the column and row naming rule:
#' \describe{
#'   \item{MXF2516}{One sample belonging to MXF sample ID MXF2516.}
#'   \item{hsa-let-7a-2*}{Gene ID.}
#' }
"data.test"

#' MiRNA Information
#'
#' @format Information of the microRNAs in test and benchmark dataset, including their names, exact sequence, length and GC content.
"data.miR.info"

#' Sample Labels of the test and benchmark data
#'
#' @format A set of labels for the test and benchmark samples, either "MXF" (myxofibrosarcoma) or "PMFH" (pleomorphic malignant fibrous histiocytoma).
"data.group"

#' Simulation Plan
#'
#' @format A dataframe including information of allocation about 20,000 simulation datasets, which contains the proportion of differential expression and median of mean difference for each dataset.
"data.simulation"


#' Simulated Data
#'
#' Function for providing simulated datasets according to different proportions of DE and medians of mean difference.
#' The range of DE proportion available is between 0 and 0.387, and the range of median of mean difference is between -2.52 and 4.59.
#'
#' @param proportion the range of proportion of DE for simulated dataset filtering
#' @param median the range of median of mean difference for simulated dataset filtering
#' @param median_R the highest median of mean difference for simulated dataset filtering
#' @param numsets number of simulated datasets requested for returning (randomly drawn from the available sets). If it exceeds the maximum of availiable datasets, all the availiable sets will be returned.
#'
#' @return list containing list of simulated benchmark data and test data
#' @import magrittr
#' @import dplyr
#' @export
#'
#' @examples
#' simulated <- simulated.data(0.0175, 0.0225, -0.5, 0.5, 10)
simulated.data <- function(proportion, median, numsets){

  proportion_L <- proportion[1]
  proportion_R <- proportion[2]
  median_L <- median[1]
  median_R <- median[2]

  if (proportion_L <=0 | proportion_R >= 0.387 | median_L <= -2.52 | median_R >= 4.59) {
    stop("The range exceeds the availability.")
  }

  benchmark_simu <- data.benchmark
  test_simu <- data.test
  colnames(benchmark_simu) <- sub(".*_", "", colnames(benchmark_simu))
  colnames(test_simu) <- colnames(benchmark_simu)
  s <- data.simulation %>%
    filter(proportion > proportion_L & proportion < proportion_R) %>%
    filter(median > median_L & median < median_R)
  rowselect <- if(nrow(s) > numsets){
    sample(nrow(s), numsets)
  }else{
    1:nrow(s)
    warning("Number of datasets that can be provided is less than the demand.")
  }

  s <- as.matrix(s[rowselect, 1:54])
  data.simulated <- list()
  for (i in 1:nrow(s)) {
    simulated_benchmark <-  benchmark_simu[, s[i,]]
    simulated_test <-  test_simu[, s[i,]]
    data.simulated[[i]] <- list(simulated_benchmark = simulated_benchmark,
                                simulated_test = simulated_test)
  }
  return(data.simulated)
}

