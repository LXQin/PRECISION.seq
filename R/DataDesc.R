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


#' Simulated Data Generation
#'
#' Function for generating simulated datasets according to different proportions of DE and medians of mean difference.
#' The range of DE proportion availiable is between 0 and 0.387, and the range of median of mean difference is between -2.52 and 4.59.
#'
#' @param proportion_L the lowest proportion of DE for simulated dataset filtering
#' @param proportion_R the highest proportion of DE for simulated dataset filtering
#' @param median_L the lowest median of mean difference for simulated dataset filtering
#' @param median_R the highest median of mean difference for simulated dataset filtering
#' @param numsets number of simulated datasets requested for returning (randomly drawn from the available sets). If it exceeds the maximum of availiable datasets, all the availiable sets will be returned.
#'
#' @return list containing list of simulated benchmark data and test data
#' @import magrittr
#' @import dplyr
#' @export
#'
#' @examples
#' simulated <- simu(0.0175, 0.0225, -0.5, 0.5, 10)
simu <- function(proportion_L, proportion_R, median_L, median_R, numsets){
  benchmark_simu <- data.benchmark
  test_simu <- data.test
  colnames(benchmark_simu) <- sub(".*_", "", colnames(benchmark_simu))
  colnames(test_simu) <- colnames(benchmark_simu)
  s <- data.simulation %>%
    filter(proportion > proportion_L & proportion < proportion_R) %>%
    filter(median > median_L & median < median_R)
  rowselect <- if(nrow(s) > numsets){sample(nrow(s), numsets)}else{1:nrow(s)}
  s <- as.matrix(s[rowselect, 1:54])
  benchmark_simued <- test_simued <- list()
  for (i in 1:nrow(s)) {
    benchmark_simued[[i]] <- benchmark_simu[, s[i,]]
    test_simued[[i]] <- test_simu[, s[i,]]
  }
  return(list(simulated_benchmark = benchmark_simued,
              simulated_test = test_simued))
}

