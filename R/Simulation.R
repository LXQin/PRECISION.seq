
#' The Algorithm for Obtaining the Simulation Plan
#' This algorithm is tailored for the data.benchmark included in this package.
#'
#' @param benchmark the data.benchmark used for the simulation
#'
#' @return
#' @export
#'
#' @examples simulation.plan <- simulation.algorithm(data.benchmark)
#' @keywords internal
simulation.algorithm <- function(benchmark){
  set.seed(12345)
  colnames(benchmark) = sub(".*_", "", colnames(benchmark))

  # Cluster 54 samples into 2 clusters
 #cluster = pam(t(benchmark), 2)$clustering
  cluster1.D = c("D9", "D7", "D17", "D19", "D13", "D8", "D20")
  cluster1.E = c("E17", "E5", "E7", "E3", "E4", "E19", "E2")
  cluster1.C = c("C17", "C2", "C19", "C5", "C18", "C10")
  cluster2.D = c("D1", "D3", "D10", "D18", "D4", "D12", "D15", "D6", "D14", "D11", "D2")
  cluster2.E = c("E12", "E13", "E8", "E18", "E20", "E11", "E16", "E6",  "E10", "E1", "E15")
  cluster2.C = c("C9", "C1", "C15", "C3", "C16", "C7", "C8", "C13", "C6", "C20", "C14", "C4")

  # Obtain random seeds
  random.seed = matrix(, nrow = 50, ncol = 18)
  colnames(random.seed) = c(rep('cluster1', 9), rep('cluster2', 9))
  for (i in 1:50) {
    random.seed[i, 1:3] = sample(cluster1.D, 3)
    random.seed[i, 4:6] = sample(cluster1.E, 3)
    random.seed[i, 7:9] = sample(cluster1.C, 3)
    random.seed[i, 10:12] = sample(cluster2.D, 3)
    random.seed[i, 13:15] = sample(cluster2.E, 3)
    random.seed[i, 16:18] = sample(cluster2.C, 3)
  }

  # Allocate remaining samples
  remain.locate = matrix(ncol = 36)
  for (i in 1:50) {
    remain.locate.temp = matrix(, nrow = 400, ncol = 36)
    candidate.D = subset(colnames(benchmark)[1:18], !colnames(benchmark)[1:18] %in% random.seed[i,])
    candidate.E = subset(colnames(benchmark)[19:36], !colnames(benchmark)[19:36] %in% random.seed[i,])
    candidate.C = subset(colnames(benchmark)[37:54], !colnames(benchmark)[37:54] %in% random.seed[i,])
    remain.locate.temp[1:200, 1:6] = t(replicate(200, sample(candidate.D, 6)))
    remain.locate.temp[1:200, 7:12] = t(replicate(200, sample(candidate.E, 6)))
    remain.locate.temp[1:200, 13:18] = t(replicate(200, sample(candidate.C, 6)))
    remain.locate.temp[1:200, 19:24] = t(apply(remain.locate.temp[1:200, 1:6], 1,
                                               function(x) subset(candidate.D, !candidate.D %in% x)))
    remain.locate.temp[1:200, 25:30] = t(apply(remain.locate.temp[1:200, 7:12], 1,
                                               function(x) subset(candidate.E, !candidate.E %in% x)))
    remain.locate.temp[1:200, 31:36] = t(apply(remain.locate.temp[1:200, 13:18], 1,
                                               function(x) subset(candidate.C, !candidate.C %in% x)))
    remain.locate.temp[201:400, 1:18] = remain.locate.temp[1:200, 19:36]
    remain.locate.temp[201:400, 19:36] = remain.locate.temp[1:200, 1:18]
    remain.locate = rbind(remain.locate, remain.locate.temp)
  }
  remain.locate = na.omit(remain.locate)

  # Obtain complete simulation plan
  random.seed.full = random.seed[rep(seq_len(nrow(random.seed)), each = 400),]
  simulate.plan = cbind(random.seed.full[,1:9], remain.locate[, 1:18],  random.seed.full[,10:18], remain.locate[, 19:36])
  colnames(simulate.plan) = c(rep('Cluster_1',27),rep('Cluster_2',27))

  return(simulate.plan)
}
