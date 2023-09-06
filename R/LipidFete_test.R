#' @title Calculate p-value of lipidomic features by permutation test
#'
#' @description Calculate p-value of lipidomic features by permutation test.
#' The test of region statistics smoothing with Gaussian kernel integrates
#' neighbor's information and provides a stable testing result under small sample size.
#'
#' @param X A \code{matrix} of lipid expression data. Row is samples, and column is
#' lipid feature such as double bond and chain length.
#' @param X.info A \code{data.frame} or \code{matrix} of lipid feature names. (Only allow 1 or 2 columns)
#' @param group A vector of observed grouping variable. (Only allow 0 and 1)
#' @param radius A value to decide the number of neighbor included
#' @param own.contri A proportion of own contribution (suggest 0.5 ~ 1)
#' @param x.distance A value of unit distance in x-axis
#' @param y.distance A value of unit distance in y-axis
#' @param dimension A value to indicate the dimension of lipid features (1 or 2)
#' @param permute.time A value of permutations
#' @param seed A value of random seed index
#'
#' @details write XXX
#'
#' @return A \code{data.frame} of containing lipidomic feature testing result.
#' \itemize{
#'     \item{X.info: }{lipidomic feature information.}
#'     \item{avg.expr: }{average expression.}
#'     \item{direction: }{the feature expression level in group 1 is higher (+) or lower (-).}
#'     \item{smoothing.pval.BH: }{p-value of lipidomic features by permutation test with BH correction.}
#'     \item{log2.FC: }{log2-transformed fold change.}
#' }
#'
#' @examples
#' data(lipid2D)
#' X <- t(as.matrix(lipid2D[,-c(1:2)]))
#' X.info <- lipid2D[,1:2]
#' group <- rep(c(0, 1), c(52,32))
#' LipidFete.test(X = X,
#'                X.info = X.info,
#'                group = group,
#'                radius = 3,
#'                own.contri = 0.5,
#'                x.distance = 2,
#'                y.distance = 1,
#'                dimension = 2,
#'                permute.time = 100000)
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{region.stat}, \link{region.plot.1D}, \link{region.plot.2D}}
#' @export
LipidFete.test <- function(X, X.info, group, radius = 3, own.contri = 0.5,
                           x.distance = 1, y.distance = 1, dimension = c(1,2),
                           permute.time = 100000, seed = 1234, ...){
  ## step 1: calculate region statistics
  region.stat.obs <- region.stat(X = X,
                                 Y = as.matrix(group))

  ## step 2: calculate distance matrix by euclidean distance
  dist.input <- as.matrix(X.info)
  dist.input[,1] <- dist.input[,1]/x.distance
  if(dimension == 2) dist.input[,2] <- dist.input[,2]/y.distance
  dist.mat <- as.matrix(dist(dist.input, method = "euclidean"))
  dist.mat[dist.mat>=radius] <- Inf   # only consider distance < 3


  ## step 3: smooth by Gaussian kernel
  dist.idx <- which.max(colSums(exp(-dist.mat^2)))
  radius.seq <- seq(0.01, 5, 0.01)
  own.contri.seq <- sapply(radius.seq, function(x) 1/sum(exp(-dist.mat[,dist.idx]^2*x)))
  radius.idx <- sum(own.contri.seq <= own.contri)
  if(radius.idx + 1 > length(radius.seq)){
    stop(paste0("Based on the data structure and 'radius' setting, \n",
                "       the 'own.contri' should smaller than ",
                round(max(own.contri.seq), 4)))
  }
  kernel.radius <- radius.seq[radius.idx + 1]
  exp.dist <- exp(-dist.mat^2*kernel.radius)

  ## smoothing statistics
  smooth.stat <- exp.dist %*% region.stat.obs

  ## step 4: permutation
  set.seed(seed)
  Y.permute <- sapply(seq_len(permute.time), function(x) sample(group, replace = F))
  region.stat.permute <- region.stat(X = X,
                                     Y = Y.permute)
  smooth.stat.permute <- exp.dist %*% region.stat.permute

  ## p-value calculation
  smooth.stat.all <- cbind(smooth.stat, smooth.stat.permute)
  smooth.stat.pval <- apply(smooth.stat.all, 1, function(x) mean(abs(x[1]) < abs(x[-1])))
  smooth.stat.pval[smooth.stat.pval == 0] <- 1/permute.time
  smooth.stat.pval.BH <- p.adjust(smooth.stat.pval, method = "BH")

  # other info for visualization
  avg.expr <- round(colMeans(X), 2)
  direction <- ifelse(sign(smooth.stat) > 0, "+", "-")
  marginal.pval <- 10^-abs(region.stat.obs)
  marginal.pval.BH <- p.adjust(marginal.pval, method = "BH")
  log2.FC <- apply(X, 2,
                   function(x) log2(mean(x[group == 1])/mean(x[group == 0])))
  output.df <- data.frame(X.info,
                          avg.expr = avg.expr,
                          direction = direction,
                          smoothing.pval.BH = smooth.stat.pval.BH,
                          marginal.pval.BH = marginal.pval.BH,
                          log2.FC = log2.FC)
  return(output.df)
}
