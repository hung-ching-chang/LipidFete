#' @title Calculate p-value of lipidomic features by Bayesian hierarchical model
#'
#' @description Calculate posterior probability of lipidomic features by Bayesian hierarchical model.
#'
#' @param X A \code{matrix} of lipid expression data. Row is samples, and column is
#' lipid feature such as double bond and chain length.
#' @param X.info A \code{data.frame} or \code{matrix} of lipid feature names. (Only allow 1 or 2 columns)
#' @param group A vector of observed grouping variable. (Only allow 0 and 1)
#' @param x.distance A value of unit distance in x-axis
#' @param y.distance A value of unit distance in y-axis
#' @param dimension A value to indicate the dimension of lipid features (1 or 2)
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
#' LipidFete.Bayes(X = X,
#'                 X.info = X.info,
#'                 group = group,
#'                 radius = 1.5,
#'                 x.distance = 2,
#'                 y.distance = 1,
#'                 dimension = 2)
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{region.plot.1D}, \link{region.plot.2D}}
#' @export
#' @import MASS
LipidFete.Bayes <- function(X, X.info, group, radius = 1.5,
                            x.distance = 1, y.distance = 1, dimension = c(1,2),
                            MC.sample = 100000, seed = 1234, ...){
  ## step 1: likelihood, prior parameters
  n.grp0 <- sum(group == 0); n.grp1 <- sum(group == 1)
  set.seed(seed)
  
  # equal variance assumption
  #sigma.grp0 <- sigma.grp1 <- apply(X, 2, var)/(n.grp0+n.grp1)
  sigma.grp0 <- apply(X[group == 0,], 2, var)/n.grp0
  sigma.grp1 <- apply(X[group == 1,], 2, var)/n.grp1
  sigma.max <- pmax(sigma.grp0, sigma.grp1)
  sigma.grp0 <- sigma.grp1 <- sigma.max
  
  # estimate sample variance from all samples when sigma.grp0 or sigma.grp0 == 0
  sigma.two.grp <- apply(X, 2, var)/(n.grp0+n.grp1)
  sigma.grp0[sigma.grp0 == 0] <- sigma.two.grp[sigma.grp0 == 0]
  sigma.grp1[sigma.grp1 == 0] <- sigma.two.grp[sigma.grp1 == 0]
  
  # prior of theta
  theta.grp0 <- colMeans(X[group == 0,])
  theta.grp1 <- colMeans(X[group == 1,])
  library(plyr)
  dist.df <- data.frame(X.info[,1]/x.distance, X.info[,2]/y.distance)
  dist.mat <- as.matrix(dist(dist.df))
  diag(dist.mat) <- Inf
  g.info <- alply(dist.mat, 1,
                  function(x) which(x <= radius)) # neighbor/graph info
  # estimated mean value
  theta.bar.grp0 <- sapply(g.info,
                           function(x) mean(theta.grp0[x]))
  theta.bar.grp1 <- sapply(g.info,
                           function(x) mean(theta.grp1[x]))
  var.grp0 <- sapply(g.info,
                     function(x) var(theta.grp0[x])/length(x))
  var.grp1 <- sapply(g.info,
                     function(x) var(theta.grp1[x])/length(x))
  var.max <- pmax(var.grp0, var.grp1)
  var.grp0 <- var.grp1 <- var.max
  
  theta.bar.grp0[is.na(theta.bar.grp0)] <- theta.bar.grp1[is.na(theta.bar.grp1)] <- 0
  var.grp0[is.na(var.grp0)] <- var.grp1[is.na(var.grp1)] <- Inf
  
  ## step 2: posterior distribution
  mu.grp0 <- (theta.grp0/sigma.grp0 + theta.bar.grp0/var.grp0)/(1/sigma.grp0 + 1/var.grp0)
  mu.grp1 <- (theta.grp1/sigma.grp1 + theta.bar.grp1/var.grp1)/(1/sigma.grp1 + 1/var.grp1)
  tau.grp0 <- 1/(1/sigma.grp0 + 1/var.grp0)
  tau.grp1 <- 1/(1/sigma.grp1 + 1/var.grp1)
  
  # check
  # message("mean0: ", round(theta.grp0[58],7),
  #         "; mean1: ", round(theta.grp1[58],7),"\n",
  #         "prior0: ", round(theta.bar.grp0[58],7),
  #         "; prior1: ", round(theta.bar.grp1[58],7),"\n",
  #         "post0: ", round(mu.grp0[58],7),
  #         "; post1: ", round(mu.grp1[58],7))
  # message("var0: ", round(sqrt(sigma.grp0[58]),7),
  #         "; var1: ", round(sqrt(sigma.grp1[58]),7),"\n",
  #         "prior0: ", round(sqrt(var.grp0[58]),7),
  #         "; prior1: ", round(sqrt(var.grp1[58]),7),"\n",
  #         "post0: ", round(sqrt(tau.grp0[58]),7),
  #         "; post1: ", round(sqrt(tau.grp1[58]),7))
  
  ## step 3: posterior probability
  post.grp0 <- mvrnorm(n = MC.sample,
                       mu = mu.grp0,
                       Sigma = diag(tau.grp0))
  post.grp1 <- mvrnorm(n = MC.sample,
                       mu = mu.grp1,
                       Sigma = diag(tau.grp1))
  
  post.mean.diff <- post.grp1 - post.grp0
  post.prob <- colMeans(post.mean.diff > 0)
  
  # other info for visualization
  avg.expr <- round(colMeans(X), 2)
  log2.FC <- apply(X, 2,
                   function(x) log2(mean(x[group == 1])/mean(x[group == 0])))
  output.df <- data.frame(X.info,
                          avg.expr = avg.expr,
                          post.prob = post.prob,
                          log2.FC = log2.FC)
  return(output.df)
}
