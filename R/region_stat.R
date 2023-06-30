#' @title Calculate region statistics
#'
#' @description Calculate region statistics of every lipidomic feature. The region statistics
#' is signed p-value, which is constituted by p-value of marginal test and sign of mean difference
#' between groups.
#'
#' @param X A \code{matrix} of lipid expression data. Row is samples, and column is
#' lipid feature such as double bond and chain length.
#' @param Y A \code{matrix} of observed or permuted grouping variable. Each column is
#' a grouping result. Only allow 0 and 1.
#'
#' @details write XXX
#'
#' @return A \code{matrix} of region statistics calculated as described above.
#' Each column is a set of region statistics corresponding with each grouping result.
#'
#' @examples
#' # one dimensional Y.
#' X0 <- matrix(rnorm(100,0,1), nrow = 10)
#' X1 <- matrix(rnorm(100,2,1), nrow = 10)
#' X <- rbind(X0,X1)
#' Y <- matrix(rep(c(0,1), each = 10), ncol = 1)
#'
#' # calculate region statistics
#' rstat <- region.stat(X = X, Y = Y)
#'
#' # high dimensional Y.
#' X0 <- matrix(rnorm(100,0,1), nrow = 10)
#' X1 <- matrix(rnorm(100,2,1), nrow = 10)
#' X <- rbind(X0,X1)
#' group <- rep(c(0,1), each = 10)
#' Y <- sapply(seq_len(20),
#'             function(i) sample(group, replace = F))
#'
#' # calculate region statistics
#' rstat <- region.stat(X = X, Y = Y)
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{LipidFete.test}}
#' @export
region.stat <- function(X, Y, ...){
  n1 <- sum(Y[,1])
  n2 <- sum(1-Y[,1])

  # Calculate sample mean/variance for each group
  G1 <- Y / rep(colSums(Y), each = nrow(Y))
  G2 <- (1-Y) / rep(colSums(1-Y), each = nrow(Y))
  mean.group1 <- t(X) %*% G1
  mean.group2 <- t(X) %*% G2
  var.group1 <- (t(X^2) %*% Y  - n1*mean.group1^2)/(n1 - 1)
  var.group2 <- (t(X^2) %*% (1-Y)  - n2*mean.group2^2)/(n2 - 1)

  # Calculate pooled variance and t-statistic with equal variance
  pooled.var <- ((n1-1)*var.group1 + (n2-1)*var.group2)/(n1+n2-2)
  t.stat <- (mean.group1-mean.group2) / sqrt((1/n1+1/n2)*pooled.var)

  # Calculate region statistic
  df <- n1 + n2 - 2
  t.test.pval <- 2*pt(abs(t.stat), df = df, lower.tail = FALSE)
  res <- -log10(t.test.pval) * sign(t.stat)
  return(res)
}
