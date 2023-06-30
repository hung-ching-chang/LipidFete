#' @title Visualize LipidFete test result of 1D feature
#'
#' @description Visualize LipidFete test result of 1D feature for identifying tendencies of lipid features.
#'
#' @param X A \code{matrix} of lipid expression data. Row is samples, and column is
#' lipid feature such as double bond and chain length.
#' @param X.info A \code{data.frame} or \code{matrix} of lipid feature names. (Only allow 1 or 2 columns)
#' @param group A vector of observed grouping variable. (Only allow 0 and 1)
#' @param direction A vector of direction of feature expression level
#' @param smoothing.pval A vector of kernel smoothing p-value
#' @param marginal.pval A vector of marginal p-value
#' @param feature.name A vector of log2-transformed fold change
#' @param cut.point A cut point for smoothing.pval
#'
#' @details write XXX
#'
#' @return A figure presenting significant lipidomic features colored by direction of feature expression level
#'
#' @examples
#' data(lipid1D)
#' X <- t(as.matrix(lipid1D[,2:85]))
#' X.info <- lipid1D[,1]
#' group <- rep(c(0, 1), c(52,32))
#' test.result <- LipidFete.test(X = X,
#'                               X.info = X.info,
#'                               group = group,
#'                               radius = 2,
#'                               own.contri = 0.5,
#'                               x.distance = 2,
#'                               y.distance = 1,
#'                               dimension = 1,
#'                               permute.time = 100000)
#'
#' region.plot.1D(X = X,
#'                X.info = X.info,
#'                group = group,
#'                direction = test.result$direction,
#'                smoothing.pval = test.result$smoothing.pval.BH,
#'                marginal.pval = test.result$marginal.pval.BH,
#'                feature.name = colnames(X)[1],
#'                cut.point = 0.05)
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{region.stat}, \link{LipidFete.test}, \link{region.plot.2D}}
#' @export
#' @import ggplot2 cowplot
region.plot.1D <- function(X, X.info, group, direction,
                           smoothing.pval, marginal.pval,
                           feature.name, cut.point = 0.05, ...){
  # selected region by cut point
  direction.int <- ifelse(direction == "+", 1, -1)
  category <- rep("None", length(smoothing.pval))
  signed.log.smoothing.pval <- direction.int * -log10(abs(smoothing.pval))
  signed.log.marginal.pval <- direction.int * -log10(marginal.pval)
  category[signed.log.smoothing.pval > -log10(cut.point)] <- "High"
  category[signed.log.smoothing.pval < log10(cut.point)] <- "Low"
  category <- factor(category,
                     levels = c("Low", "None", "High"))

  # boxplot
  boxplot.df <- as.matrix(X)
  colnames(boxplot.df) <- X.info
  boxplot.df <- reshape2::melt(boxplot.df)
  boxplot.df$Var2 <- factor(boxplot.df$Var2,
                            levels = X.info)
  boxplot.df$group <- rep(as.character(group), ncol(X))

  barchart.df <- data.frame(xlabel = X.info,
                            log.marginal.pval = signed.log.marginal.pval,
                            log.smoothing.pval = signed.log.smoothing.pval,
                            category = category)
  barchart.df$xlabel <-  factor(barchart.df$xlabel,
                                levels = X.info)
  bar.ymin <- min(signed.log.smoothing.pval, signed.log.marginal.pval)
  bar.ymax <- max(signed.log.smoothing.pval, signed.log.marginal.pval)

  boxplot <- ggplot(boxplot.df)  +
    geom_boxplot(aes(x = Var2, y = value, fill=group)) +
    scale_fill_manual(values=c("cadetblue2", "#F8766D")) +
    xlab(feature.name)

  barchart1 <- ggplot(barchart.df, aes(x=xlabel, y = log.marginal.pval)) +
    geom_bar(stat = "identity") +
    xlab(feature.name) +
    ylim(bar.ymin, bar.ymax)

  barchart2 <- ggplot(barchart.df, aes(x=xlabel, y = log.smoothing.pval, fill = category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Low" = "cadetblue2",
                                 "None" = "grey",
                                 "High" = "#F8766D")) +
    xlab(feature.name) +
    ylim(bar.ymin, bar.ymax)

  # combine
  result <- cowplot::plot_grid(boxplot, barchart1, barchart2,
                               ncol = 1, rel_heights = c(3, 2, 2),
                               align = 'v', axis = 'lr')
  return(result)
}
