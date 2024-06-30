#' @title Visualize LipidFete test result of 2D feature
#'
#' @description Visualize LipidFete test result of 2D feature for identifying tendencies of lipid features.
#'
#' @param test.result Result from LipidFete.test function
#' @param pval.thres A threshold for smoothing.pval
#' @param x.distance A value of unit distance in x-axis
#' @param y.distance A value of unit distance in y-axis
#' @param log2.FC.thres A threshold for log2-transformed fold change
#'
#' @details write XXX
#'
#' @return A heatmap presenting significant lipidomic features, colored by log2 fold change, and annotated by marginal p-value.
#'
#' @examples
#' data(lipid2D)
#' X <- t(as.matrix(lipid2D[,-c(1:2)]))
#' X.info <- lipid2D[,1:2]
#' group <- rep(c(0, 1), c(52,32))
#' test.result <- LipidFete.test(X = X,
#'                               X.info = X.info,
#'                               group = group,
#'                               radius = 3,
#'                               own.contri = 0.5,
#'                               x.distance = 2,
#'                               y.distance = 1,
#'                               dimension = 2,
#'                               permute.time = 100000)
#'
#' region.plot.2D(test.result = test.result,
#'                pval.thres = 0.05,
#'                x.distance = 2,
#'                y.distance = 1)
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{LipidFete.test}, \link{build.wall}, \link{pval.annotation}}
#' @export
#' @import ggplot2
region.plot.2D <- function(test.result, pval.thres = 0.05, x.distance = 1, y.distance = 1,
                           log2.FC.thres = 3, ...){
  # input
  X.info <- test.result[,1:2]
  avg.expr <- test.result$avg.expr
  direction <- test.result$direction
  smoothing.pval <- test.result$smoothing.pval.BH
  marginal.pval <- test.result$marginal.pval.BH
  log2.FC <- test.result$log2.FC

  # selected region by cut point
  direction.int <- ifelse(direction == "+", 1, -1)
  selected.region.high <- selected.region.low <- rep("None", length(smoothing.pval))
  selected.region.high[smoothing.pval < pval.thres & direction.int > 0] <- "High"
  selected.region.low[smoothing.pval < pval.thres & direction.int < 0] <- "Low"

  ## plot border
  # divided total length by 2 if  only consider even length
  dist.input <- as.matrix(X.info)
  dist.input[,1] <- dist.input[,1]/x.distance
  dist.input[,2] <- dist.input[,2]/y.distance

  ## build wall (high group)
  wall.high <- lapply(seq_len(nrow(X.info)),
                      function(x) build.wall(feature.idx = x, X.info, selected.region.high,
                                             x.distance, y.distance))
  wall.high <- as.data.frame(do.call(rbind, wall.high))
  top.wall.high <- wall.high[wall.high$pos == "top", seq_len(2)]
  bottom.wall.high <- wall.high[wall.high$pos == "bottom", seq_len(2)]
  right.wall.high <- wall.high[wall.high$pos == "right", seq_len(2)]
  left.wall.high <- wall.high[wall.high$pos == "left", seq_len(2)]
  colnames(top.wall.high) <- colnames(bottom.wall.high) <- colnames(right.wall.high) <- colnames(left.wall.high) <- c("x","y")

  ## build wall (low group)
  wall.low <- lapply(seq_len(nrow(X.info)),
                     function(x) build.wall(feature.idx = x, X.info, selected.region.low,
                                            x.distance, y.distance))
  wall.low <- as.data.frame(do.call(rbind, wall.low))
  top.wall.low <- wall.low[wall.low$pos == "top", seq_len(2)]
  bottom.wall.low <- wall.low[wall.low$pos == "bottom", seq_len(2)]
  right.wall.low <- wall.low[wall.low$pos == "right", seq_len(2)]
  left.wall.low <- wall.low[wall.low$pos == "left", seq_len(2)]
  colnames(top.wall.low) <- colnames(bottom.wall.low) <- colnames(right.wall.low) <- colnames(left.wall.low) <- c("x","y")

  # p-value annotation
  pval.annotate <- sapply(marginal.pval,
                          function(pval) pval.annotation(pval))

  # heatmap
  heatmap.df <- as.data.frame(X.info)   # x- and y-axis
  var.name <- colnames(heatmap.df)
  colnames(heatmap.df) <- c("v1","v2")
  heatmap.df$avg.expr <- avg.expr
  heatmap.df$log2.FC <- log2.FC    # color
  heatmap.df$log2.FC[log2.FC > log2.FC.thres] <- log2.FC.thres
  heatmap.df$log2.FC[log2.FC < -log2.FC.thres] <- -log2.FC.thres
  heatmap.df$pval.annotate <- pval.annotate      # label
  heatmap.df$signed.logp.smooth <- direction.int * smoothing.pval
  x.label <- seq(min(X.info[,1]), max(X.info[,1]), x.distance)
  y.label <- seq(min(X.info[,2]), max(X.info[,2]), y.distance)

  result <- ggplot(heatmap.df, aes(x = v1, y = v2))  +
    geom_point(aes(size = avg.expr, color = log2.FC)) +
    scale_colour_gradient2(high = "#F8766D", mid = "white", low = "blue", midpoint = 0) +
    geom_text(aes(x = v1, y = v2, label = pval.annotate),
              color = "black", size = 4) +
    scale_x_continuous(breaks=x.label) +
    scale_y_continuous(breaks=y.label) +
    labs(x = var.name[1], y = var.name[2]) +
    scale_size(range = c(5, 10))

  # add high group wall
  result <-  result +
    geom_segment(data=top.wall.high, aes(x=x-x.distance/2, xend=x+x.distance/2,
                                         y=y+y.distance/2, yend=y+y.distance/2),
                 linewidth = 1, colour = "red")+
    geom_segment(data=right.wall.high, aes(x=x+x.distance/2, xend=x+x.distance/2,
                                      y=y-y.distance/2, yend=y+y.distance/2),
                 linewidth = 1, colour = "red") +
    geom_segment(data=bottom.wall.high, aes(x=x-x.distance/2, xend=x+x.distance/2,
                                       y=y-y.distance/2, yend=y-y.distance/2),
                 linewidth = 1, colour = "red")+
    geom_segment(data=left.wall.high, aes(x=x-x.distance/2, xend=x-x.distance/2,
                                     y=y-y.distance/2, yend=y+y.distance/2),
                 linewidth = 1, colour = "red")

  # add low group wall
  result <-  result +
    geom_segment(data=top.wall.low, aes(x=x-x.distance/2, xend=x+x.distance/2,
                                         y=y+y.distance/2, yend=y+y.distance/2),
                 linewidth = 1, colour = "blue")+
    geom_segment(data=right.wall.low, aes(x=x+x.distance/2, xend=x+x.distance/2,
                                           y=y-y.distance/2, yend=y+y.distance/2),
                 linewidth = 1, colour = "blue") +
    geom_segment(data=bottom.wall.low, aes(x=x-x.distance/2, xend=x+x.distance/2,
                                            y=y-y.distance/2, yend=y-y.distance/2),
                 linewidth = 1, colour = "blue")+
    geom_segment(data=left.wall.low, aes(x=x-x.distance/2, xend=x-x.distance/2,
                                          y=y-y.distance/2, yend=y+y.distance/2),
                 linewidth = 1, colour = "blue")
  return(result)
}
