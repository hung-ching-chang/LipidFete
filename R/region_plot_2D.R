#' @title Visualize LipidFete test result of 2D feature
#'
#' @description Visualize LipidFete test result of 2D feature for identifying tendencies of lipid features.
#'
#' @param X.info A \code{data.frame} or \code{matrix} of lipid feature names. (Only allow 1 or 2 columns)
#' @param direction A vector of direction of feature expression level
#' @param smoothing.pval A vector of kernel smoothing p-value
#' @param marginal.pval A vector of marginal p-value
#' @param log2.FC A vector of log2-transformed fold change
#' @param cut.point A cut point for smoothing.pval
#' @param x.distance A value of unit distance in x-axis
#' @param y.distance A value of unit distance in y-axis
#' @param log2.FC.cut.point A cut point for log2-transformed fold change
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
#' region.plot.2D(X.info = X.info,
#'                direction = test.result$direction,
#'                smoothing.pval = test.result$smoothing.pval.BH,
#'                marginal.pval = test.result$marginal.pval.BH,
#'                log2.FC = test.result$log2.FC,
#'                cut.point = 0.05,
#'                x.distance = 2,
#'                y.distance = 1)
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{LipidFete.test}, \link{build.wall}, \link{pval.annotation}}
#' @export
#' @import ggplot2
region.plot.2D <- function(X.info, direction, smoothing.pval, marginal.pval,
                           log2.FC, cut.point = 0.05, x.distance = 1, y.distance = 1,
                           log2.FC.cut.point = 3, ...){
  # selected region by cut point
  direction.int <- ifelse(direction == "+", 1, -1)
  selected.region <- rep("None", length(smoothing.pval))
  selected.region[smoothing.pval < cut.point & direction.int > 0] <- "High"
  selected.region[smoothing.pval < cut.point & direction.int < 0] <- "Low"

  ## plot border
  # divided total length by 2 if  only consider even length
  dist.input <- as.matrix(X.info)
  dist.input[,1] <- dist.input[,1]/x.distance
  dist.input[,2] <- dist.input[,2]/y.distance

  ## build wall
  wall <- lapply(seq_len(nrow(X.info)),
                 function(x) build.wall(feature.idx = x, X.info, selected.region, x.distance, y.distance))
  wall <- as.data.frame(do.call(rbind, wall))
  top.wall <- wall[wall$pos == "top", seq_len(2)]
  bottom.wall <- wall[wall$pos == "bottom", seq_len(2)]
  right.wall <- wall[wall$pos == "right", seq_len(2)]
  left.wall <- wall[wall$pos == "left", seq_len(2)]
  colnames(top.wall) <- colnames(bottom.wall) <- colnames(right.wall) <- colnames(left.wall) <- c("x","y")

  # p-value annotation
  pval.annotate <- sapply(marginal.pval,
                          function(pval) pval.annotation(pval))

  # heatmap
  heatmap.df <- as.data.frame(X.info)   # x- and y-axis
  var.name <- colnames(heatmap.df)
  colnames(heatmap.df) <- c("v1","v2")
  heatmap.df$log2.FC <- log2.FC    # color
  heatmap.df$log2.FC[log2.FC > log2.FC.cut.point] <- log2.FC.cut.point
  heatmap.df$log2.FC[log2.FC < -log2.FC.cut.point] <- -log2.FC.cut.point
  heatmap.df$pval.annotate <- pval.annotate      # label
  heatmap.df$signed.logp.smooth <- direction.int * smoothing.pval
  x.label <- seq(min(X.info[,1]), max(X.info[,1]), x.distance)
  y.label <- seq(min(X.info[,2]), max(X.info[,2]), y.distance)

  result <- ggplot(heatmap.df, aes(x = v1, y = v2))  +
    geom_tile(aes(fill=log2.FC)) +
    scale_fill_gradient2(high="#F8766D",mid="white",low="blue", midpoint = 0) +
    geom_text(aes(x = v1, y = v2, label = pval.annotate),
              color = "black", size = 4) +
    scale_x_continuous(breaks=x.label) +
    scale_y_continuous(breaks=y.label) +
    labs(x = var.name[1], y = var.name[2]) +
    geom_segment(data=top.wall, aes(x=x-x.distance/2, xend=x+x.distance/2,
                                    y=y+y.distance/2, yend=y+y.distance/2), linewidth = 1)+
    geom_segment(data=right.wall, aes(x=x+x.distance/2, xend=x+x.distance/2,
                                      y=y-y.distance/2, yend=y+y.distance/2), linewidth = 1) +
    geom_segment(data=bottom.wall, aes(x=x-x.distance/2, xend=x+x.distance/2,
                                       y=y-y.distance/2, yend=y-y.distance/2), linewidth = 1)+
    geom_segment(data=left.wall, aes(x=x-x.distance/2, xend=x-x.distance/2,
                                     y=y-y.distance/2, yend=y+y.distance/2), linewidth = 1)
  return(result)
}
