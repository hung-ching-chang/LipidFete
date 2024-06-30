#' @title Visualize LipidFete Baysian hierarchical model result of 2D feature
#'
#' @description Visualize LipidFete Baysian hierarchical model result of 2D feature for identifying tendencies of lipid features.
#'
#' @param test.result Result from LipidFete.Bayes function
#' @param prob.thres A threshold for posterior probability
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
#' test.result <-  LipidFete.Bayes(X = X,
#'                 X.info = X.info,
#'                 group = group,
#'                 x.distance = 2,
#'                 y.distance = 1,
#'                 dimension = 2)
#'
#' region.plot.2D(test.result = test.result,
#'                prob.thres = 0.8,
#'                x.distance = 2,
#'                y.distance = 1)
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{LipidFete.Bayes}, \link{build.wall}, \link{pval.annotation}}
#' @export
#' @import ggplot2
region.plot.2D <- function(test.result, prob.thres = 0.8, x.distance = 1, y.distance = 1,
                           log2.FC.thres = 3, ...){
  # input
  X.info <- test.result[,1:2]
  avg.expr <- test.result$avg.expr
  post.prob <- test.result$post.prob
  log2.FC <- test.result$log2.FC
  
  #direction.int <- ifelse(post.prob >=0.5, 1, -1)
  
  # selected region by cut point
  selected.region <- rep("None", length(post.prob))
  selected.region[post.prob > prob.thres] <- "High"
  selected.region[post.prob < (1-prob.thres)] <- "Low"
  
  ## plot border
  # divided total length by 2 if  only consider even length
  dist.input <- as.matrix(X.info)
  dist.input[,1] <- dist.input[,1]/x.distance
  dist.input[,2] <- dist.input[,2]/y.distance
  
  ## build wall
  wall <- lapply(seq_len(nrow(X.info)),
                 function(x) build.wall(feature.idx = x, X.info, 
                                        selected.region, 
                                        x.distance, y.distance))
  wall <- as.data.frame(do.call(rbind, wall))
  top.wall <- wall[wall$pos == "top", seq_len(2)]
  bottom.wall <- wall[wall$pos == "bottom", seq_len(2)]
  right.wall <- wall[wall$pos == "right", seq_len(2)]
  left.wall <- wall[wall$pos == "left", seq_len(2)]
  colnames(top.wall) <- colnames(bottom.wall) <- colnames(right.wall) <- colnames(left.wall) <- c("x","y")
  
  # heatmap
  heatmap.df <- as.data.frame(X.info)   # x- and y-axis
  var.name <- colnames(heatmap.df)
  colnames(heatmap.df) <- c("v1","v2")
  heatmap.df$avg.expr <- avg.expr
  heatmap.df$log2.FC <- log2.FC    # color
  heatmap.df$log2.FC[log2.FC > log2.FC.thres] <- log2.FC.thres
  heatmap.df$log2.FC[log2.FC < -log2.FC.thres] <- -log2.FC.thres
  x.label <- seq(min(X.info[,1]), max(X.info[,1]), x.distance)
  y.label <- seq(min(X.info[,2]), max(X.info[,2]), y.distance)
  
  result <- ggplot(heatmap.df, aes(x = v1, y = v2))  +
    geom_point(aes(size = avg.expr, color = log2.FC)) +
    scale_colour_gradient2(high = "#F8766D", mid = "grey70", low = "blue", midpoint = 0) +
    scale_x_continuous(breaks=x.label) +
    scale_y_continuous(breaks=y.label) +
    labs(x = var.name[1], y = var.name[2]) +
    scale_size(range = c(2, 10)) +
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
