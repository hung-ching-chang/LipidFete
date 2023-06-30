#' @title Identify all locations of walls
#'
#' @description Identify all locations of walls
#'
#' @param feature.idx A value of feature index
#' @param X.info A \code{data.frame} or \code{matrix} of lipid feature names.
#' (Only allow 1 or 2 columns)
#' @param selected.region A vector of description of feature expression level.
#' Either "Low", "None", or "High".
#' @param x.distance A value of unit distance in x-axis
#' @param y.distance A value of unit distance in y-axis
#'
#' @details write XXX
#'
#' @return An \code{data.frame} of wall position:
#' \itemize{
#'     \item{column 1: }{lipid feature 1.}
#'     \item{column 2: }{lipid feature 2.}
#'     \item{pos: }{desciption of wall position. Either "right", "left", "top", or "bottom"}
#' }
#'
#' @examples
#' X.info <- data.frame(x = c(1,1,1),
#'                      y = c(1,2,3))
#' feature.idx <- 1
#' selected.region <- c("low", "None", "None")
#' x.distance <- y.distance <- 1
#'
#' build.wall(feature.idx, X.info, selected.region, x.distance, y.distance)
#' # result indicates that we shall build 4 walls around (1,1) coordinates
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{region.plot.2D}, \link{find.neighbor}}
#' @export
build.wall <- function(feature.idx, X.info, selected.region,
                       x.distance, y.distance, ...){
  ## find neighbors
  dist.input <- as.matrix(X.info)
  dist.input[,1] <- dist.input[,1]/x.distance
  dist.input[,2] <- dist.input[,2]/y.distance
  right.neighbor <- apply(dist.input, 1,
                          function(x) find.neighbor(dist.input[feature.idx,],x, type = "right"))
  left.neighbor <- apply(dist.input, 1,
                         function(x) find.neighbor(dist.input[feature.idx,],x, type = "left"))
  top.neighbor <- apply(dist.input, 1,
                        function(x) find.neighbor(dist.input[feature.idx,],x, type = "top"))
  bottom.neighbor <- apply(dist.input, 1,
                           function(x) find.neighbor(dist.input[feature.idx,],x, type = "bottom"))

  ## build wall
  wall.out <- data.frame(matrix(ncol = 3, nrow = 0))
  X.info.df <- as.data.frame(X.info)
  ## right wall
  if(sum(right.neighbor) == 1){ # if has right neighbor
    if(selected.region[feature.idx] != selected.region[right.neighbor]){  # if they are assigned to diff group
      wall.out <- rbind(wall.out, c(X.info.df[feature.idx,], pos = "right"))
    }
  }else{
    if(selected.region[feature.idx] != "None"){  # if they are assigned to diff group
      wall.out <- rbind(wall.out, c(X.info.df[feature.idx,], pos = "right"))
    }
  }
  ## top wall
  if(sum(top.neighbor) == 1){  # if has top neighbor
    if(selected.region[feature.idx] != selected.region[top.neighbor]){
      wall.out <- rbind(wall.out, c(X.info.df[feature.idx,], pos = "top"))
    }
  }else{
    if(selected.region[feature.idx] != "None"){  # if they are assigned to diff group
      wall.out <- rbind(wall.out, c(X.info.df[feature.idx,], pos = "top"))
    }
  }
  ## left wall
  ## (only check selected region & no-left neighbor)
  if(selected.region[feature.idx] != "None" & sum(left.neighbor) == 0){
    wall.out <- rbind(wall.out, c(X.info.df[feature.idx,], pos = "left"))
  }
  ## bottom wall
  ## (only check selected region & no-bottom neighbor)
  if(selected.region[feature.idx] != "None" & sum(bottom.neighbor) == 0){
    wall.out <- rbind(wall.out, c(X.info.df[feature.idx,], pos = "bottom"))
  }
  colnames(wall.out) <- c(colnames(X.info.df), "pos")
  return(wall.out)
}
