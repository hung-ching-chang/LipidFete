#' @title Identify all neighbors of input feature
#'
#' @description Identify all neighbors of input feature
#'
#' @param self A vector of (x,y) coordinates
#' @param other A vector of (x,y) coordinates
#' @param type Either "right", "left", "top", or "bottom"
#'
#' @details write XXX
#'
#' @return An indicator of whether self and other are neighbors.
#'
#' @examples
#' # two dimensional location
#' find.neighbor(c(1,1),c(1,2), type = "top")
#' find.neighbor(c(1,1),c(1,2), type = "bottom")
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{region.plot.2D}, \link{build.wall}}
#' @export
find.neighbor <- function(self, other,
                          type = c("right", "left", "top", "bottom")){
  if(type == "right"){
    result <- (other[2] == self[2]) & (other[1] - self[1] == 1)
    return(result)   # other is to the right of self
  }
  if(type == "left"){
    result <- (other[2] == self[2]) & (other[1] - self[1] == -1)
    return(result)   # other is to the right of self
  }
  if(type == "top"){
    result <- (other[1] == self[1]) & (other[2] - self[2] == 1)
    return(result)  # other is on top of self
  }
  if(type == "bottom"){
    result <- (other[1] == self[1]) & (other[2] - self[2] == -1)
    return(result)  # other is on top of self
  }
}
