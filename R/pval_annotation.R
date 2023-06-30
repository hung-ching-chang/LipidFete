#' @title Annotate p-value with four thresholds
#'
#' @description Annotate p-value with four thresholds.
#'
#' @param pval A vector of p-value
#'
#' @details write XXX
#'
#' @return An \code{character} of p-value annotation:
#' \itemize{
#'     \item{*** }{represents p-value < 0.001.}
#'     \item{** }{represents 0.001 < p-value < 0.01.}
#'     \item{* }{represents 0.01 < p-value < 0.05.}
#'     \item{∆ }{represents 0.05 < p-value < 0.1.}
#' }
#'
#' @examples
#' pval <- 0.04
#' pval.annotation(pval)
#'
#' @author Hung-Ching Chang
#' @seealso \code{\link{region.plot.2D}}
#' @export
pval.annotation <- function(pval){
  if(is.na(pval)){
    return(NA)
  }else if(pval < 0.001){
    return("***")
  }else if(pval < 0.01){
    return("**")
  }else if(pval < 0.05){
    return("*")
  }else if(pval < 0.1){
    return("∆")
  }
  return(NA)
}
