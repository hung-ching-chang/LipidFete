% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_wall.R
\name{build.wall}
\alias{build.wall}
\title{Identify all locations of walls}
\usage{
build.wall(feature.idx, X.info, selected.region, x.distance, y.distance, ...)
}
\arguments{
\item{feature.idx}{A value of feature index}

\item{X.info}{A \code{data.frame} or \code{matrix} of lipid feature names.
(Only allow 1 or 2 columns)}

\item{selected.region}{A vector of description of feature expression level.
Either "Low", "None", or "High".}

\item{x.distance}{A value of unit distance in x-axis}

\item{y.distance}{A value of unit distance in y-axis}
}
\value{
An \code{data.frame} of wall position:
\itemize{
    \item{column 1: }{lipid feature 1.}
    \item{column 2: }{lipid feature 2.}
    \item{pos: }{desciption of wall position. Either "right", "left", "top", or "bottom"}
}
}
\description{
Identify all locations of walls
}
\details{
write XXX
}
\examples{
X.info <- data.frame(x = c(1,1,1),
                     y = c(1,2,3))
feature.idx <- 1
selected.region <- c("low", "None", "None")
x.distance <- y.distance <- 1

build.wall(feature.idx, X.info, selected.region, x.distance, y.distance)
# result indicates that we shall build 4 walls around (1,1) coordinates

}
\seealso{
\code{\link{region.plot.2D}, \link{find.neighbor}}
}
\author{
Hung-Ching Chang
}
