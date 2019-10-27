\name{carPalette}
\alias{carPalette}
\title{
Set or Retrieve \pkg{car} Package Color Palette
}
\description{
This function is used to set or retrieve colors to be used in \pkg{car} package graphics functions.
}
\usage{
carPalette(palette)
}
\arguments{
  \item{palette}{
if missing, returns the colors that will be used in \pkg{car} graphics; if present, the colors to be used in graphics will be set. The \code{palette} argument may also be one of \code{"car"} or \code{"default"} to use the default \pkg{car} palette (defined below), \code{"R"} to use the default R palette, or \code{"colorblind"} to use a colorblind-friendly palette (from \url{https://jfly.uni-koeln.de/color/}). 
}
}
\details{
This function sets or returns the value of \code{options(carPalette=pallete)} that will be use in \pkg{car} graphics functions to determine colors.  The default is \code{c("black", "blue", "magenta", "cyan", "orange", "gray", "green3", "red"))}, which is nearly a permutation of the colors returned by the standard \code{palette} function that minimizes the use of red and green in the same graph, and that substitutes orange for the often hard to see yellow.
}
\value{
Invisibly returns the previous value of the \pkg{car} palette.
}

\references{
  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{Sanford Weisberg and John Fox}

\seealso{
\code{\link{palette}}, \code{\link{colors}}
}
\examples{
# Standard color palette
palette()
# car standard color palette
carPalette()
# set colors to all black
carPalette(rep("black", 8))
# Use a custom color palette with 12 distinct colors
carPalette(sample(colors(distinct=TRUE), 12))
# restore default
carPalette("default")
}


\keyword{ color }
