\name{panel.car}
\alias{panel.car}

\title{Panel Function for Coplots}
\description{
  a panel function for use with \code{coplot} that plots points, a lowess line,
  and a regression line.
}
\usage{
panel.car(x, y, col, pch, cex=1, span=0.5, lwd=2, 
  reg.line=lm, lowess.line=TRUE, ...)
}

\arguments{
  \item{x}{vector giving horizontal coordinates.}
  \item{y}{vector giving vertical coordinates.}
  \item{col}{point color.}
  \item{pch}{plotting character for points.}
  \item{cex}{character expansion factor for points.}
  \item{span}{span for lowess smoother.}
  \item{lwd}{line width, default is \code{2}.}
  \item{reg.line}{function to compute coefficients of regression line,
    or \code{FALSE} for no line.}
  \item{lowess.line}{if \code{TRUE} plot lowess smooth.}
  \item{\dots}{other arguments to pass to functions \code{lines}
    and \code{regLine}.}
}

\value{
  \code{NULL}. This function is used for its side effect: producing
  a panel in a coplot.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[graphics]{coplot}}, \code{\link{regLine}}}

\examples{
coplot(prestige ~ income|education, panel=panel.car, 
  col="red", data=Prestige)
}

\keyword{aplot}
