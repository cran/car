\name{regLine}
\alias{regLine}

\title{Plot Regression Line}
\description{
  Plots a regression line on a scatterplot; the line is plotted 
  between the minimum and maximum x-values.
}
\usage{
regLine(mod, col=carPalette()[2], lwd=2, lty=1,...)
}

\arguments{
  \item{mod}{a model, such as produced by \code{lm},
    that responds to the \code{coef} function
    by returning a 2-element vector, whose elements are
    interpreted respectively as the
    intercept and slope of a regresison line.}
  \item{col}{color for points and lines; the default is the \emph{second} entry
    in the current \pkg{car} palette (see \code{\link{carPalette}}
    and \code{\link[graphics]{par}}).}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{lty}{line type; default is \code{1}, a solid line (see \code{\link[graphics]{par}}).}
  \item{\dots}{optional arguments to be passed to the
    \code{lines} plotting function.}
}
\details{
  In contrast to \code{abline}, this function plots only over
  the range of the observed x-values. The x-values are extracted from
  \code{mod} as the second column of the model matrix. 
}
\value{
  \code{NULL}. This function is used for its side effect: adding
  a line to the plot.
}

\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{\code{\link[graphics]{abline}}, \code{\link[graphics]{lines}}}

\examples{
plot(repwt ~ weight, pch=c(1,2)[sex], data=Davis)
regLine(lm(repwt~weight, subset=sex=="M", data=Davis))
regLine(lm(repwt~weight, subset=sex=="F", data=Davis), lty=2)
}

\keyword{aplot}
