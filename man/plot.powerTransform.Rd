\name{plot.powerTransform}
\Rdversion{1.1}
\alias{plot.powerTransform}
\title{ plot Method for powerTransform Objects
}
\description{
This function provides a simple function for plotting data using power 
transformations.

}
\usage{
\method{plot}{powerTransform}(x, z = NULL, round = TRUE, plot = pairs, ...)
}

\arguments{
  \item{x}{name of the power transformation object
}
  \item{z}{
  Additional variables of the same length as those used to get the transformation
  to be plotted, default is \code{NULL}.
}
  \item{round}{
  If \code{TRUE}, the default, use rounded transforms, if \code{FALSE} use the MLEs.
}
  \item{plot}{
  Plotting method.  Default is \code{pairs}.  Another possible choice is
  \code{scatterplot.matrix} from the \code{car} package.
}
  \item{\dots}{
  Optional arguments passed to the plotting method
}
}
\details{
The data used to estimate transformations using \code{powerTransform} are plotted
in the transformed scale.
}
\value{
None.  Produces a graph as a side-effect.
}
\references{
Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.

Fox, J. and Weisberg, S. (2011) \emph{An R Companion to Applied Linear Regression},
Second Edition, Sage.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }


\seealso{
\code{\link{powerTransform}}
}
\examples{
summary(a3 <- powerTransform(cbind(len, ADT, trks, sigs1) ~ hwy, Highway1))
with(Highway1, plot(a3, z=rate, col=as.numeric(hwy)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ regression }% __ONLY ONE__ keyword per line
