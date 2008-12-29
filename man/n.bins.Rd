\name{n.bins}
\alias{n.bins}

\title{Number of Bins for Histogram}
\description{
  Several rules for calculating the number of bins to use for
  a histogram.
}
\usage{
n.bins(x, rule=c("freedman.diaconis", "sturges", "scott", "simple"))
}

\arguments{
  \item{x}{numeric vector, variable for histogram}
  \item{rule}{see below.}
}
\details{
  \describe{
    \item{\code{"freedman.diaconis"}:}{\eqn{n^{1/3}range/2\times IQR}{(n^1/3 * range)/(2 * IQR)}.}
    \item{\code{"sturges"}:}{\eqn{ceiling(\log _{2}n+1)}{ceiling(log2n + 1)}.}
    \item{\code{"scott"}:}{\eqn{ceiling(n^{1/3}range/3.5s)}{ceiling[n^1/3 * range/(3.5 * s)]}.}
    \item{\code{"simple"}:}{\eqn{floor(10\log _{10}n)}{floor[10 * log10(n)]} for \eqn{n>100}{n > 100},
     or \eqn{floor(2*sqrt{n})}{floor[2*sqrt(n)]} for \eqn{n\leq 100}{n <= 100}.}
    }
  where \eqn{n} is the number of observations, \eqn{range} is the range of \code{x},
  \eqn{IQR} is the inter-quartile range of \code{x}, and \eqn{s} is the standard
  deviation of \code{x}.  
}

\value{
  the number of bins.
}

\references{
  Freedman, D. and Diaconis, P. (1981)
  On the histogram as a density estimator.
  \emph{Zeitschrift fur Wahrscheinlichkeitstheorie und verwandte Gebiete}
  \bold{57}, 453--476.
  
  Scott, D. W. (1979)
  On optimal and data based-histograms.
  \emph{Biometrika} \bold{66}, 605--610.
  
  Venables, W. N. and Ripley, B. D. (1999)
  \emph{Modern Applied Statistics with S-PLUS, Third Edition}, Springer.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[graphics]{hist}}}

\examples{
attach(UN)
n.bins(gdp)
## [1] 19
n.bins(gdp, "sturges")
## [1] 9
n.bins(gdp, "scott")
## [1] 8
n.bins(gdp, "simple")
## [1] 22
hist(gdp)
hist(gdp, nclass=n.bins(gdp))
}

\keyword{dplot}
