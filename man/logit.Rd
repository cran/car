\name{logit}
\alias{logit}

\title{Logit Transformation}
\description{
  Compute the logit transformation of proportions or percentages.
}
\usage{
logit(p, percents=range.p[2] > 1, adjust)
}

\arguments{
  \item{p}{numeric vector or array of proportions or percentages.}
  \item{percents}{\code{TRUE} for percentages.}
  \item{adjust}{adjustment factor to avoid proportions of 0 or 1;
    defaults to \code{0} if there are no such proportions in the data,
    and to \code{.025} if there are.}
}
\details{
  Computes the logit transformation logit \eqn{=\log [p/(1-p)]}{= log[p/(1 - p)]}
  for the proportion \eqn{p}. 
  
  If \eqn{p=0}{p = 0} or \eqn{1}, then the logit
  is undefined. \code{logit} can remap the proportions to the interval 
  \code{(adjust, 1 - adjust)} prior to the transformation. If it adjusts the data
  automatically, \code{logit} will print a warning message.
}
\value{
  a numeric vector or array of the same shape and size as \code{p}.
}

\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{probabilityAxis}}}

\examples{
options(digits=4)
logit(.1*0:10)
##  [1] -3.6636 -1.9924 -1.2950 -0.8001 -0.3847  0.0000  0.3847
##  [8]  0.8001  1.2950  1.9924  3.6636
## Warning message: 
## proportions remapped to (0.025, 0.975) in: logit(0.1 * 0:10) 

logit(.1*0:10, adjust=0)
##  [1]    -Inf -2.1972 -1.3863 -0.8473 -0.4055  0.0000  0.4055
##  [8]  0.8473  1.3863  2.1972     Inf
}

\keyword{manip}
