\name{logit}
\alias{logit}

\title{Logit Transformation}
\description{
  Compute the logit transformation of proportions or percentages.
}
\usage{
logit(p, percents, adjust)
}

\arguments{
  \item{p}{numeric vector or array of proportions or percentages.}
  \item{percents}{\code{TRUE} for percentages, \code{FALSE} for proportions. If the argument is missing and the largest value of \code{p} > 1, \code{percents} is set to \code{TRUE}, otherwise to \code{FALSE}.}
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
save.opt <- options(digits=4)
logit(.1*0:10)
logit(.1*0:10, adjust=0)
options(save.opt)
}

\keyword{manip}
