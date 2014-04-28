\name{testTransform}
\alias{testTransform}
\alias{testTransform.powerTransform}

\title{Likelihood-Ratio Tests for Univariate or Multivariate Power Transformations
to Normality}
\description{
  \code{testTransform} computes likelihood ratio tests for particular transformations
based on \code{powerTransform} objects.
}
\usage{
testTransform(object, lambda)

\S3method{testTransform}{powerTransform}(object, lambda=rep(1, dim(object$y)[2]))
}

\arguments{
  \item{object}{An object created by a call to 
\code{estimateTransform} or \code{powerTransform}.}
  \item{lambda}{A vector of values of length equal to the number of
variables to be transformed.}
}

\details{The function \code{\link{powerTransform}} is used to estimate a power
transformation for a univariate or multivariate sample or multiple linear
regression problem, using the method of Box and Cox (1964).  It is usual to
round the estimates to nearby convenient values, and this function is use to
compulte a likelihood ratio test for values of the transformation parameter
other than the ml estimate.  This is a generic function, but with only one
method, for objects of class \code{powerTransform}.}

\value{
A data frame with one row giving the value of the test statistic, its
degrees of freedom, and a p-value.  The test is the likelihood ratio test, 
comparing the value of the log-likelihood at the hypothesized value to 
the value of the log-likelihood at the maximum likelihood estimate.
}
\references{Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations. \emph{Journal
of the Royal Statisistical Society, Series B}. 26 211-46.

Cook, R. D. and Weisberg, S. (1999)  \emph{Applied Regression Including Computing
and Graphics}.  Wiley.

Fox, J. and Weisberg, S. (2011) 
\emph{An R Companion to Applied Regression}, Second Edition, Sage.

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{\code{\link{powerTransform}}. }
\examples{
summary(a3 <- powerTransform(cbind(len, ADT, trks, sigs1) ~ hwy, Highway1))
# test lambda = (0 0 0 -1)
testTransform(a3, c(0, 0, 0, -1))
}
\keyword{ regression}% at least one, from doc/KEYWORDS
