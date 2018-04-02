\name{testTransform}
\alias{testTransform}
\alias{testTransform.powerTransform}
\alias{testTransform.lmerModpowerTransform}
\alias{testTransform.bcnPowerTransformlmer}

\title{Likelihood-Ratio Tests for Univariate or Multivariate Power Transformations
to Normality}
\description{
  \code{testTransform} computes likelihood ratio tests for particular values of the power parameter
based on \code{powerTransform} objects.
}
\usage{
testTransform(object, lambda)

\S3method{testTransform}{powerTransform}(object, lambda=rep(1, dim(object$y)[2]))

\S3method{testTransform}{lmerModpowerTransform}(object, lambda=1)

\S3method{testTransform}{bcnPowerTransformlmer}(object, lambda=1)

}

\arguments{
  \item{object}{An object created by a call to \code{powerTransform}.}
  \item{lambda}{A vector of powers of length equal to the number of
variables transformed.}
}

\details{The function \code{\link{powerTransform}} is used to estimate a power
transformation for a univariate or multivariate sample or multiple linear
regression problem, using the method of Box and Cox (1964).  It is usual to
round the estimates to nearby convenient values, and this function is use to
compulte a likelihood ratio test for values of the transformation parameter
other than the ml-type estimate.

For one-parameter families of transformations, namely the Box-Cox power family \code{\link{bcPower}} and the Yeo-Johnson power family \code{\link{yjPower}}, this function computes a test based on twice the difference in the log-likelihood between the maximum likelihood-like estimate and the log-likelihood evaluated at the value of \code{lambda} specified.

For the \code{\link{bcnPower}} Box-Cox power with negatives allowed, the test is based on the profile loglikelihood maximizing over the location (or \code{gamma}) parameter(s).  Thus, \code{gamma} is treated as a nusiance parameter.}

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

Fox, J. and Weisberg, S. (2019)
\emph{An R Companion to Applied Regression}, Third Edition, Sage.

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{\code{\link{powerTransform}} and \code{\link{bcnPower}} for examples of the use of this function and other tests that might be of interest in some circumstances. }
\examples{
summary(a3 <- powerTransform(cbind(len, adt, trks, sigs1) ~ htype, Highway1))
# test lambda = (0 0 0 -1)
testTransform(a3, c(0, 0, 0, -1))
summary(q1 <- powerTransform(lm(cbind(LoBD$I1L2, LoBD$I1L1) ~ pool, LoBD), family="bcnPower"))
testTransform(q1, c(.3, .8))

}
\keyword{ regression}% at least one, from doc/KEYWORDS
