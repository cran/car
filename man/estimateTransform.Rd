\name{estimateTransform}
\alias{estimateTransform}

\title{Finding Univariate or Multivariate Power Transformations}
\description{
  \code{estimateTransform} computes members of families of transformations
indexed by one
parameter, the Box-Cox power family, or the Yeo and Johnson (2000) family, or
the
basic power family, interpreting zero power as logarithmic. 
The family can be modified to have Jacobian one, or not, except for the basic
power family.  Most users will use the function \code{\link{powerTransform}}, which
is a front-end for this function.
}
\usage{
estimateTransform(X, Y, weights=NULL, family="bcPower", start=NULL,
         method="L-BFGS-B", ...)
}

\arguments{
  \item{X}{A matrix or data.frame giving the \dQuote{right-side variables}.}
  \item{Y}{A vector or matrix or data.frame giving the \dQuote{left-side
variables.}}
  \item{weights}{Weights as in \code{lm}.}
  \item{family}{The transformation family to use.  This is the quoted name of a
  function for computing the transformed values.  The default is
\code{bcPower}
  for the Box-Cox power family and the most likely alternative is
  \code{yjPower} for the Yeo-Johnson family of transformations.}
  \item{start}{Starting values for the computations.  It is usually adequate
  to leave this at its default value of NULL.}
  \item{method}{The computing alogrithm used by \code{\link{optim}} for the
  maximization.  The default \code{"L-BFGS-B"} appears to work well.}
  \item{...}{Additional arguments that are passed to the \code{\link{optim}} 
  function that does the maximization.  Needed only if there are convergence problems.}
}

\details{See the documentation for the function \code{\link{powerTransform}}.}

\value{An object of class \code{powerTransform} with components
\item{value}{The value of the loglikelihood at the mle.}
\item{counts}{See \code{\link{optim}}.}
\item{convergence}{See \code{\link{optim}}.}
\item{message}{See \code{\link{optim}}.}
\item{hessian}{The hessian matrix.}
\item{start}{Starting values for the computations.}
\item{lambda}{The ml estimate}
\item{roundlam}{Convenient rounded values for the estimates.  These rounded
values will  often be the desirable transformations.}
\item{family}{The transformation family} 
\item{xqr}{QR decomposition of the predictor matrix.}
\item{y}{The responses to be transformed}
\item{x}{The predictors}
\item{weights}{The weights if weighted least squares.}
}
\references{Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations. \emph{Journal
of the Royal Statisistical Society, Series B}. 26 211-46.

Cook, R. D. and Weisberg, S. (1999)  \emph{Applied Regression Including Computing
and Graphics}.  Wiley.

Fox, J. and Weisberg, S. (2011) 
\emph{An R Companion to Applied Regression}, Second Edition, Sage.

Velilla, S. (1993)  A note on the multivariate Box-Cox transformation to 
normality.  \emph{Statistics and Probability Letters}, 17, 259-263.

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.

Yeo, I. and Johnson, R. (2000) A new family of
power transformations to improve normality or symmetry.  
\emph{Biometrika}, 87, 954-959.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{\code{\link{powerTransform}}, \code{\link{testTransform}},
\code{\link{optim}}. }
\examples{
data(trees,package="MASS")
summary(out1 <- powerTransform(Volume~log(Height)+log(Girth),trees))
# multivariate transformation:
summary(out2 <- powerTransform(cbind(Volume,Height,Girth)~1,trees))
testTransform(out2,c(0,1,0))
# same transformations, but use lm objects
m1 <- lm(Volume~log(Height)+log(Girth),trees)
(out3 <- powerTransform(m1))
# update the lm model with the transformed response
update(m1,basicPower(out3$y,out3$roundlam)~.)
}
\keyword{ regression}% at least one, from doc/KEYWORDS
