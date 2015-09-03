\name{powerTransform}
\alias{powerTransform}
\alias{powerTransform.default}
\alias{powerTransform.lm}
\alias{powerTransform.formula}
\alias{estimateTransform}
\alias{estimateTransform.default}
\alias{estimateTransform.skewPower}

\title{Finding Univariate or Multivariate Power Transformations}
\description{
  \code{powerTransform} uses the maximum likelihood-like approach of Box and Cox (1964) to select a transformatiion of a univariate or multivariate response for normality, linearity and/or constant variance.  Available families are the default Box-Cox power family, and the Yeo-Johnson and skew power familes that may be useful when a response is not strictly positive.  \code{powerTransform} passes information to \code{estimateTransform}, and only the former will be of interest for most users.
}
\usage{
powerTransform(object, ...)

\S3method{powerTransform}{default}(object, family="bcPower", ...)

\S3method{powerTransform}{lm}(object, family="bcPower", ...)

\S3method{powerTransform}{formula}(object, data, subset, weights, na.action, family="bcPower",
  ...)
  
estimateTransform(X, Y, weights=NULL, family="bcPower", start=NULL,
         method="L-BFGS-B", ...)
         
\method{estimateTransform}{default}(X, Y, weights=NULL, family="bcPower", start=NULL, 
         method="L-BFGS-B", ...)
         
\method{estimateTransform}{skewPower}(X, Y, weights=NULL, ...)
}

\arguments{
  \item{object}{This can either be an object of class
\code{lm}, a formula, or a matrix or vector; see below.}
  \item{data}{A data frame or environment, as in \code{\link{lm}}.}
  \item{subset}{Case indices to be used, as in \code{\link{lm}}.}
  \item{weights}{Weights as in \code{\link{lm}}.}
  \item{na.action}{Missing value action, as in \sQuote{lm}.}
  \item{family}{The quoted name of a family of transformations.  The available options are \code{"bcPower"} the default for the Box-Cox power family; \code{"yjpower"} for the Yeo-Johnson family, and \code{"skewPower"} for the two-parameter skew power family.  The families are documented at \code{\link{bcPower}} and \code{\link{skewPower}}.}
  \item{...}{Additional arguments that are passed to \code{estimateTransform}
  which does the actual computing, or to the \code{\link{optim}} 
  function, which does the maximization.}
  \item{X}{A matrix or data.frame giving the \dQuote{right-side variables}, including a column of ones if the intercept is present.}
  \item{Y}{A vector or matrix or data.frame giving the \dQuote{left-side
variables.}}
  \item{start}{Starting values for the computations.  The default value of NULL is usually adequate.}
  \item{method}{The computing alogrithm used by \code{\link{optim}} for the
  maximization.  The default \code{"L-BFGS-B"} appears to work well.}
}

\details{The function \code{powerTransform} is used to estimate normalizing/linearizing/variance stabilizing transformations
of a univariate or a multivariate response in a linear regression.  For a univariate response,
a formula like \code{z~x1+x2+x3} will estimate a transformation for the response
\code{z} from a family of transformations indexed by one parameter for Box-Cox and Yeo-Johnson transformations,
or two parameters for the skew power family,
that makes the residuals from the regression of the transformed \code{z} on the predictors
as closed to normally distributed as possible.  

For a formula like \code{cbind(y1,y2,y3)~x1+x2+x3}, the three variables on
the left-side are all transformed, generally with different transformations
to make all the residuals as close to
normally distributed as possible.  This is not the same as three univariate transformations becuase the variables transformed are allowed to be correlated.  \code{cbind(y1,y2,y3)~1} would specify transformations
to multivariate normality with no predictors.  This generalizes the multivariate
power transformations suggested by Velilla (1993) by allowing for different
families of transformations, and by allowing for predictors.  Cook and Weisberg (1999)
and Weisberg (2014) suggest the usefulness of transforming
a set of predictors \code{z1, z2, z3} for multivariate normality and for transforming
for multivariate normality conditional on levels of a factor, which is equivalent
to setting the predictors to be indicator variables for that factor.  

Specifying the first argument as a vector, for example
\code{powerTransform(ais$LBM)}, is equivalent to
\code{powerTransform(LBM ~ 1, ais)}.  Similarly,
\code{powerTransform(cbind(ais$LBM, ais$SSF))}, where the first argument is a matrix
rather than a formula is equivalent to specification of a mulitvariate linear model
\code{powerTransform(cbind(LBM, SSF) ~ 1, ais)}. 

Three families of power transformations are available.
The Box-Cox pwoer family of power transformations,
\code{family="bcPower"}, 
equals \eqn{(U^{\lambda}-1)/\lambda}{(U^(lambda)-1)/lambda} 
for \eqn{\lambda}{lambda} \eqn{\neq}{not equal to} 0, and 
\eqn{\log(U)}{log(U)} if \eqn{\lambda =0}{lambda = 0}.  A scaled version of this transformation is used in computing with all the families to make the Jacobian of the transformation equal to 1.
  
If \code{family="yjPower"} then the Yeo-Johnson transformations are used.  
This is is Box-Cox transformation of \eqn{U+1} for nonnegative values, 
and of \eqn{|U|+1} with parameter \eqn{2-\lambda}{2-lambda} for \eqn{U}
negative.

If \code{family="skewPower"} then the skew power family of transformations suggested
by Hawkins and Weisberg (2015) is used.  This is a two-parameter family that would
generally be applied with a response with occasional negative values; see \code{\link{skewPower}}
for the details and examples.  This family has a power parameter \eqn{\lambda}{lambda} and a non-negative start parameter \eqn{\gamma}{gamma}, with \eqn{\gamma = 0}{gamma=0} equal to the Box-Cox transformation.

The function \code{\link{testTransform}} is used to obtain likelihood ratio
tests for
any specified value for the transformation parameter(s).  
  
}
\value{
An object of class \code{powerTransform} or class \code{skewpowerTransform} if \code{family="skewPower"} that
inherits from \code{powerTransfrom} is returned, including the components listed below.

Several methods are available for use with \code{powerTransform} objects.  The \code{coef} method returns 
the estimated transformation parameters, while \code{coef(object, round=TRUE)}
will
return the transformations rounded to nearby convenient values within 1.96
standard errors of the mle, if any exist.
The \code{vcov} function returns the estimated covariance matrix of the
estimated
transformation parameters.  A \code{print} method is used to
print the estimates and \code{summary} method provides more information including  likelihood ratio type
tests that all power parameters equal one and that all transformation
parameters equal zero, for log transformations, and for a convenient rounded value
not far from the mle.  In the case of the skew power family, these tests are based on the profile log-likelihood obtained by maximizing over the start parameter, thus treating the start as a nusiance parameter of lesser interest than the pwoer parameter.  \code{testTransform} can be called
directly to test any other value for \eqn{\lambda}{lambda} or for skew power \eqn{\lambda}{lambda} and \eqn{\gamma}{gamma}.  There is a \code{\link{plot.powerTransform}} method for plotting the transformed values, and also a \code{\link{contour.skewpowerTransform}} method to obtain a contour plot of the two-dimensional log-likelihood for the skew power parameters when the response in univariate.  Finally, the \code{\link{boxCox}} method can be used to plot the univariate log-likleihood for the Box-Cox or Yeo-Johnson power families, or the profile log-likelihood of each of the parameters in the skew power family.

The components of the returned object are
\item{value}{The value of the loglikelihood at the mle.}
\item{counts}{See \code{\link{optim}}.}
\item{convergence}{See \code{\link{optim}}.}
\item{message}{See \code{\link{optim}}.}
\item{hessian}{The hessian matrix.}
\item{start}{Starting values for the computations.}
\item{lambda}{The ml estimate for the power parameter}
\item{gamma}{The ml estimate for the start parameter for the skew power family}
\item{roundlam}{Convenient rounded values for the estimates.  These rounded
values will  often be the desirable transformations.}
\item{family}{The transformation family} 
\item{xqr}{QR decomposition of the predictor matrix.}
\item{y}{The responses to be transformed}
\item{x}{The predictors}
\item{weights}{The weights if weighted least squares.
}
}
\references{Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations.
\emph{Journal
of the Royal Statisistical Society, Series B}. 26 211-46.

Cook, R. D. and Weisberg, S. (1999)  \emph{Applied Regression Including
Computing
and Graphics}.  Wiley.

Fox, J. and Weisberg, S. (2011) 
\emph{An R Companion to Applied Regression}, Second Edition, Sage.

Hawkins, D. and Weisberg, S. (2015)
Combining the Box-Cox Power and Genralized Log Transformations to Accomodate Negative Responses,
submitted for publication.

Velilla, S. (1993)  A note on the multivariate Box-Cox transformation to 
normality.  \emph{Statistics and Probability Letters}, 17, 259-263.

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.

Yeo, I. and Johnson, R. (2000) A new family of
power transformations to improve normality or symmetry.  
\emph{Biometrika}, 87, 954-959.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{\code{\link{testTransform}},
\code{\link{optim}}, \code{\link{bcPower}}, \code{\link{skewPower}}, \code{\link{transform}}, \code{\link{boxCox}}. Documentation for
\code{\link{skewPower}} includes examples of the use of the skew power family.}
\examples{
# Box Cox Method, univariate
summary(p1 <- powerTransform(cycles ~ len + amp + load, Wool))

# fit linear model with transformed response:
coef(p1, round=TRUE)
summary(m1 <- lm(bcPower(cycles, p1$roundlam) ~ len + amp + load, Wool))

# Multivariate Box Cox
summary(powerTransform(cbind(len, ADT, trks, sigs1) ~ 1, Highway1))

# Multivariate transformation to normality within levels of 'hwy'
summary(a3 <- powerTransform(cbind(len, ADT, trks, sigs1) ~ hwy, Highway1))

# test lambda = (0 0 0 -1)
testTransform(a3, c(0, 0, 0, -1))

# save the rounded transformed values, plot them with a separate
# color for each highway type
transformedY <- bcPower(with(Highway1, cbind(len, ADT, trks, sigs1)),
                coef(a3, round=TRUE))
\dontrun{pairs(transformedY, col=as.numeric(Highway1$hwy)) }
}
\keyword{ regression}% at least one, from doc/KEYWORDS
