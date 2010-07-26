\name{deltaMethod}
\alias{deltaMethod}
\alias{deltaMethod.default}
\alias{deltaMethod.lm} 
\alias{deltaMethod.nls}
\alias{deltaMethod.multinom} 
\alias{deltaMethod.polr} 
\alias{deltaMethod.survreg}
\alias{deltaMethod.coxph}

\title{Estimate and Standard Error of a Nonlinear Function of Estimated Regression 
Coefficients}
\description{
\code{deltaMethod} is a generic function that uses the delta method to get a 
first-order approximate 
standard error for a nonlinear function of a vector of random variables
with known or estimated covariance matrix.
}
\usage{
deltaMethod(object, ...)

\method{deltaMethod}{default}(object, g, vcov., func=g, ...)
\method{deltaMethod}{lm} (object, g, vcov.=vcov, parameterPrefix="b", ...)
\method{deltaMethod}{nls}(object, g, vcov.=vcov, ...)
\method{deltaMethod}{multinom} (object, g, vcov.=vcov, parameterPrefix="b", ...)
\method{deltaMethod}{polr} (object, g, vcov.=vcov, ...)
\method{deltaMethod}{survreg} (object, g, vcov.=vcov, ...)
\method{deltaMethod}{coxph} (object, g, vcov.=vcov, ...)

}
\arguments{
  \item{object}{For the default method, \code{object} is a named vector of \code{p} 
  elements.  This means that the call \code{names(object)} would return a list
  of \code{p} character strings that are the names of the elements of
  \code{object}.  For the other methods, \code{object} is a 
  regression object for which \code{coef(object)} returns a vector of parameter
  estimates.}
  \item{g}{A quoted string that is the function of the parameter estimates to be
  evaluated; see the details below. }
  \item{vcov.}{The (estimated) covariance matrix of the coefficient
  estimates.  For the default method, this argument is required.  For all
  other methods, this argument must either provide the estimated covariance
  matrix or a function that when applied
  to \code{object} returns a covariance matrix.  The default is 
  to use the function \code{vcov}.}
  \item{func}{A quoted string used to annotate output.  The default of
  \code{func = g} is usually appropriate.}
  \item{parameterPrefix}{Typically a single letter with default
  \code{"b"} giving the prefix of the
  names of the parameter names used in the argument \code{g} for some
  regression models; see details.}
  \item{...}{Additional arguments; not currently used.}
}

\details{
Suppose \eqn{x} is a random vector of length \eqn{p} that is at least approximately 
normally distributed with mean \eqn{\beta} and estimated covariance matrix 
\eqn{C}.   Then any function \eqn{g(\beta)} of \eqn{\beta}, is estimated by
\eqn{g(x)}, which is in large samples normally distributed with mean \eqn{g(\beta)}
and estimated variance \eqn{h'Ch}, where \eqn{h} is the first derivative of
\eqn{g(\beta)} with respect to \eqn{\beta} evaluated at \eqn{x}.  This function
returns both \eqn{g(x)} and its standard error, the square root of the estimated
variance.

The default method requires that you provide \eqn{x} in the argument
\code{object}, \eqn{C} in the argument \code{vcov.}, and a text expression
in argument \code{g} that when evaluated gives the function \eqn{g}.  

Since
the delta method is often applied to functions of regression parameter 
estimates, the argument \code{object} may be the name of a regression
object from which the vector \eqn{x} will be taken from \code{coef(object)},
and \eqn{C} will be taken from \code{vcov(object)} unless you provide
some other estimate of variance, for example, using a sandwich estimator.  Methods
have been provided for several common regression models.   

For regression models for which methods are not provided, you must extract
the named vector of coefficient estimates and and estimate of its covariance
matrix and then apply the default \code{deltaMethod} function. 

In the argument \code{g} you must provide a quoted character string
that gives the function of interest.  For example, if you set 
\code{m2 <- lm(Y ~ X1 + X2)}, then \code{deltaMethod(m2,"X1/X2")} applies the
delta method to the ratio of the coefficient estimates for \code{X1} and
\code{X2}.  For the product of the interecpt and the coefficient of \code{X2},
use \code{deltaMethod(m2, "(Intercept)*X2")}, since the name of the interecpt
parameter estimate is \code{(Intercept)}. 

For models of type \code{lm}, \code{glm} and \code{polr}, you can replace the
names in the coefficient vector by \code{b0}, \code{b1}, ..., \code{bp}; you can
also change the prefix \code{b} using the \code{parameterPrefix} argument.
For \code{multinom} objects you 
can only use the parameter names starting with the prefix, not the names of
the corresponding regressors. 

For nonlinear regression objects of type nls, the call \code{coef(object)} 
returns the estimated
coefficient vectors with names corresponding to parameter names.  
For example,           
\code{m2 <- nls(y ~ theta/(1 + gamma * x), start = list(theta=2, gamma=3))} will
have parameters named  \code{c("theta", "gamma")}.   
In many other familiar regression methods, such as lm and glm, the names of
the coefficient estimates are the corresponding variable names, not
parameter names.    
}

\value{
A data.frame with two components
named \code{Estimate} for the estimate, \code{SE} for its standard error.
The value of \code{g} is given as a row label.  
}

\seealso{First derivatives of \code{g} are computed using symbolic differentiation
by the function \code{\link{D}}.}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.

  S. Weisberg (2005) \emph{Applied
  Linear Regression}, Third Edition, Wiley, Section 6.1.2.
}

\author{Sanford Weisberg, \email{sandy@stat.umn.edu},
and John Fox \email{jfox@mcmaster.ca}}

\examples{
m1 <- lm(time ~ t1 + t2, data = Transact) 
deltaMethod(m1, "b1/b2") # ratio of coefficients
deltaMethod(m1, "t1/t2") # use names of preds. rather than coefs.
deltaMethod(m1, "t1/t2", vcov=hccm) # use hccm function to est. vars.
# The next example calls the default method by extracting the
# vector of estimates and covariance matrix explicitly
deltaMethod(coef(m1), "t1/t2", vcov.=vcov(m1)) 
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ regression }
