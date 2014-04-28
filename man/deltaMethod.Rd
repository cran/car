\name{deltaMethod}
\alias{deltaMethod}
\alias{deltaMethod.default}
\alias{deltaMethod.lm} 
\alias{deltaMethod.nls}
\alias{deltaMethod.multinom} 
\alias{deltaMethod.polr} 
\alias{deltaMethod.survreg}
\alias{deltaMethod.coxph}
\alias{deltaMethod.mer}
\alias{deltaMethod.merMod}
\alias{deltaMethod.lme}
\alias{deltaMethod.lmList}

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

\method{deltaMethod}{default}(object, g, vcov., func=g, constants, ...)
\method{deltaMethod}{lm} (object, g, vcov.=vcov, 
           parameterNames=names(coef(object)), ...)
\method{deltaMethod}{nls}(object, g, vcov.=vcov, ...)
\method{deltaMethod}{multinom} (object, g, vcov. = vcov, 
           parameterNames = if (is.matrix(coef(object))) 
           colnames(coef(object)) else names(coef(object)), ...) 
\method{deltaMethod}{polr} (object, g, vcov.=vcov, ...)
\method{deltaMethod}{survreg} (object, g, vcov. = vcov, 
           parameterNames = names(coef(object)), ...)
\method{deltaMethod}{coxph} (object, g, vcov. = vcov, 
           parameterNames = names(coef(object)), ...)
\method{deltaMethod}{mer} (object, g, vcov. = vcov,
           parameterNames = names(fixef(object)), ...)
\method{deltaMethod}{merMod} (object, g, vcov. = vcov,
           parameterNames = names(fixef(object)), ...)
\method{deltaMethod}{lme} (object, g, vcov. = vcov,
           parameterNames = names(fixef(object)), ...)
\method{deltaMethod}{lmList} (object, g,  ...)

}
\arguments{
  \item{object}{For the default method, \code{object} is either (1) a vector of \code{p} 
  named elements, so \code{names(object)} returns a list
  of \code{p} character strings that are the names of the elements of
  \code{object}; or (2) a model object for which there are \code{\link[stats]{coef}} and \code{\link[stats]{vcov}} methods,
  and for which the named coefficient vector returned by \code{coef} is asymptotically normally distributed with asymptotic
  covariance matrix returned by \code{vcov}.
  For the other methods, \code{object} is a 
  regression object for which \code{coef(object)} or \code{fixef(object)} returns a vector of parameter
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
  \item{parameterNames}{A character vector of length \code{p} that gives the 
  names of the parameters in the same order as they appear in the vector of
  estimates.  This argument will be useful if some of the names in the
  vector of estimates include special characters, like \code{I(x2^2)}, or 
  \code{x1:x2} that will confuse the numerical differentiation function.  See
  details below.}
  \item{constants}{This argument is a named vector whose elements are constants that
  are used in the \code{f} argument.  This is needed only when the function is
  called from within another function to comply to R scoping rules.  See
  example below.}
  \item{...}{Used to pass arguments to the generic method.}
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
in argument \code{g} that when evaluated gives the function \eqn{g}.  The
call \code{names(object)} must return the names of the elements of \code{x}
that are used in the expression \code{g}.  

Since
the delta method is often applied to functions of regression parameter 
estimates, the argument \code{object} may be the name of a regression
object from which the the estimates and their estimated variance matrix can
be extracted.  In most regression models, estimates are returned by the
\code{coef(object)} and the variance matrix from \code{vcov(object)}.
You can provide an alternative function for computing the sample variance
matrix, for example to use a sandwich estimator.
         
For mixed models using \code{lme4} or \code{nlme}, the coefficient estimates
are returned by the \code{fixef} function, while for \code{multinom},
\code{lmList} and \code{nlsList} coefficient estimates are returned by 
\code{coef} as a matrix.  Methods for these models are provided to get the
correct estimates and variance matrix.

The argument \code{g} must be a quoted character string
that gives the function of interest.  For example, if you set 
\code{m2 <- lm(Y ~ X1 + X2 + X1:X2)}, then \code{deltaMethod(m2,"X1/X2")} applies the
delta method to the ratio of the coefficient estimates for \code{X1} and
\code{X2}.  The argument \code{g} can consist of constants and names 
associated with the elements of the vector of coefficient estimates.
       
In some cases the names may include characters including such as the colon 
\code{:} used in interactions, or mathematical symbols like \code{+} or 
\code{-} signs that would confuse the 
function that computes numerical derivatives, and for this case you can replace
the names of the estimates with the \code{parameterNames} argument.  For
example, the ratio of the
\code{X2} main effect to the interaction term could be computed using
\code{deltaMethod(m2, "b1/b3", parameterNames=c("b0", "b1", "b2", "b3"))}.
The name \dQuote{\code{(Intercept)}} used for the intercept in linear and generalized
linear models is an exception, and it will be correctly interpreted by
\code{deltaMethod}.

For \code{multinom} objects, the \code{coef} function returns a matrix of
coefficients, with each row giving the estimates for comparisons of one category
to the baseline.  The \code{deltaMethod} function applies the delta method to
each row of this matrix.  Similarly, for \code{lmList} and \code{nlsList}
objects, the delta method is computed for each element of the list of
models fit. 

For nonlinear regression objects of type nls, the call \code{coef(object)} 
returns the estimated
coefficient vectors with names corresponding to parameter names.  
For example,           
\code{m2 <- nls(y ~ theta/(1 + gamma * x), start = list(theta=2, gamma=3))} will
have parameters named  \code{c("theta", "gamma")}.   
In many other familiar regression methods, such as lm and glm, the names of
the coefficient estimates are the corresponding variable names, not
parameter names. 

For mixed-effects models fit with \code{lmer} and \code{nlmer} from the 
\code{lme4} package or \code{lme} and \code{nlme} from the \code{nlme} package, 
only fixed-effect coefficients are considered.

For regression models for which methods are not provided, you can extract
the named vector of coefficient estimates and and estimate of its covariance
matrix and then apply the default \code{deltaMethod} function. 

Earlier versions of \code{deltaMethod} included an argument 
\code{parameterPrefix} that implemented the same functionality as the
\code{parameterNames} argument, but it caused several unintended bugs that
were not easily fixed without the change in syntax.  
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

  Weisberg, S. (2014) \emph{Applied
  Linear Regression}, Fourth Edition, Wiley, Section 6.1.2.
}

\author{Sanford Weisberg, \email{sandy@umn.edu},
and John Fox \email{jfox@mcmaster.ca}}

\examples{
m1 <- lm(time ~ t1 + t2, data = Transact) 
deltaMethod(m1, "b1/b2", parameterNames= paste("b", 0:2, sep="")) 
deltaMethod(m1, "t1/t2") # use names of preds. rather than coefs.
deltaMethod(m1, "t1/t2", vcov=hccm) # use hccm function to est. vars.
# to get the SE of 1/intercept, rename coefficients
deltaMethod(m1, "1/b0", parameterNames= paste("b", 0:2, sep=""))
# The next example calls the default method by extracting the
# vector of estimates and covariance matrix explicitly
deltaMethod(coef(m1), "t1/t2", vcov.=vcov(m1))
# the following works:
a <- 5
deltaMethod(m1, "(t1 + a)/t2")
# ...but embedded in a function this will fail
f1 <- function(mod, ...) {
 z <- 3
 deltaMethod(m1, "(t1+z)/t2", ...)
 }
\dontrun{f1(m1)}
# if z is defined globally f1 could even return the wrong answer.
# the following function works
f2 <- function(mod, ...) {
 deltaMethod(m1, "(t1+z)/t2", ...)
 }
f2(m1, constants=c(z=3))
# as does
f3 <- function(mod) {
 a <- 3
 deltaMethod(m1, "(t1+z)/t2", constants=c(z=a))
 }
f3(m1)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ regression }
