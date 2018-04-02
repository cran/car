\name{compareCoefs}
\alias{compareCoefs}
\title{
Print estimated coefficients and their standard errors in a table
for several regression models.
}
\description{
This function extracts estimates of regression parameters and
their standard errors from one or more models and prints them in a
table.
}
\usage{
compareCoefs(..., se = TRUE, zvals = FALSE, pvals = FALSE, vcov.,
    print = TRUE, digits = 3)
}

\arguments{
  \item{\dots}{
One or more regression-model objects.  These may be of class \code{lm}, \code{glm},
\code{nlm}, or any other regression method for which the functions \code{coef} and
\code{vcov} return appropriate values, or if the object inherits from the \code{mer}
class created by the \code{lme4} package or \code{lme} in the \code{nlme}
package.
}
  \item{se}{
If \code{TRUE}, the default, show standard errors as well as estimates.
}
\item{zvals}{If \code{TRUE} (the default is \code{FALSE}), print Wald statistics, the
ratio of each coefficient to its standard error.}
\item{pvals}{If code{TRUE} (the default is \code{FALSE}), print two-sided p-values
from the standard normal distribution corresponding to the Wald statistics.}
\item{vcov.}{an optional argument, specifying a function to be applied to all of the models, returning a coefficient covariance matrix
for each, or a list with one element for each model, with each element either containing a function to
be applied to the corresponding model or a coefficient covariance matrix for that model. If omitted, \code{vcov} is applied to
each model.}
\item{print}{
If \code{TRUE}, the defualt, the results are printed in a nice format
using \code{\link{printCoefmat}}.  If
\code{FALSE}, the results are returned as a matrix
}
  \item{digits}{
Passed to the \code{\link{printCoefmat}} function for printing the result.
}
}

\value{
This function is mainly used for its side-effect of printing the result.  It
also  invisibly returns a matrix of estimates, standard errors, Wald statistics, and p-values.
}
\references{
Fox, J. and Weisberg, S. (2019) \emph{An R Companion to Applied Regression},
Third Edition, Sage.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
mod1 <- lm(prestige ~ income + education, data=Duncan)
mod2 <- update(mod1, subset=-c(6,16))
mod3 <- update(mod1, . ~ . + type)
mod4 <- update(mod1, . ~ . + I(income + education)) # aliased coef.
compareCoefs(mod1)
compareCoefs(mod1, mod2, mod4)
compareCoefs(mod1, mod2, mod3, zvals=TRUE, pvals=TRUE)
compareCoefs(mod1, mod2, se=FALSE)
compareCoefs(mod1, mod1, vcov.=list(vcov, hccm))
}

\keyword{print }

