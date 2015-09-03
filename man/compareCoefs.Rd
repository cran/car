\name{compareCoefs}
\alias{compareCoefs}
\title{
Print estimated coefficients and their standard errors in a table
for several regression models.
}
\description{
This simple function extracts estimates of regression parameters and
their standard errors from one or more models and prints them in a 
table.
}
\usage{
compareCoefs(..., se = TRUE, print=TRUE, digits = 3)
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
If \code{TRUE}, the default, show standard errors as well as estimates, if \code{FALSE},
show only estimates.
}
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
This function is used for its side-effect of printing the result.  It
returns a matrix of estimates and standard errors.
}
\references{
Fox, J. and Weisberg, S. (2011) \emph{An R Companion to Applied Regression},
Second Edition, Sage.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
mod1 <- lm(prestige ~ income + education, data=Duncan)
mod2 <- update(mod1, subset=-c(6,16))
mod3 <- update(mod1, . ~ . + type)
compareCoefs(mod1)
compareCoefs(mod1, mod2)
compareCoefs(mod1, mod2, mod3)
compareCoefs(mod1, mod2, se=FALSE)
}

\keyword{print }

