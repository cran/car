\name{coefTable}
\alias{coefTable}
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
coefTable(..., se = TRUE, digits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
The names of one or more regression models.  These may be of class \code{lm}, \code{glm},
\code{nlm}, or any other regression method for which the functions \code{coef} and 
\code{vcov} return appropriate values.
}
  \item{se}{
If \code{TRUE}, the default, show standard errors as well as estimates, if \code{FALSE},
show only estimates.
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
coefTable(mod1)
coefTable(mod1, mod2)
coefTable(mod1, mod2, mod3)
coefTable(mod1, mod2, se=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print }

