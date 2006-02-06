\name{ncv.test}
\alias{ncv.test}
\alias{ncv.test.lm}
\alias{ncv.test.glm}

\title{Score Test for Non-Constant Error Variance}
\description{
  Computes a score test of the hypothesis of constant error variance
  against the alternative that the error variance changes with the
  level of the response (fitted values), or with a linear combination
  of predictors.
}
\usage{
ncv.test(model, ...)

\method{ncv.test}{lm}(model, var.formula, data=NULL, subset, na.action, ...)

\method{ncv.test}{glm}(model, ...)
}

\arguments{
  \item{model}{a weighted or unweighted linear model, produced by \code{lm}.}
  \item{var.formula}{a one-sided formula for the error variance; if omitted,
    the error variance depends on the fitted values.}
  \item{data}{an optional data frame containing the variables in the model. 
    By default the variables are taken from the environment from which 
    \code{ncv.test} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{na.action}{a function that indicates what should happen when the data contain \code{NA}s. 
    The default is set by the \code{na.action} setting of \code{options}.}
  \item{\dots}{arguments passed down to methods functions.}
}
\details{
  This test is often called the Breusch-Pagan test; it was independently
  suggested by Cook and Weisberg (1983).
  
  \code{ncv.test.glm} is a dummy function to generate an error when a \code{glm}
  model is used.
}
\value{
  The function returns a \code{chisq.test} object, which is usually just printed.
}
\references{
  Breusch, T. S. and Pagan, A. R. (1979)
  A simple test for heteroscedasticity and random coefficient variation.
  \emph{Econometrica} \bold{47}, 1287--1294.
  
  Cook, R. D. and Weisberg, S. (1983)
  Diagnostics for heteroscedasticity in regression.
  \emph{Biometrika} \bold{70}, 1--10.
  
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage. 
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{hccm}}, \code{\link{spread.level.plot}} }

\examples{
ncv.test(lm(interlocks~assets+sector+nation, data=Ornstein))
## Non-constant Variance Score Test 
## Variance formula: ~ fitted.values 
## Chisquare = 46.98537    Df = 1     p = 7.151835e-12 
ncv.test(lm(interlocks~assets+sector+nation, data=Ornstein), 
    ~ assets+sector+nation, data=Ornstein)
## Non-constant Variance Score Test 
## Variance formula: ~ assets + sector + nation 
## Chisquare = 74.73535    Df = 13     p = 1.066320e-10 
}

\keyword{htest}
\keyword{regression}
