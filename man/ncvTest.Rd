\name{ncvTest}
\alias{ncvTest}
\alias{ncvTest.lm}
\alias{ncvTest.glm}

\title{Score Test for Non-Constant Error Variance}
\description{
  Computes a score test of the hypothesis of constant error variance
  against the alternative that the error variance changes with the
  level of the response (fitted values), or with a linear combination
  of predictors.
}
\usage{
ncvTest(model, ...)

\method{ncvTest}{lm}(model, var.formula, data=NULL, subset, na.action, ...)

\method{ncvTest}{glm}(model, ...) # to report an error
}

\arguments{
  \item{model}{a weighted or unweighted linear model, produced by \code{lm}.}
  \item{var.formula}{a one-sided formula for the error variance; if omitted,
    the error variance depends on the fitted values.}
  \item{data}{an optional data frame containing the variables in the model. 
    By default the variables are taken from the environment from which 
    \code{ncvTest} is called. The \code{data} argument may therefore need to be specified even when
    the \code{data} argument was specified in the call to \code{lm} when the model was fit 
    (see the second example below).}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{na.action}{a function that indicates what should happen when the data contain \code{NA}s. 
    The default is set by the \code{na.action} setting of \code{options}.}
  \item{\dots}{arguments passed down to methods functions.}
}
\details{
  This test is often called the Breusch-Pagan test; it was independently
  suggested by Cook and Weisberg (1983).
  
  \code{ncvTest.glm} is a dummy function to generate an error when a \code{glm}
  model is used.
}
\value{
  The function returns a \code{chisqTest} object, which is usually just printed.
}
\references{
  Breusch, T. S. and Pagan, A. R. (1979)
  A simple test for heteroscedasticity and random coefficient variation.
  \emph{Econometrica} \bold{47}, 1287--1294.
  
  Cook, R. D. and Weisberg, S. (1983)
  Diagnostics for heteroscedasticity in regression.
  \emph{Biometrika} \bold{70}, 1--10.
  
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Weisberg, S. (2005) \emph{Applied Linear Regression}, Third Edition, Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{hccm}}, \code{\link[car]{spreadLevelPlot}} }

\examples{
ncvTest(lm(interlocks ~ assets + sector + nation, data=Ornstein))

ncvTest(lm(interlocks ~ assets + sector + nation, data=Ornstein), 
    ~ assets + sector + nation, data=Ornstein)
}

\keyword{htest}
\keyword{regression}
