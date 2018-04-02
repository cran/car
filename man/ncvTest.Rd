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

\method{ncvTest}{lm}(model, var.formula, ...)

\method{ncvTest}{glm}(model, ...) # to report an error
}

\arguments{
  \item{model}{a weighted or unweighted linear model, produced by \code{lm}.}
  \item{var.formula}{a one-sided formula for the error variance; if omitted,
    the error variance depends on the fitted values.}
  \item{\dots}{arguments passed down to methods functions; not currently used.}
}
\details{
  This test is often called the Breusch-Pagan test; it was independently
  suggested with some extension by Cook and Weisberg (1983).
  
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
  
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
  
  Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}, Sandy Weisberg \email{sandy@umn.edu}}

\seealso{\code{\link{hccm}}, \code{\link[car]{spreadLevelPlot}} }

\examples{
ncvTest(lm(interlocks ~ assets + sector + nation, data=Ornstein))

ncvTest(lm(interlocks ~ assets + sector + nation, data=Ornstein), 
    ~ assets + sector + nation, data=Ornstein)
}

\keyword{htest}
\keyword{regression}
