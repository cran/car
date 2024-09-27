\name{poTest}
\alias{poTest}
\alias{poTest.polr}
\alias{print.poTest}

\title{Test for Proportional Odds in the Proportional-Odds Logistic-Regression Model}

\description{
The \code{poTest} function implements tests proposed by Brant (1990) for proportional odds for
logistic models fit by the \code{\link[MASS]{polr}} function in the MASS package.
}

\usage{
poTest(model, ...)
\method{poTest}{polr}(model, ...)
\method{print}{poTest}(x, digits=3, ...)
}

\arguments{
  \item{model}{a proptional-odds logit model fit by \code{\link[MASS]{polr}}.}
  \item{x}{an object produced by \code{poTest}.}
  \item{digits}{number of significant digits to print.}
  \item{\dots}{ignored.}
}

\value{
    \code{poTest} returns an object meant to be printed showing the results of the tests.
}

\references{
    R. Brant, "Assessing proportionality in the proportional odds model for ordinal logistic regression."
    Biometrics 46: 1171--1178, 1990.
}
\author{
John Fox \email{jfox@mcmaster.ca}
}

\examples{
if (require("MASS")){
    .W <- Womenlf
    .W$partic <- factor(.W$partic, levels=c("not.work", "parttime", "fulltime"))
    poTest(polr(partic ~ hincome + children + region, data=.W))
}
}

\keyword{models}
\keyword{htest}
