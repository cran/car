\name{durbinWatsonTest}
\alias{durbinWatsonTest}
\alias{dwt}
\alias{durbinWatsonTest.lm}
\alias{durbinWatsonTest.default}
\alias{print.durbinWatsonTest}

\title{Durbin-Watson Test for Autocorrelated Errors}
\description{
  Computes residual autocorrelations and generalized Durbin-Watson
  statistics and their bootstrapped p-values. \code{dwt} is an
  abbreviation for \code{durbinWatsonTest}.
}
\usage{
durbinWatsonTest(model, ...)

dwt(...)

\method{durbinWatsonTest}{lm}(model, max.lag=1, simulate=TRUE, reps=1000,
    method=c("resample","normal"),
    alternative=c("two.sided", "positive", "negative"), ...)

\method{durbinWatsonTest}{default}(model, max.lag=1, ...)

\method{print}{durbinWatsonTest}(x, ...)
}

\arguments{
  \item{model}{a linear-model object, or a vector of residuals from a linear model.}
  \item{max.lag}{maximum lag to which to compute residual autocorrelations
    and Durbin-Watson statistics.}
  \item{simulate}{if \code{TRUE} p-values will be estimated by bootstrapping.}
  \item{reps}{number of bootstrap replications.}
  \item{method}{bootstrap method: \code{"resample"} to resample from the observed
    residuals; \code{"normal"} to sample normally distributed errors with 0 mean
    and standard deviation equal to the standard error of the regression.}
  \item{alternative}{sign of autocorrelation in alternative hypothesis; specify
    only if \code{max.lag = 1}; if \code{max.lag > 1}, then \code{alternative} is
    taken to be \code{"two.sided"}.}
  \item{\dots}{arguments to be passed down.}
  \item{x}{\code{durbinWatsonTest} object.}
}

\value{
    Returns an object of type \code{"durbinWatsonTest"}.
}

\note{
p-values are available only from the \code{lm} method.
}

\references{ 
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models}, Third Edition. Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
durbinWatsonTest(lm(fconvict ~ tfr + partic + degrees + mconvict, data=Hartnagel))
}

\keyword{regression}
\keyword{ts}
