\name{durbin.watson}
\alias{durbin.watson}
\alias{durbin.watson.lm}
\alias{durbin.watson.default}
\alias{print.durbin.watson}

\title{Durbin-Watson Test for Autocorrelated Errors}
\description{
  Computes residual autocorrelations and generalized Durbin-Watson
  statistics and their bootstrapped p-values.
}
\usage{
durbin.watson(model, ...)

\method{durbin.watson}{lm}(model, max.lag=1, simulate=TRUE, reps=1000,
    method=c("resample","normal"),
    alternative=c("two.sided", "positive", "negative"), ...)

\method{durbin.watson}{default}(model, max.lag=1, ...)

\method{print}{durbin.watson}(x, ...)
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
  \item{\dots}{arguments to be passed down to method functions.}
  \item{x}{\code{durbin.watson} object.}
}

\value{
    Returns an object of type \code{"durbin.watson"}.
}

\note{
p-values are available only from the \code{lm} method.
}

\references{ 
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
durbin.watson(lm(fconvict ~ tfr + partic + degrees + mconvict, data=Hartnagel))
##  lag Autocorrelation D-W Statistic p-value 
##    1        0.688345     0.6168636       0 
##  Alternative hypothesis: rho != 0
}

\keyword{regression}
\keyword{ts}
