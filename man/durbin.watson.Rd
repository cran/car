\name{durbin.watson}
\alias{durbin.watson}
\alias{durbin.watson.lm}
\alias{durbin.watson.default}

\title{Durbin-Watson Test for Autocorrelated Errors}
\description{
  Computes residual autocorrelations and generalized Durbin-Watson
  statistics and their bootstrapped p-values.
}
\usage{
durbin.watson(model, ...)

durbin.watson.lm(model, max.lag=1, simulate=T, reps=1000,
    method=c("resample","normal"))

durbin.watson.default(residuals, max.lag=1)
}

\arguments{
  \item{model}{a linear-model object.}
  \item{max.lag}{maximum lag to which to compute residual autocorrelations
    and Durbin-Watson statistics.}
  \item{simulate}{if \code{TRUE} p-values will be estimated by bootstrapping.}
  \item{reps}{number of bootstrap replications.}
  \item{method}{bootstrap method: \code{"resample"} to resample from the observed
    residuals; \code{"normal"} to sample normally distributed errors with 0 mean
    and standard deviation equal to the standard error of the regression.}
  \item{residuals}{vector of residuals from a linear model.}
  \item{\dots}{arguments to be passed down to method functions.}
}

\value{
    Returns an object of type \code{"durbin.watson"}.
}

\references{ 
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
data(Hartnagel)
durbin.watson(lm(ftheft ~ tfr + partic + degrees + mtheft, data=Hartnagel))
## lag Autocorrelation D-W Statistic p-value
##   1         0.24419        1.4236   0.018
}

\keyword{regression}
\keyword{ts}
