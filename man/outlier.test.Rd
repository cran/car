\name{outlier.test}
\alias{outlier.test}
\alias{outlier.test.lm}
\alias{outlier.test.glm}
\alias{print.outlier.test}

\title{Bonferroni Outlier Test}

\description{
  Reports the Bonferroni p-value for the most extreme observation.
  At present, there are methods for studentized residuals in
  linear and generalized linear models.
}

\usage{
outlier.test(model, ...)

\method{outlier.test}{lm}(model, labels=names(rstud), ...)

\method{outlier.test}{glm}(model, labels=names(rstud), ...)

\method{print}{outlier.test}(x, digits=options("digits")[[1]], ...)
}

\arguments{
  \item{model}{a suitable model object.}
  \item{labels}{an optional vector of observation names.}
  \item{...}{arguments passed down to methods functions.}
  \item{x}{\code{outlier.test} object.}
  \item{digits}{number of digits for printed output.}
}

\details{
  For a linear model, the p-value reported is for the largest absolute
  studentized residual, using the \eqn{t} distribution with degrees of
  freedom one less than the residual df for the model. For a generalized
  linear model, the largest absolute studentized residual is also used,
  but with the standard-normal distribution. The Bonferroni
  adjustment multiplies the usual two-sided p-value by the number of 
  observations.
}

\value{
  an object of class \code{outlier.test}, which is normally just
  printed.
}

\references{
  Cook, R. D. and Weisberg, S. (1984)
  \emph{Residuals and Influence in Regression.} Wiley.
  
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
  
  Williams, D. A. (1987)
  Generalized linear model diagnostics using the deviance and single
  case deletions. \emph{Applied Statistics} \bold{36}, 181--191.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
outlier.test(lm(prestige~income+education, data=Duncan))
##  max|rstudent| df unadjusted p Bonferroni p
##       3.134519 41  0.003177202    0.1429741
## 
## Observation: minister 
}

\keyword{regression}
\keyword{htest}
