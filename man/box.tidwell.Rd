\name{box.tidwell}
\alias{box.tidwell}
\alias{box.tidwell.formula}
\alias{box.tidwell.default}
\alias{print.box.tidwell}

\title{Box-Tidwell Transformations}

\description{
  Computes the Box-Tidwell power transformations of the predictors in a
  linear model.
}

\usage{
box.tidwell(y, ...)

\method{box.tidwell}{formula}(formula, other.x=NULL, data=NULL, subset, 
  na.action=options()$na.action, verbose=FALSE, tol=0.001, 
  max.iter=25, ...)

\method{box.tidwell}{default}(y, x1, x2=NULL, max.iter=25, tol=0.001, 
  verbose=FALSE, ...)
  
\method{print}{box.tidwell}(x, digits, ...)
}

\arguments{
  \item{formula}{two-sided formula, the right-hand-side of which gives the
    predictors to be transformed.}
  \item{other.x}{one-sided formula giving the predictors that are \emph{not}
    candidates for transformation, including (e.g.) factors.}
  \item{data}{an optional data frame containing the variables in the model. 
    By default the variables are taken from the environment from which 
    \code{box.tidwell} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{na.action}{a function that indicates what should happen when the data contain \code{NA}s. 
    The default is set by the \code{na.action} setting of \code{options}.}
  \item{verbose}{if \code{TRUE} a record of iterations is printed.}
  \item{tol}{if maximum relative change in coefficients is less than \code{tol} then
    convergence is declared.}
  \item{max.iter}{maximum number of iterations.}
  \item{y}{response variable.}
  \item{x1}{matrix of predictors to transform.}
  \item{x2}{matrix of predictors that are \emph{not} candidates for transformation.}
  \item{\dots}{not for the user.}
  \item{x}{\code{box.tidwell} object.}
  \item{digits}{number of digits for rounding.}
}

\details{
  The maximum-likelihood estimates of the transformation parameters are computed
  by Box and Tidwell's (1962) method, which is usually more efficient than using
  a general nonlinear least-squares routine for this problem. Score tests for the
  transformations are also reported.
}

\value{
  an object of class \code{box.tidwell}, which is normally just printed.
}

\references{ 
  Box, G. E. P. and Tidwell, P. W. (1962)
  Transformation of the independent variables.
  \emph{Technometrics} \bold{4}, 531-550.

  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage. 
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
box.tidwell(prestige~income+education, ~ poly(women,2), data=Prestige)
##                   income education
## Initial Power   -0.91030   2.24354
## Score Statistic -5.30129   2.40556
## p-value          0.00000   0.01615
## MLE of Power    -0.03777   2.19283
}

\keyword{regression}
