\name{boxTidwell}
\alias{boxTidwell}
\alias{boxTidwell.formula}
\alias{boxTidwell.default}
\alias{print.boxTidwell}

\title{Box-Tidwell Transformations}

\description{
  Computes the Box-Tidwell power transformations of the predictors in a
  linear model.
}

\usage{
boxTidwell(y, ...)

\method{boxTidwell}{formula}(formula, other.x=NULL, data=NULL, subset, 
  na.action=getOption("na.action"), verbose=FALSE, tol=0.001, 
  max.iter=25, ...)

\method{boxTidwell}{default}(y, x1, x2=NULL, max.iter=25, tol=0.001, 
  verbose=FALSE, ...)
  
\method{print}{boxTidwell}(x, digits=getOption("digits") - 2, ...)
}

\arguments{
  \item{formula}{two-sided formula, the right-hand-side of which gives the
    predictors to be transformed.}
  \item{other.x}{one-sided formula giving the predictors that are \emph{not}
    candidates for transformation, including (e.g.) factors.}
  \item{data}{an optional data frame containing the variables in the model. 
    By default the variables are taken from the environment from which 
    \code{boxTidwell} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{na.action}{a function that indicates what should happen when the data contain \code{NA}s. 
    The default is set by the \code{na.action} setting of \code{options}.}
  \item{verbose}{if \code{TRUE} a record of iterations is printed; default is \code{FALSE}.}
  \item{tol}{if the maximum relative change in coefficients is less than \code{tol} then
    convergence is declared.}
  \item{max.iter}{maximum number of iterations.}
  \item{y}{response variable.}
  \item{x1}{matrix of predictors to transform.}
  \item{x2}{matrix of predictors that are \emph{not} candidates for transformation.}
  \item{\dots}{not for the user.}
  \item{x}{\code{boxTidwell} object.}
  \item{digits}{number of digits for rounding.}
}

\details{
  The maximum-likelihood estimates of the transformation parameters are computed
  by Box and Tidwell's (1962) method, which is usually more efficient than using
  a general nonlinear least-squares routine for this problem. Score tests for the
  transformations are also reported.
}

\value{
  an object of class \code{boxTidwell}, which is normally just printed.
}

\references{ 
  Box, G. E. P. and Tidwell, P. W. (1962)
  Transformation of the independent variables.
  \emph{Technometrics} \bold{4}, 531-550.

  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
boxTidwell(prestige ~ income + education, ~ type + poly(women, 2), data=Prestige)
}

\keyword{regression}
