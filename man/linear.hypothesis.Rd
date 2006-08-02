\name{linear.hypothesis}
\alias{linear.hypothesis}
\alias{lht}
\alias{linear.hypothesis.lm}
\alias{linear.hypothesis.glm}
\alias{linear.hypothesis.default}

\title{Test Linear Hypothesis}

\description{
Generic function for testing a linear hypothesis, and methods
for linear models, generalized linear models, and other models
that have methods for \code{coef} and \code{vcov}.
}

\usage{
linear.hypothesis(model, ...)

lht(model, ...)

\method{linear.hypothesis}{default}(model, hypothesis.matrix, rhs=NULL, 
    test=c("Chisq", "F"), vcov.=NULL, verbose=FALSE, ...)  

\method{linear.hypothesis}{lm}(model, hypothesis.matrix, rhs=NULL,
    test=c("F", "Chisq"), vcov.=NULL, white.adjust=FALSE, ...)

\method{linear.hypothesis}{glm}(model,  ...)  
}

\arguments{
  \item{model}{fitted model object. The default method works for models
   for which the estimated parameters can be retrieved by \code{coef} and
   the corresponding estimated covariance matrix by \code{vcov}. See the
   \emph{Details} for more information.}
  \item{hypothesis.matrix}{matrix (or vector) giving linear combinations
    of coefficients by rows, or a character vector giving the hypothesis
    in symbolic form (see \emph{Details}).}
  \item{rhs}{right-hand-side vector for hypothesis, with as many entries as
    rows in the hypothesis matrix; can be omitted, in which case it defaults
    to a vector of zeroes.}
  \item{test}{character specifying whether to compute the finite sample
    F statistic (with approximate F distribution) or the large sample
    Chi-squared statistic (with asymptotic Chi-squared distribution).}
  \item{vcov.}{a function for estimating the covariance matrix of the regression
    coefficients, e.g., \code{\link{hccm}}, or an estimated covariance matrix
    for \code{model}. See also \code{white.adjust}.}        
  \item{white.adjust}{logical or character. Convenience interface to \code{hccm}
    (instead of using the argument \code{vcov}). Can be set either to a character
    specifying the \code{type} argument of \code{\link{hccm}} or \code{TRUE},
    in which case \code{"hc3"} is used implicitly. For backwards compatibility.}
  \item{verbose}{If \code{TRUE}, the hypothesis matrix and right-hand-side
    vector are printed to standard output; if \code{FALSE} (the default),
    the hypothesis is only printed in symbolic form.}
  \item{...}{aruments to pass down.}
}

\details{
  Computes either a finite sample F statistic or asymptotic Chi-squared
  statistic for carrying out a Wald-test-based comparison between a model
  and a linearly restricted model. The default method will work with any
  model object for which the coefficient vector can be retrieved by 
  \code{coef} and the coefficient-covariance matrix by \code{vcov} (otherwise
  the argument \code{vcov.} has to be set explicitely). For computing the
  F statistic (but not the Chi-squared statistic) a \code{df.residual} 
  method needs to be available. If a \code{formula} method exists, it is
  used for pretty printing.

  The method for \code{"lm"} objects calls the default method, but it
  changes the default test to \code{"F"}, supports the convenience argument
  \code{white.adjust} (for backwards compatibility), and enhances the output
  by residual sums of squares. For \code{"glm"} objects just the default
  method is called (bypassing the \code{"lm"} method).

  The function \code{lht} also dispatches to \code{linear.hypothesis}.

  The hypothesis matrix can be supplied as a numeric matrix (or vector), 
  the rows of which specify linear combinations of the model coefficients, 
  which are tested equal to the corresponding entries in the righ-hand-side 
  vector, which defaults to a vector of zeroes. 

  Alternatively, the
  hypothesis can be specified symbolically as a character vector with one
  or more elements, each of which gives either a linear combination of
  coefficients, or a linear equation in the coefficients (i.e., with both
  a left and right side separated by an equals sign). Components of a
  linear expression or linear equation can consist of numeric constants, or
  numeric constants multiplying coefficient names (in which case the number
  precedes the coefficient, and may be separated from it by spaces or an
  asterisk); constants of 1 or -1 may be omitted.  Spaces are always optional.
  Components are separated by positive or negative signs. See the examples below.}
}

\value{
  An object of class \code{"anova"} which contains the residual degrees of freedom
  in the model, the difference in degrees of freedom, Wald statistic
  (either \code{"F"} or \code{"Chisq"}) and corresponding p value.
}



\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}

\author{Achim Zeileis and John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[stats]{anova}}, \code{\link{Anova}}, \code{\link[lmtest]{waldtest}},
  \code{\link{hccm}}, \code{\link[sandwich]{vcovHC}}, \code{\link[sandwich]{vcovHAC}}, 
  \code{\link[stats]{coef}}, \code{\link[stats]{vcov}}

}

\examples{
mod.davis <- lm(weight~repwt, data=Davis)

## the following are equivalent:
linear.hypothesis(mod.davis, diag(2), c(0,1))
linear.hypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"))
linear.hypothesis(mod.davis, c("(Intercept)", "repwt"), c(0,1))
linear.hypothesis(mod.davis, c("(Intercept)", "repwt = 1"))

## use asymptotic Chi-squared statistic
linear.hypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"), test = "Chisq")


## the following are equivalent:
  ## use HC3 standard errors via white.adjust option
linear.hypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"), 
    white.adjust = TRUE)
  ## covariance matrix *function*
linear.hypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"), vcov = hccm)
  ## covariance matrix *estimate*
linear.hypothesis(mod.davis, c("(Intercept) = 0", "repwt = 1"), 
    vcov = hccm(mod.davis, type = "hc3"))

mod.duncan <- lm(prestige ~ income + education, data=Duncan)

## the following are all equivalent:
linear.hypothesis(mod.duncan, "1*income - 1*education = 0")
linear.hypothesis(mod.duncan, "income = education")
linear.hypothesis(mod.duncan, "income - education")
linear.hypothesis(mod.duncan, "1income - 1education = 0")
linear.hypothesis(mod.duncan, "0 = 1*income - 1*education")
linear.hypothesis(mod.duncan, "income-education=0")
linear.hypothesis(mod.duncan, "1*income - 1*education + 1 = 1")
linear.hypothesis(mod.duncan, "2income = 2*education")

mod.duncan.2 <- lm(prestige ~ type*(income + education), data=Duncan)
coefs <- names(coef(mod.duncan.2))

## test against the null model (i.e., only the intercept is not set to 0)
linear.hypothesis(mod.duncan.2, coefs[-1]) 

## test all interaction coefficients equal to 0
linear.hypothesis(mod.duncan.2, coefs[grep(":", coefs)], verbose=TRUE) 
}

\keyword{htest}
\keyword{models}
\keyword{regression}

