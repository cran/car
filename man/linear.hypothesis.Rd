\name{linear.hypothesis}
\alias{linear.hypothesis}
\alias{lht}
\alias{linear.hypothesis.lm}
\alias{linear.hypothesis.glm}

\title{Test Linear Hypothesis}
\description{
Generic function for testing a linear hypothesis, and methods
for fitted linear or generalized linear models.
}
\usage{
linear.hypothesis(model, ...)

lht(...)

\method{linear.hypothesis}{lm}(model, hypothesis.matrix, rhs=0, 
  summary.model=summary(model, corr = FALSE), 
  test=c("F", "Chisq"), vcov=NULL,
  white.adjust=FALSE, error.SS, error.df, ...)

\method{linear.hypothesis}{glm}(model, hypothesis.matrix, rhs=0, 
  summary.model=summary(model, corr = FALSE),
  test=c("Chisq", "F"), vcov=NULL, error.df, ...)
  
}
\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.}
  \item{hypothesis.matrix}{matrix (or vector) giving linear combinations
    of coefficients by rows.}
  \item{rhs}{right-hand-side vector for hypothesis, with as many entries as
    rows in \code{hypothesis.matrix}.}.
  \item{summary.model}{a \code{summary} object for the model; usually specified
    only when \code{linear.hypothesis} is called from another function that has
    already computed the summary.}
  \item{test}{character specifying wether to compute the finite sample
    F statistic (with approximate F distribution) or the large sample
    Chi-squared statistic (with asymptotic Chi-squared distribution).}
  \item{vcov}{a function for estimating the covariance matrix of the regression
    coefficients, e.g., \code{\link{hccm}} or an estimated covariance matrix
    for \code{model}. See also \code{white.adjust}.}        
  \item{white.adjust}{logical or character. Convenience interface to \code{hccm}
    (instead of using the argument \code{vcov}). Can be set either to a character
    specifying the \code{type} argument of \code{\link{hccm}} or \code{TRUE},
    then \code{"hc3"} is used implicitly.}
  \item{error.SS}{error sum of squares for the hypothesis; if not specified, will be
    taken from \code{summary.model}.}
  \item{error.df}{error degrees of freedom for the hypothesis; if not specified,
    will be taken from \code{summary.model}.}
  \item{...}{aruments to pass down.}
}
\details{
  Computes either a finite sample F statistic (default for \code{"lm"} objects)
  or asymptotic Chi-squared statistic (default for \code{"glm"} objects) for
  carrying out a Wald-test-based comparison between a model and a linearly 
  restricted model.
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
\author{John Fox \email{jfox@mcmaster.ca} and Achim Zeleis}

\seealso{\code{\link[base]{anova}}, \code{\link{Anova}}, \code{\link{hccm}},
         \code{\link[lmtest]{waldtest}}}

\examples{
data(Davis)
mod<-lm(weight~repwt, data=Davis)
linear.hypothesis(mod, diag(2), c(0,1))

## use asymptotic Chi-squared statistic
linear.hypothesis(mod, diag(2), c(0,1), test = "Chisq")

## use HC3 standard errors via
## white.adjust option
linear.hypothesis(mod, diag(2), c(0,1), white.adjust = TRUE)
## covariance matrix *function*
linear.hypothesis(mod, diag(2), c(0,1), vcov = hccm)
## covariance matrix *estimate*
linear.hypothesis(mod, diag(2), c(0,1), vcov = hccm(mod, type = "hc3"))
}

\keyword{htest}
\keyword{models}
\keyword{regression}
