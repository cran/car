\name{linear.hypothesis}
\alias{linear.hypothesis}
\alias{lht}
\alias{linear.hypothesis.lm}
\alias{linear.hypothesis.glm}
\alias{print.F.test}
\alias{print.chisq.test}

\title{Test Linear Hypothesis}
\description{
Test a linear hypothesis for a linear or generalized linear model.
}
\usage{
linear.hypothesis(model, ...)

lht(...)

linear.hypothesis.lm(model, hypothesis.matrix, rhs=0, 
  summary.model=summary(model, corr = FALSE), 
  white.adjust=FALSE, error.SS, error.df, ...)

linear.hypothesis.glm(model, hypothesis.matrix, rhs=0, 
  summary.model=summary(model, corr = FALSE), ...)
  
\method{print}{chisq.test}(x, ...)

\method{print}{F.test}(x, ...)
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
  \item{white.adjust}{if \code{TRUE} use heteroscedasticity-corrected covariance matrix.}
  \item{error.SS}{error sum of squares for the hypothesis; if not specified, will be
    taken from \code{model}.}
  \item{error.df}{error degrees of freedom for the hypothesis; if not specified,
    will be taken from \code{model}.}
  \item{x}{\code{chisq.test} or \code{F.test} object.}
  \item{...}{aruments to pass down.}
}
\details{
  Computes an F-test for the hypothesis in a linear model, or a Wald test in a
  generalized linear model.
}
\value{
  Returns an \code{F.test} or \code{chisq.test} object, with components:
  \item{SSH}{sum of squares for hypothesis (for a linear model).}
  \item{SSE}{error sum of squares (for a linear model).}
  \item{f}{F-statistic for the hypothesis (for a linear model.)}
  \item{Df}{degrees of freedom for F or chisquare.}
  \item{p}{p-value for the hypothesis.}
  \item{ChiSquare}{chisquare statistic for the hypothesis (for a generalized linear model).}
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[base]{anova}}, \code{\link{Anova}}, \code{\link{hccm}}}

\examples{
data(Davis)
mod<-lm(weight~repwt, data=Davis)
linear.hypothesis(mod, diag(2), c(0,1))
## F-Test 
## SS = 245.9738     SSE = 12828.03     F = 1.735312  Df = 2 and 181     p = 0.179266 
}

\keyword{htest}
\keyword{models}
\keyword{regression}
