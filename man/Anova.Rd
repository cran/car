\name{Anova}
\alias{Anova}
\alias{Anova.lm}
\alias{Anova.aov}
\alias{Anova.II.lm}
\alias{Anova.III.lm}
\alias{Anova.glm}
\alias{Anova.II.F.glm}
\alias{Anova.II.LR.glm}
\alias{Anova.II.Wald.glm}
\alias{Anova.III.F.glm}
\alias{Anova.III.LR.glm}
\alias{Anova.III.Wald.glm}
\title{Anova Tables for Linear and Generalized Linear Models}
\description{
  Calculates type-II or type-III analysis-of-variance tables for
  model objects produced by \code{lm} and \code{glm}. For linear
  models, F-tests are calculated; for generalized linear models, 
  likelihood-ratio chisquare, Wald chisquare, or F-tests are calculated.
}
\usage{
Anova(mod, ...)

Anova.lm(mod, error, type=c("II", "III"), ...)

Anova.aov(mod, ...)

Anova.glm(mod, type=c("II", "III"), test.statistic=c("LR", "Wald", "F"), 
    error, error.estimate=c("pearson", "dispersion", "deviance"))
}

\arguments{
  \item{mod}{\code{lm} or \code{glm} model object.}
  \item{error}{for a linear model, an \code{lm} model object from which the
    error sum of squares and degrees of freedom are to be calculated. For 
    F-tests for a generalized linear model, a \code{glm} object from which the
    dispersion is to be estimated. If not specified, \code{mod} is used.}
  \item{type}{type of test, \code{"II"} or \code{"III"}.}
  \item{test.statistic}{for a generalized linear model, whether to calculate 
    \code{"LR"} (likelihood-ratio), \code{"Wald"}, or \code{"F"} tests.}
  \item{error.estimate}{for F-tests for a generalized linear model, base the
    dispersion estimate on the Pearson residuals (\code{pearson}, the default); use the
    dispersion estimate in the model object (\code{dispersion}), which, e.g., is
    fixed to 1 for binomial and Poisson models; or base the dispersion estimate on
    the residual deviance (\code{deviance}).}
  \item{\dots}{arguments to be passed to \code{linear.hypothesis}; only use
    \code{white.adjust} for a linear model.}
}
\details{
  So-called type-II tests are calculated according to the principle of marginality,
  testing each term after all others, except ignoring the term's higher-order relatives;
  so-called type-III tests violate marginality, testing 
  each term in the model after all of the others. Be
  very careful in formulating the model for type-III tests, or the hypotheses tested
  will not make sense. As implemented here, type-II Wald tests 
  for generalized linear models are actually
  \emph{differences} of Wald statistics.
  
  For all but type-II likelihood-ratio and \emph{F} tests for generalized linear models, 
  \code{Anova} finds the test statistics without refitting the model.
  
  The standard R \code{anova} function calculates sequential (type-I) tests.
  These rarely test meaningful hypotheses.
}
\value{
  An object of class \code{anova}, usually printed.
}
\references{ 
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\section{Warning}{Be careful of type-III tests.}

\seealso{ \code{\link{linear.hypothesis}}, \code{\link[base]{anova}}}

\examples{
data(Moore)
mod<-lm(conformity~fcategory*partner.status, data=Moore, 
  contrasts=list(fcategory=contr.sum, partner.status=contr.sum))
Anova(mod)
## Anova Table (Type II tests)
##
## Response: conformity
##                         Sum Sq Df F value   Pr(>F)
## fcategory                 11.61  2  0.2770 0.759564
## partner.status           212.21  1 10.1207 0.002874
## fcategory:partner.status 175.49  2  4.1846 0.022572
## Residuals                817.76 39                 
Anova(mod, type="III")
## Anova Table (Type III tests)
##
## Response: conformity
##                          Sum Sq Df  F value    Pr(>F)
## (Intercept)              5752.8  1 274.3592 < 2.2e-16
## fcategory                  36.0  2   0.8589  0.431492
## partner.status            239.6  1  11.4250  0.001657
## fcategory:partner.status  175.5  2   4.1846  0.022572
## Residuals                 817.8 39                   
}
\keyword{htest}
\keyword{models}
\keyword{regression}
