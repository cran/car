\name{S}
\alias{S}
\alias{Confint}
\alias{S.lm}
\alias{S.glm}
\alias{S.default}
\alias{S.multinom}
\alias{S.polr}
\alias{S.lme}
\alias{S.lmerMod}
\alias{S.glmerMod}
\alias{print.S.lm}
\alias{print.S.glm}
\alias{print.S.multinom}
\alias{print.S.polr}
\alias{print.S.lme}
\alias{print.S.lmerMod}
\alias{print.S.glmerMod}
\alias{Confint.lm}
\alias{Confint.glm}
\alias{Confint.multinom}
\alias{Confint.polr}
\alias{Confint.lme}
\alias{Confint.lmerMod}
\alias{Confint.glmerMod}
\alias{Confint.default}


\title{
Modified Functions for Summarizing Linear, Generalized Linear, and Some Other Models
}

\description{
\pkg{car} package replacements for the \code{\link{summary}} (\code{S}) and \code{\link{confint}} (\code{Confint}) functions for \code{\link{lm}}, \code{\link{glm}}, \code{\link[nnet]{multinom}}, and \code{\link[MASS]{polr}} objects, with additional arguments but the same defaults as the original functions. The \code{Confint} method for \code{"polr"} objects profiles the likelihood to get confidence intervals for the regression parameters but uses Wald intervals for the thresholds.
Default methods that call the standard R \code{\link{summary}} and \code{\link{confint}} functions are provided for the \code{S} and \code{Confint} generics, so the \pkg{car} functions should be safe to use in general. The default method for \code{Confint} also assumes that there is an appropriate \code{\link{coef}} method. For briefer model summaries, see \code{\link{brief}}.
}

\usage{

S(object, brief, ...)

\method{S}{lm}(object, brief=FALSE,
    correlation = FALSE, symbolic.cor = FALSE,
    vcov. = vcov(object, complete=FALSE), header = TRUE,
    resid.summary = FALSE, adj.r2 = FALSE,
    ...)

\method{S}{glm}(object, brief=FALSE,
    exponentiate, dispersion, correlation = FALSE, symbolic.cor = FALSE,
    vcov. = vcov(object, complete=FALSE), header = TRUE,
    resid.summary = FALSE, ...)

\method{S}{multinom}(object, brief=FALSE, exponentiate=FALSE, ...)

\method{S}{polr}(object, brief=FALSE, exponentiate=FALSE, ...)

\method{S}{lme}(object, brief=FALSE, correlation=FALSE, ...)

\method{S}{lmerMod}(object, brief=FALSE, KR=FALSE, correlation=FALSE, ...)

\method{S}{glmerMod}(object, brief=FALSE, correlation=FALSE, exponentiate, ...)

\method{S}{default}(object, brief, ...)

\method{print}{S.lm}(x, digits = max(3, getOption("digits") - 3),
    symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"), ...)

\method{print}{S.glm}(x, digits = max(3L, getOption("digits") - 3L),
    symbolic.cor = x$symbolic.cor, signif.stars = getOption("show.signif.stars"), ...)

\method{print}{S.multinom}(x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)

\method{print}{S.polr}(x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)

\method{print}{S.lme}(x, digits=max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)

\method{print}{S.lmerMod}(x, digits=max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)

\method{print}{S.glmerMod}(x, digits=max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)

Confint(object, ...)

\method{Confint}{lm}(object, estimate=TRUE,
    parm, level=0.95, vcov.=vcov(object, complete=FALSE), ...)

\method{Confint}{glm}(object, estimate=TRUE, exponentiate=FALSE,
    vcov., dispersion, type=c("LR", "Wald"), ...)

\method{Confint}{polr}(object, estimate=TRUE, exponentiate=FALSE,
    thresholds=!exponentiate, ...)

\method{Confint}{multinom}(object, estimate=TRUE, exponentiate=FALSE, ...)

\method{Confint}{lme}(object, estimate=TRUE, level=0.95, ...)

\method{Confint}{lmerMod}(object, estimate=TRUE, level=0.95, ...)

\method{Confint}{glmerMod}(object, estimate=TRUE, level=0.95,
    exponentiate=FALSE, ...)

\method{Confint}{default}(object, estimate=TRUE, level=0.95, vcov., ...)
}

\arguments{
  \item{object}{
a model object, e.g., of class \code{"lm"} as produced by a call to \code{\link{lm}}.
}
  \item{exponentiate}{for a \code{"glm"} or \code{"glmerMod"} model using the \code{log} or \code{logit} link, or a
    \code{"polr"} or \code{"multinom"} model, show exponentiated coefficient estimates and confidence bounds.}
  \item{correlation, symbolic.cor}{see \code{\link{summary.lm}}}
  \item{x, digits, signif.stars}{see \code{\link{summary.lm}}}
  \item{dispersion}{see \code{\link{summary.glm}}}
  \item{vcov.}{either a matrix giving the estimated covariance matrix of the estimates, or a function that
   when called with \code{object} as an argument returns an estimated covariance matrix of the estimates.
   The default of \code{vcov. = vcov} uses the usual estimated covariance matrix.  Other choices include
   the functions documented at \code{\link{hccm}}.  See example below for using a bootstrap to estimate  the
   covariance matrix.  For the \code{glm} methods of \code{Confint} and \code{S}, if the \code{vcov.} or
   \code{dispersion} argument is specified, then Wald-based confidence limits are computed;
   otherwise the reported confidence limits are computed by profiling the likelihood.
   NOTE:  The \code{dispersion} and \code{vcov.} arguments may not \emph{both} be specified.
}
  \item{header}{
if \code{TRUE}, print the header for the summary output, default is \code{TRUE}
  }
  \item{resid.summary}{
if \code{TRUE}, print the five-number summary of the residuals in the summary, defaults to \code{FALSE}
  }
  \item{adj.r2}{
if \code{TRUE}, print the adjusted r-squared in the summary, default is \code{FALSE}
  }
  \item{brief}{
if \code{TRUE}, set \code{header}, \code{resid.summary} and \code{adj.r.squared} to \code{FALSE}, and suppress exponeniated coefficients for GLMs with log or logit link.
  }
  \item{KR}{if \code{TRUE} (default is \code{FALSE}), compute Kenward-Roger standard errors and
  Satterthwaite degrees of freedom for t-tests. \emph{Warning:} This computation can be very time-consuming.}
  \item{parm, level}{see \code{\link{confint}}}
  \item{estimate}{show the estimated coefficients in the confidence-interval table; default is \code{TRUE}.}
  \item{thresholds}{show confidence intervals for the estimated thresholds in the \code{"polr" model.}}
  \item{type}{if \code{"LR"} (the default) compute confidence intervals based on the LR statistics by profiling the likelihood; if \code{"Wald"} base confidence intervals on the Wald statistic using the coefficient standard error and the normal distribution.}
  \item{\dots}{
additional arguments to be passed down, for consistency with \code{summary} and \code{confint} methods
}
}
\details{
All these functions mimic functions in the \pkg{stats} and other standard R packages for summarizing aspects of linear, generalized linear, and some other statistical models.

The \code{S} and \code{Confint} functions add support for the \code{vcov.} argument for linear models, which allows specifying a covariance matrix for the regression coefficients other than the usual covariance matrix returned by the function \code{\link{vcov}}.  This argument may be either the name of a function, so that the call to \code{vcov.(object)} returns a covariance matrix, or else \code{vcov.} is set equal to a covariance matrix.  For example, setting \code{vcov.=hccm} uses 'proposal 3' described by Long and Ervin (2000) for a sandwich coefficient-variance estimator that may be robust against nonconstant variance (see \code{\link{hccm}}).  Setting \code{vcov. = hccm(object, type = "hc2")} would use the matrix returned by the \code{hccm} function using proposal 2.  For use with a bootstrap, see the examples below. The overall F-test in the \code{S.lm} output uses the supplied covariance matrix in a call to the \code{\link{linearHypothesis}} function.

The supplied \code{print} method for \code{S.lm} (and for other \code{S} methods) has additional arguments to customize the standard \code{summary.lm} output.  Standard output is obtained by setting \code{resid.summary=TRUE, adj.r2=TRUE}.

Using a heterscedasticy-corrected covariance matrix computed using \code{\link{hccm}} with GLMs other than Gaussian is not justified; see the article by Freedman (2006).

The \code{Summary.glm} method for models fit with the log or logit link by default prints a table of exponentiated coefficients and their confidence limits; \code{Summary.multinom} and \code{Summary.polr} optionally print tables of exponentiated coefficients.
}

\value{
The \code{S.lm} and \code{S.glm} functions return a list with all the elements shown at \code{\link{summary.lm}} and \code{\link{summary.glm}}. The \code{S.multinom} and \code{S.polr} functions return a list with all the elements shown at \code{\link{summary.multinom}} and \code{summary.polr} plus potentially a table of exponentiated coefficients and confidence bounds.

The \code{Confint.lm} function returns either the output from \code{\link{confint.lm}} if
\code{vcov. = vcov} or Wald-type confidence intervals using the supplied covariance matrix for any other choice of \code{vcov.}.

Finally, \code{Confint} applied to any object that does not inherit from \code{"lm"}, \code{"multinom"}, or \code{"polr"} simply calls \code{confint}, along with, by default, using \code{\link{coef}} to add a column of estimates to the confidence limits.
}

\references{

  Freedman, David A. (2006). On the so-called Huber sandwich estimator and robust standard errors.
  \emph{The American Statistician}, \bold{60}, 299-302.

  Long, J. S. and Ervin, L. H. (2000)
  Using heteroscedasity consistent standard errors in the linear regression model.
  \emph{The American Statistician} \bold{54}, 217--224.

  White, H. (1980)
  A heteroskedastic consistent covariance matrix estimator and a direct test of heteroskedasticity.
  \emph{Econometrica} \bold{48}, 817--838.
}
\author{
Sanford Weisberg \email{sandy@umn.edu}
}


\seealso{\code{\link{brief}}, \code{\link{summary}}, \code{\link{confint}}, \code{\link{coef}},
\code{\link{summary.lm}}, \code{\link{confint}}, \code{\link{vcov.lm}}, \code{\link{hccm}},
\code{\link{Boot}}, \code{\link{linearHypothesis}}}

\examples{
mod.prestige <- lm(prestige ~ education + income + type, Prestige)
S(mod.prestige, vcov.=hccm)
S(mod.prestige, brief=TRUE)
Confint(mod.prestige, vcov.=hccm)

# A logit model
mod.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
S(mod.mroz)

\dontrun{
# Using the bootstrap for standard errors
b1 <- Boot(mod.prestige)
S(mod.prestige, vcov.= vcov(b1))
Confint(b1) # run with the boot object to get corrected confidence intervals
}
}


\keyword{models}
\keyword{regression}
