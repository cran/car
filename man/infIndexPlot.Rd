\name{infIndexPlot}
\alias{infIndexPlot}
\alias{influenceIndexPlot}
\alias{infIndexPlot.lm}
\alias{infIndexPlot.lmerMod}
\alias{infIndexPlot.influence.merMod}
\alias{infIndexPlot.influence.lme}
\title{Influence Index Plot}
\description{
Provides index plots of influence and related diagnostics for a regression model.
}
\usage{
infIndexPlot(model, ...)

influenceIndexPlot(model, ...)

\method{infIndexPlot}{lm}(model, vars=c("Cook", "Studentized", "Bonf", "hat"),
    id=TRUE, grid=TRUE, main="Diagnostic Plots", ...)

\method{infIndexPlot}{influence.merMod}(model,
    vars = c("dfbeta", "dfbetas", "var.cov.comps",
    "cookd"), id = TRUE, grid = TRUE, main = "Diagnostic Plots", ...)
\method{infIndexPlot}{influence.lme}(model,
    vars = c("dfbeta", "dfbetas", "var.cov.comps",
    "cookd"), id = TRUE, grid = TRUE, main = "Diagnostic Plots", ...)
}
\arguments{
  \item{model}{A regression object of class \code{lm}, \code{glm}, or \code{lmerMod}, or an influence
  object for a \code{lmer}, \code{glmer}, or \code{lme} object (see
  \code{\link{influence.mixed.models}}). The \code{"lmerMod"} method calls the \code{"lm"} method and can take the same arguments.}
  \item{vars}{
All the quantities listed in this argument are plotted.  Use \code{"Cook"}
for Cook's distances, \code{"Studentized"} for Studentized
residuals, \code{"Bonf"} for Bonferroni p-values for an outlier test, and
and \code{"hat"} for hat-values (or leverages) for a linear or generalized
linear model, or \code{"dfbeta"}, \code{"dfbetas"}, \code{"var.cov.comps"}, and
\code{"cookd"} for an influence object derived from a mixed model.   Capitalization is optional.
All but \code{"dfbeta"} and \code{"dfbetas"} may be abbreviated by the first one or more letters.
}
  \item{main}{main title for graph}
  \item{id}{a list of named values controlling point labelling. The default, \code{TRUE}, is
    equivalent to \code{id=list(method="y", n=2, cex=1, col=carPalette()[1], location="lr")};
    \code{FALSE} suppresses point labelling. See \code{\link{showLabels}} for details.}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the graph.}
  \item{\dots}{Arguments passed to \code{plot}}
}

\value{
Used for its side effect of producing a graph.  Produces index plots
of diagnostic quantities.
}
\references{
  Cook, R. D. and Weisberg, S. (1999)
  \emph{Applied Regression, Including Computing and Graphics.} Wiley.

  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.
  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.

  Weisberg, S. (2014)
  \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}

\author{Sanford Weisberg \email{sandy@umn.edu} and John Fox}

\seealso{ \code{\link{cooks.distance}}, \code{\link{rstudent}},
  \code{\link{outlierTest}}, \code{\link{hatvalues}}, \code{\link{influence.mixed.models}}. }
\examples{
influenceIndexPlot(lm(prestige ~ income + education + type, Duncan))

\dontrun{ # a little slow
  if (require(lme4)){
      print(fm1 <- lmer(Reaction ~ Days + (Days | Subject),
          sleepstudy)) # from ?lmer
      infIndexPlot(influence(fm1, "Subject"))
      infIndexPlot(influence(fm1))
      }
      
  if (require(lme4)){
      gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
          data = cbpp, family = binomial) # from ?glmer
      infIndexPlot(influence(gm1, "herd", maxfun=100))
      infIndexPlot(influence(gm1, maxfun=100))
      gm1.11 <- update(gm1, subset = herd != 11) # check deleting herd 11
      compareCoefs(gm1, gm1.11)
      }
    }
}

\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
