\name{influence.mixed.models}
\alias{influence.mixed.models}
\alias{influence.merMod}
\alias{influence.lme}
\alias{cooks.distance.influence.merMod}
\alias{cooks.distance.influence.lme}
\alias{dfbeta.influence.merMod}
\alias{dfbeta.influence.lme}
\alias{dfbetas.influence.merMod}
\alias{dfbetas.influence.lme}

\title{
Influence Diagnostics for Mixed-Effects Models
}

\description{
These functions compute deletion influence diagnostics for linear mixed-effects models fit by \code{\link[lme4]{lmer}} in the \pkg{lme4} package and
\code{\link[nlme]{lme}} in the \pkg{nlme} package and for generalized linear mixed-effects models fit by  \code{\link[lme4]{glmer}} in the
\pkg{lme4} package. The main functions are methods for the \code{\link{influence}} generic function. Other functions are provided for
computing \code{\link{dfbeta}}, \code{\link{dfbetas}}, \code{\link{cooks.distance}}, and influence on variance-covariance components based
on the objects computed by \code{influence.merMod} and \code{influence.lme}.
}


\usage{
\method{influence}{merMod}(model, groups, data, maxfun=1000, ncores=1, ...)
\method{influence}{lme}(model, groups, data, ncores=1, ...)

\method{cooks.distance}{influence.merMod}(model, ...)
\method{cooks.distance}{influence.lme}(model, ...)
\method{dfbeta}{influence.merMod}(model, which = c("fixed", "var.cov"), ...)
\method{dfbeta}{influence.lme}(model, which = c("fixed", "var.cov"), ...)
\method{dfbetas}{influence.merMod}(model, ...)
\method{dfbetas}{influence.lme}(model, ...)
}

\arguments{
  \item{model}{in the case of \code{influence.merMod} or \code{influence.lme}, a model of class \code{"merMod"} or \code{"lme"};
  in the case of \code{cooks.distance}, \code{dfbeta}, or \code{dfbetas}, an object returned by
  \code{influence.merMod} or \code{influence.lme}.}
  \item{groups}{a character vector containing the name of a grouping factor or names of grouping factors; if more than one name
  is supplied, then groups are defined by all combinations of levels of the grouping factors that appear in
  the data. If omitted, then each individual row of the data matrix is treated as a "group" to be deleted in turn.}
  \item{data}{an optional data frame with the data to which \code{model} was fit; \code{influence.merMod} can usually reconstruct
  the data, and \code{influence.lme} can access the data unless \code{keep.data=FALSE} was specified in the call to \code{lme},
  so it's usually unnecessary to supply the \code{data} argument.}
  \item{maxfun}{The maximum number of function evaluations (for \code{influence.merMod})
  to perform after deleting each group; the defaults are large enough so that the iterations will typically continue to convergence.
  Setting to \code{maxfun=20} for an \code{lmer} model or \code{100} for a \code{glmer} model will typically produce a faster reasonable approximation. 
  An even smaller value can be used if interest is only in influence on the fixed effects.}
  \item{ncores}{number of cores for parallel computation of diagnostics; if \code{1}
  (the default), the computation isn't parallelized; if \code{Inf}, all of the available \emph{physical} cores 
  (not necessarily \emph{logical} cores --- see \code{\link[parallel]{detectCores}}) on the computer will be used.}
  \item{which}{if \code{"fixed.effects"} (the default), return influence on the fixed effects; if \code{"var.cov"}, return influence on the variance-covariance components.}
  \item{\dots}{ignored.}
}

\details{
\code{influence.merMod} and \code{influence.lme} start with the estimated variance-covariance components from \code{model} and then refit
the model omitting each group in turn, not necessarily iterating to completion. For example, \code{maxfun=20} takes up to 20 function evaluations
step away from the ML or REML solution for the full data, which usually provides decent approximations to the fully iterated estimates.

The other functions are methods for the \code{\link{dfbeta}}, \code{\link{dfbetas}}, and \code{\link{cooks.distance}} generics, to be applied to the
\code{"influence.merMod"} or \code{"influence.lme"} object produced by the \code{influence} function; the \code{dfbeta} methods can also return
influence on the variance-covariance components.
}

\value{\code{influence.merMod} and \code{influence.lme} return objects of class \code{"influence.merMod"} and \code{"influence.lme"} respectively,
each of which contains the following elements:
\describe{
    \item{\code{"fixed.effects"}}{the estimated fixed effects for the model.}
    \item{\code{"fixed.effects[-groups]"}}{a matrix with columns corresponding to the fixed-effects coefficients and rows corresponding to groups, giving the
    estimated fixed effects with each group deleted in turn; \emph{groups} is formed from the name(s) of the grouping factor(s).}
    \item{\code{"var.cov.comps"}}{the estimated variance-covariance parameters for the model.}
    \item{\code{"var.cov.comps[-groups]"}}{a matrix with the estimated covariance parameters (in columns) with each group deleted in turn.}
    \item{\code{"vcov"}}{The estimated covariance matrix of the fixed-effects coefficients.}
    \item{\code{"vcov[-groups]"}}{a list each of whose elements is the estimated covariance matrix of the fixed-effects coefficients with one group deleted.}
    \item{\code{"groups"}}{a character vector giving the names of the grouping factors.}
    \item{\code{"deleted"}}{the possibly composite grouping factor, each of whose elements is deleted in turn.}
    \item{\code{"converged"}}{for \code{influence.merMod}, a logical vector indicating whether the computation converged for each group.}
    \item{\code{"function.evals"}}{for \code{influence.merMod}, a vector of the number of function evaluations performed for each group.}
}
For plotting \code{"influence.merMod"} and \code{"influence.lme"} objects, see \code{\link{infIndexPlot}}.
}

\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{
J. Fox \email{jfox@mcmaster.ca}
}

\seealso{
\code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}}, \code{\link[nlme]{lme}}, \code{\link{infIndexPlot}}.
}

\examples{
if (require("lme4")){
    print(fm1 <- lmer(Reaction ~ Days + (Days | Subject),
        sleepstudy)) # from ?lmer
    infIndexPlot(influence(fm1, "Subject"))
    }
if (require("lme4")){
    gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
        data = cbpp, family = binomial)
    infIndexPlot(influence(gm1, "herd", maxfun=100))
    gm1.11 <- update(gm1, subset = herd != 11) # check deleting herd 11
    compareCoefs(gm1, gm1.11)
    }
}

\keyword{models}
