\name{influence.mixed.models}
\alias{influence.mixed.models}
\alias{influence.lme}
\alias{cooks.distance.influence.lme}
\alias{dfbeta.influence.lme}
\alias{dfbetas.influence.lme}

\title{
Influence Diagnostics for Mixed-Effects Models
}

\description{
These functions compute deletion influence diagnostics for linear mixed-effects models fit by \code{\link[nlme]{lme}} in the \pkg{nlme} package.
The main function is a method for the \code{\link{influence}} generic function. Other functions are provided for
computing \code{\link{dfbeta}}, \code{\link{dfbetas}}, \code{\link{cooks.distance}}, and influence on variance-covariance components based
on the object computed by \code{influence.lme}.
}


\usage{
\method{influence}{lme}(model, groups, data, ncores=1, ...)

\method{cooks.distance}{influence.lme}(model, ...)
\method{dfbeta}{influence.lme}(model, which = c("fixed", "var.cov"), ...)
\method{dfbetas}{influence.lme}(model, ...)
}

\arguments{
  \item{model}{in the case \code{influence}, a model of class  \code{"lme"};
  in the case of \code{cooks.distance}, \code{dfbeta}, or \code{dfbetas}, an object returned by \code{influence.lme}.}
  \item{groups}{a character vector containing the name of a grouping factor or names of grouping factors; if more than one name
  is supplied, then groups are defined by all combinations of levels of the grouping factors that appear in
  the data. If omitted, then each individual row of the data matrix is treated as a "group" to be deleted in turn.}
  \item{data}{an optional data frame with the data to which \code{model} was fit; 
  \code{influence.lme} can access the data unless \code{keep.data=FALSE} was specified in the call to \code{lme},
  so it's usually unnecessary to supply the \code{data} argument.}
  \item{ncores}{number of cores for parallel computation of diagnostics; if \code{1}
  (the default), the computation isn't parallelized; if \code{Inf}, all of the available \emph{physical} cores 
  (not necessarily \emph{logical} cores --- see \code{\link{detectCores}}) on the computer will be used.}
  \item{which}{if \code{"fixed.effects"} (the default), return influence on the fixed effects; if \code{"var.cov"}, return influence on the variance-covariance components.}
  \item{\dots}{ignored.}
}

\details{
\code{influence.lme} starts with the estimated variance-covariance components from \code{model} and then refits
the model omitting each group in turn.

The other functions are methods for the \code{\link{dfbeta}}, \code{\link{dfbetas}}, and \code{\link{cooks.distance}} generics, to be applied to the
\code{"influence.lme"} object produced by the \code{influence} function; the \code{dfbeta} methods can also return
influence on the variance-covariance components.
}

\value{
\code{influence.lme} returns an object of class %\code{"influence.merMod"} and 
\code{"influence.lme"},% respectively,
which contains the following elements:
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
}
For plotting \code{"influence.lme"} objects, see \code{\link{infIndexPlot}}.
}

\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{
J. Fox \email{jfox@mcmaster.ca}
}

\seealso{
\code{\link[nlme]{lme}}, \code{\link{infIndexPlot}}.
}

\examples{

if (require(nlme)){
    print(fm1 <- lme(distance ~ age, data = Orthodont))
    infIndexPlot(influence(fm1, "Subject"))
    infIndexPlot(influence(fm1))
    }
}

\keyword{models}
