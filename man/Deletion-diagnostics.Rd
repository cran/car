\name{Deletion Diagnostics}

\alias{rstudent}
\alias{rstudent.lm}
\alias{rstudent.glm}

\alias{hatvalues}
\alias{hatvalues.lm}

\alias{cookd}
\alias{cookd.lm}
\alias{cookd.glm}

\alias{dfbeta}
\alias{dfbeta.lm}

\alias{dfbetas}
\alias{dfbetas.lm}

\alias{influence}
\alias{influence.lm}
\alias{influence.glm}

\title{Deletion Diagnostics for Linear and Generalized Linear Models}
\description{
  These functions calculate a variety of leave-one-out deletion diagnostics for
  linear and generalized linear models, including studentized residuals (for
  outlier detection), hatvalues (for detecting high-leverage observations), and
  Cook's distances, dfbeta, and dfbetas (for detecting influential observations).
}

\usage{

rstudent(model, ...)

rstudent.lm(model, infl=influence(model))

rstudent.glm(model, infl=influence(model))

hatvalues(model, ...)

hatvalues.lm(model, infl=influence(model))

cookd(model, ...)

cookd.lm(model, infl=influence(model), sumry=summary(model))

cookd.glm(model, infl=influence(model), sumry=summary(model))

dfbeta(model, ...)

dfbeta.lm(model, infl=influence(model))

dfbetas(model, ...)

dfbetas.lm(model, infl=influence(model), sumry=summary(model))

influence(model, ...)

influence.lm(model)

influence.glm(model)

}

\arguments{
  \item{model}{\code{lm} or \code{glm} model object.}
  \item{infl}{optionally, an influence-object precomputed for the \code{model} by \code{influence}.}
  \item{sumry}{optionally, a summary-object precomputed for the \code{model} by \code{summary}.}
  \item{...}{arguments to be passed down from generic functions to method functions.}
}

\details{
  Basic quantities are computed by \code{influence.lm} or \code{influence.glm}, which are slightly
  modified versions of \code{lm.influence} from the base package. Values for generalized linear
  models are approximations, as described in Williams (1987) (except that Cook's distances are
  scaled as \emph{F} rather than as chi-square values).
  
  Normally, the generic versions of these functions are the ones to be used directly. For 
  \code{hatvalues}, \code{dfbeta}, and \code{dfbetas}, the method for linear models
  also works for generalized linear models. 
  
  The following diagnostics are provided:
  
  \describe{
    \item{\code{rstudent}} {studentized residuals.}
    \item{\code{hatvalues}} {observation leverages.}
    \item{\code{cookd}} {Cook's distance influence measure for observations.}
    \item{\code{dfbeta}} {change in each coefficient upon deleting observations.}
    \item{\code{dfbetas}} {standardized change in each coefficient for deleting observations.}
    }
}

\value{
  \code{rstudent}, \code{hatvalues}, and \code{cookd} return vectors with one entry for
  each observation; \code{dfbeta} and \code{dfbetas} return matrices with rows for 
  observations and columns for coefficients.
  
  \code{influence} returns a list with entries:
  
  \item{hat}{hat-values.}
  \item{sigma}{leave-one-out estimates of linear-model standard error or 
    generalized-linear-model scale.}
  \item{coefficients}{dfbeta values.}
  \item{wt.res}{weighted residuals (for a linear model).}
  \item{dev.res}{deviance residuals (for a generalized linear model).}
  \item{pear.res}{Pearson residuals (for a generalized linear model).}
  
}

\references{
  Belsley, D. A. and Kuh, E. and Welsch, R. E. (1980)
  \emph{Regression Diagnostics.} Wiley.
  
  Cook, R. D. and Weisberg, S. (1984)
  \emph{Residuals and Influence in Regression.) Wiley.
  
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
  
  Williams, D. A. (1987)
  Generalized linear model diagnostics using the deviance and single
  case deletions. \emph{Applied Statistics} \bold{36}, 181--191.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[base:influence.measures]{influence.measures}}}

\examples{
data(Duncan)
attach(Duncan)
mod <- lm(prestige ~ income + education)
qq.plot(rstudent(mod), distribution="t", df=41)
plot(hatvalues(mod))
plot(cookd(mod))
plot(dfbeta(mod)[,2])
plot(dfbetas(mod)[,2])
}

\keyword{regression}
\keyword{models}
