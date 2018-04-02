\name{subsets}
\alias{subsets}
\alias{subsets.regsubsets}
\title{Plot Output from regsubsets Function in leaps package}
\description{
  The \code{\link[leaps]{regsubsets}} function in the \pkg{leaps} package finds
  optimal subsets of predictors based on some criterion statistic. This function plots a measure of fit against subset size.
}
\usage{
subsets(object, ...)

\method{subsets}{regsubsets}(object, 
    names=abbreviate(object$xnames, minlength = abbrev), 
    abbrev=1, min.size=1, max.size=length(names), 
    legend="interactive", 
    statistic=c("bic", "cp", "adjr2", "rsq", "rss"), 
    las=par('las'), cex.subsets=1, ...)
}

\arguments{
  \item{object}{a \code{regsubsets} object produced by the \code{regsubsets} function
    in the \pkg{leaps} package.}
  \item{names}{a vector of (short) names for the predictors, excluding the
    regression intercept, if one is present; if missing, these are
    derived from the predictor names in \code{object}.}
  \item{abbrev}{minimum number of characters to use in abbreviating predictor names.}
  \item{min.size}{minimum size subset to plot; default is \code{1}.}
  \item{max.size}{maximum size subset to plot; default is number of predictors.}
  \item{legend}{If not \code{FALSE}, in which case the legend is suppressed,
    the coordinates at which to place a legend of the abbreviated predictor names on the plot, 
    in a form recognized by
    the \code{\link{legend}} function. If \code{"interactive"}, the legend is placed on the
    plot interactively with the mouse. By expanding the left or right plot margin, you can place
    the legend in the margin, if you wish (see \code{\link{par}}).}
  \item{statistic}{statistic to plot for each predictor subset; one of: 
    \code{"bic"}, Bayes Information Criterion; 
    \code{"cp"}, Mallows's \eqn{C_{p}}{Cp};
    \code{"adjr2"}, \eqn{R^{2}}{R^2} adjusted for degrees of freedom;
    \code{"rsq"}, unadjusted \eqn{R^{2}}{R^2};
    \code{"rss"}, residual sum of squares.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{cex.subsets}{can be used to change the relative size of the characters used to
    plot the regression subsets; default is \code{1}.}
  \item{\dots}{arguments to be passed down to 
    \code{subsets.regsubsets} and \code{plot}.}
}

\value{
  \code{NULL} if the \code{legend} is \code{TRUE}; otherwise a data frame with the legend.
}
 
\author{John Fox}

\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\seealso{\code{\link[leaps]{regsubsets}}}

\examples{
if (require(leaps)){
    subsets(regsubsets(undercount ~ ., data=Ericksen),
            legend=c(3.5, -37))
}
}

\keyword{hplot}
\keyword{regression}
