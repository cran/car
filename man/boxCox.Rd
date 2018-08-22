\name{boxCox}
\alias{boxCox}
\alias{boxCox2d}
\alias{boxCox.lm}
\alias{boxCox.default}
\alias{boxCox.formula}
\alias{boxCox.bcnPowerTransform}

\title{
Graph the profile log-likelihood for Box-Cox transformations in 1D, or in 2D with the bcnPower family.
}
\description{
  Computes and optionally plots profile log-likelihoods for the parameter of the
  Box-Cox power family, the Yeo-Johnson power family, or for either of the parameters in a bcnPower family.   This is a slight generalization of the
  \code{boxcox} function in the \pkg{MASS} package that allows for families of 
  transformations other than the Box-Cox power family.  the \code{boxCox2d} function
  produces a  contour
  plot of the two-dimensional likelihood profile for the bcnPower family.
}
\usage{
boxCox(object, ...)

\method{boxCox}{default}(object,
        lambda = seq(-2, 2, 1/10), plotit = TRUE,
        interp = plotit, eps = 1/50,
        xlab=NULL, ylab=NULL,
        family="bcPower",
        param=c("lambda", "gamma"), gamma=NULL,
        grid=TRUE, ...)

\method{boxCox}{formula}(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, family = "bcPower",
    param = c("lambda", "gamma"), gamma = NULL, grid = TRUE,
    ...)

\method{boxCox}{lm}(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, ...)

boxCox2d(x, ksds = 4, levels = c(0.5, 0.95, 0.99, 0.999),
                 main = "bcnPower Log-likelihood", grid=TRUE, ...)
}

\arguments{
  \item{object}{
    a formula or fitted model object of class \code{lm} or \code{aov}.
  }
  \item{lambda}{
    vector of values of \eqn{\lambda}{lambda}, with default (-2, 2) in steps of 0.1, where the profile
    log-likelihood will be evaluated.
  }
  \item{plotit}{
    logical which controls whether the result should be plotted; default \code{TRUE}.
  }
  \item{interp}{
    logical which controls whether spline interpolation is used. Default to
    \code{TRUE} if plotting with lambda of length less than 100.
  }
  \item{eps}{
    Tolerance for lambda = 0; defaults to 0.02.
  }
  \item{xlab}{
    defaults to \code{"lambda"} or \code{"gamma"}.
  }
  \item{ylab}{
    defaults to \code{"log-Likelihood"} or for bcnPower family to the appropriate label.
  }

  \item{family}{
    Defaults to \code{"bcPower"} for the Box-Cox power family of
    transformations.  If set to \code{"yjPower"} the Yeo-Johnson family, which
    permits negative responses, is used.  If set to \code{bcnPower} the function gives the profile
    log-likelihood for the parameter selected via \code{param}.
  }
  \item{param}{Relevant only to \code{family="bcnPower"}, produces a profile log-likelihood for the parameter selected, maximizing over the remaining parameter.}
 \item{gamma}{
    For use when the \code{family="bcnPower", param="gamma"}.  If this is a vector of positive values, then the profile log-likelihood for the location (or start) parameter in the  \code{bcnPower} family is evaluated at these values of gamma.  If gamma is \code{NULL}, then evaulation is done at 100 equally spaced points between \code{min(.01, gmax - 3*sd)} and \code{gmax + 3*sd}, where \code{gmax} is the maximimum likelihood estimate of gamma, and \code{sd} is the sd of the response.  See \code{\link{bcnPower}} for the definition of \code{gamma}.
  }
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
              graph.
  }
  \item{\dots}{
  additional arguments passed to the \code{lm} method with \code{boxCox.formula} or passed to \code{contour} with \code{boxCox2d}.
  }
    \item{x}{
  An object created by a call to \code{\link{powerTransform}} using \code{family="bcnPower"}.
}
  \item{ksds}{
Contour plotting of the log-likelihood surface will cover plus of minus \code{ksds} standard deviations on each axis.
}
  \item{levels}{
Contours will be drawn at the values of levels.  For example, \code{levels=c(.5, .99)} would display two contours, at the 50\% level and at the 99\% level.
}
  \item{main}{
Title for the contour plot
}
}

\details{
The \code{boxCox} function is an elaboration of the \code{\link{boxcox}} function in the
  \pkg{MASS} package.  The first 7 arguments are the same as in \code{boxcox}, and if the argument \code{family="bcPower"} is used, the result is essentially identical to the function in \pkg{MASS}.  Two additional families are the \code{yjPower} and \code{bcnPower} families that allow a few values of the response to be non-positive.
The bcnPower family has two parameters: a power \eqn{\lambda}{lambda} and a start or location parameter \eqn{\gamma}{gamma}, and the \code{boxCox} function can be used to obtain a profile log-likelihood for either parameter with \eqn{\lambda}{lambda} as the default. Alternatively, the \code{boxCox2d} function can be used to get a contour plot of the profile log-likelihood.
}

\value{
  Both functions ae designed for their side effects of drawing a graph.  The \code{boxCox} function returns a list of the lambda (or possibly, gamma) vector and the computed profile log-likelihood vector,
  invisibly if the result is plotted.  If \code{plotit=TRUE} plots log-likelihood vs
  lambda and indicates a 95\% confidence interval about the maximum observed value of
  lambda. If \code{interp=TRUE}, spline interpolation is used to give a smoother plot.

}
\references{

  Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations.
  \emph{Journal
        of the Royal Statisistical Society, Series B}. 26 211-46.

  Cook, R. D. and Weisberg, S. (1999)  \emph{Applied Regression Including
                                             Computing
                                             and Graphics}.  Wiley.

  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.

  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.

Hawkins, D. and Weisberg, S. (2017)
Combining the Box-Cox Power and Generalized Log Transformations to Accomodate Nonpositive  Responses In Linear and Mixed-Effects Linear Models \emph{South African Statistics Journal}, 51, 317-328.

  Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.

  Yeo, I. and Johnson, R. (2000)  A new family of
  power transformations to improve normality or symmetry.
  \emph{Biometrika}, 87, 954-959.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{
  \code{\link{boxcox}}, \code{\link{yjPower}}, \code{\link{bcPower}}, \code{\link{bcnPower}},
  \code{\link{powerTransform}}, \code{\link{contour}}
}
\examples{
  with(trees, boxCox(Volume ~ log(Height) + log(Girth), data = trees,
         lambda = seq(-0.25, 0.25, length = 10)))

  data("quine", package = "MASS")
  with(quine, boxCox(Days ~ Eth*Sex*Age*Lrn, 
         lambda = seq(-0.05, 0.45, len = 20), family="yjPower"))
}
\keyword{ regression}
