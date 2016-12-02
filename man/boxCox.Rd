\name{boxCox}
\alias{boxCox}
\alias{boxCox.lm}
\alias{boxCox.default}
\alias{boxCox.formula}
\alias{boxCox.skewpowerTransform}
\alias{boxCox.lmerMod}

\title{
  Box-Cox Transformations for Linear Models 
}
\description{
  Computes and optionally plots profile log-likelihoods for the parameter of the 
  Box-Cox power family, the Yeo-Johnson power family, or for either of the parameters in a skew power family.   This is a slight generalization of the
  \code{boxcox} function in the \pkg{MASS} package that allows for families of transformations
  other than the Box-Cox power family.
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
  
\method{boxCox}{formula}(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, ...)
  
\method{boxCox}{lm}(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, ...)

\method{boxCox}{lmerMod}(object,
          lambda = seq(-2, 2, 1/10), plotit = TRUE, interp = plotit, eps = 1/50,
          xlab=NULL, ylab=NULL, family="bcPower", 
          param=c("lambda", "gamma"), gamma=NULL, grid=TRUE, ...)
}

\arguments{
  \item{object}{
    a formula or fitted model object of class \code{lm} or \code{aov} or \code{lmerMod}.
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
    defaults to \code{"log-Likelihood"} or for skewPower family to the appropriate label. 
  }
  
  \item{family}{
    Defaults to \code{"bcPower"} for the Box-Cox power family of
    transformations.  If set to \code{"yjPower"} the Yeo-Johnson family, which
    permits negative responses, is used.  If set to \code{skewPower} the function gives the profile
    log-likelihood for the parameter selected via \code{param}.
  }  
  \item{param}{Relevant only to \code{family="skewPower"}, produces a profile log-likelihood for the parameter selected, maximizing over the remaining parameter.}
 \item{gamma}{
    For use when the \code{family="skewPower", param="gamma"}.  If this is a vector of positive values, then the profile log-likelihood for the location (or start) parameter in the skew power family is evaluated at these values of gamma.  If gamma is \code{NULL}, then evaulation is done at 100 equally spaced points between \code{min(.01, gmax - 3*se)} and \code{gmax + 3*se}, where \code{gmax} is the maximimul likelihood estimate of gamma, and \code{se} is its estimated standard error.  See \code{\link{skewPower}} for the definition of \code{gamma}.
  } 
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
              graph.
  }
  \item{\dots}{
    passes arguments from methods to the default, or from the default to \code{plot}.
  }
}
\details{
  This routine is an elaboration of the \code{\link{boxcox}} function in the
  \pkg{MASS} package.  The first 7 arguments are the same as in \code{boxcox}, and if the argument \code{family="bcPower"} is used, the result is essentially identical to the function in \pkg{MASS}.  Two additional families are the \code{yjPower} and \code{skewPower} families that allow a few values of the response to be non-positive. 
The skew power family has two parameters: a power \eqn{\lambda}{lambda} and a start or location parameter \eqn{\gamma}{gamma}, and this function can be used to obtain a profile log-likelihood for either parameter with \eqn{\lambda}{lambda} as the default.
}

\value{
  A list of the lambda (or possibly, gamma) vector and the computed profile log-likelihood vector, 
  invisibly if the result is plotted.  If \code{plotit=TRUE} plots log-likelihood vs 
  lambda and indicates a 95% confidence interval about the maximum observed value of 
  lambda. If \code{interp=TRUE}, spline interpolation is used to give a smoother plot. 
  
}
\references{
  
  Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations.
  \emph{Journal
        of the Royal Statisistical Society, Series B}. 26 211-46.
  
  Cook, R. D. and Weisberg, S. (1999)  \emph{Applied Regression Including
                                             Computing
                                             and Graphics}.  Wiley.
  
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Hawkins, D. and Weisberg, S. (2015)
  Combining the Box-Cox Power and Genralized Log Transformations to Accomodate Negative Responses,
  submitted for publication.
  
  Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.
  
  Yeo, I. and Johnson, R. (2000)  A new family of
  power transformations to improve normality or symmetry.  
  \emph{Biometrika}, 87, 954-959.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{
  \code{\link{boxcox}}, \code{\link{yjPower}}, \code{\link{bcPower}}, \code{\link{skewPower}},
  \code{\link{powerTransform}}
}
\examples{
  boxCox(Volume ~ log(Height) + log(Girth), data = trees,
         lambda = seq(-0.25, 0.25, length = 10))
  
  data("quine", package = "MASS")
  boxCox(Days ~ Eth*Sex*Age*Lrn, data = quine,
         lambda = seq(-0.05, 0.45, len = 20), family="yjPower")
}
\keyword{ regression}
