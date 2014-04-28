\name{boxCox}
\alias{boxCox}
\alias{boxCox.lm}
\alias{boxCox.default}
\alias{boxCox.formula}

\title{
Box-Cox Transformations for Linear Models 
}
\description{
Computes and optionally plots profile log-likelihoods for the parameter of the 
Box-Cox power transformation.   This is a slight generalization of the
\code{boxcox} function in the \pkg{MASS} package that allows for families of transformations
other than the Box-Cox power family.
}
\usage{
boxCox(object, ...)

\method{boxCox}{default}(object, lambda = seq(-2, 2, 1/10), plotit = TRUE,
       interp = (plotit && (m < 100)), eps = 1/50, 
       xlab = expression(lambda),
       ylab = "log-Likelihood", family="bcPower", grid=TRUE, ...)

\method{boxCox}{formula}(object, lambda = seq(-2, 2, 1/10), plotit = TRUE,
       interp = (plotit && (m < 100)), eps = 1/50, 
       xlab = expression(lambda),
       ylab = "log-Likelihood", family="bcPower", ...)

\method{boxCox}{lm}(object, lambda = seq(-2, 2, 1/10), plotit = TRUE,
       interp = (plotit && (m < 100)), eps = 1/50, 
       xlab = expression(lambda),
       ylab = "log-Likelihood", family="bcPower", ...)
}

\arguments{
  \item{object}{
a formula or fitted model object. Currently only \code{lm} and \code{aov} objects are handled.
}
  \item{lambda}{
vector of values of lambda, with default (-2, 2) in steps of 0.1, where the profile
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
defaults to \code{"lambda"}.
} 
  \item{ylab}{
defaults to \code{"log-Likelihood"}. 
}
  \item{family}{
Defaults to \code{"bcPower"} for the Box-Cox power family of
transformations.  If set to \code{"yjPower"} the Yeo-Johnson family, which
permits negative responses, is used.
}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph.
}
  \item{\dots}{
additional parameters to be used in the model fitting.
}
}
\details{
This routine is an elaboration of the \code{\link{boxcox}} function in the
\pkg{MASS} package.  All arguments except for \code{family} and \code{grid} are 
identical, and if the arguments 
\code{family = "bcPower", grid=FALSE} is set it gives an identical graph.  If 
\code{family = "yjPower"} then the Yeo-Johnson power transformations, which
allow nonpositive responses, will be used.
}

\value{
A list of the lambda vector and the computed profile log-likelihood vector, 
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

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.

Yeo, I. and Johnson, R. (2000)  A new family of
power transformations to improve normality or symmetry.  
\emph{Biometrika}, 87, 954-959.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{
\code{\link{boxcox}}, \code{\link{yjPower}}, \code{\link{bcPower}},
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
