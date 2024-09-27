\name{boxCoxVariable}
\alias{boxCoxVariable}

\title{Constructed Variable for Box-Cox Transformation}
\description{
  Computes a constructed variable for the Box-Cox transformation of the
  response variable in a linear model.
}
\usage{
boxCoxVariable(y)
}

\arguments{
  \item{y}{response variable.}
}

\details{
  The constructed variable is defined as \eqn{y[\log(y/\widetilde{y}) - 1]}{y[log(y/y') -1]}, where 
  \eqn{\widetilde{y}}{y'} is the geometric mean of \code{y}. 
  
  The constructed variable is meant to be
  added to the right-hand-side of the linear model. The t-test for the
  coefficient of the constructed variable is an approximate score test for whether a
  transformation is required. 
  
  If \eqn{b} is the coefficient of the constructed variable,
  then an estimate of the normalizing power transformation based on the score statistic
  is \eqn{1 - b}{1 - b}. An added-variable plot for the constructed
  variable shows leverage and influence on the decision to transform \code{y}. 
}
\value{
  a numeric vector of the same length as \code{y}.
}
\references{
  Atkinson, A. C. (1985)
  \emph{Plots, Transformations, and Regression}. Oxford.
  
  Box, G. E. P. and Cox, D. R. (1964)
  An analysis of transformations.
  \emph{JRSS B} \bold{26} 211--246.
  
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[MASS]{boxcox}}, \code{\link{powerTransform}}, \code{\link{bcPower}}} % , \code{\link{avPlots}}

\examples{
mod <- lm(interlocks + 1 ~ assets, data=Ornstein)
mod.aux <- update(mod, . ~ . + boxCoxVariable(interlocks + 1))
summary(mod.aux)
# avPlots(mod.aux, "boxCoxVariable(interlocks + 1)")
}
\keyword{manip}
\keyword{regression}
