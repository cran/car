\name{box.cox.var}
\alias{box.cox.var}

\title{Constructed Variable for Box-Cox Transformation}
\description{
  Computes a constructed variable for the Box-Cox transformation of the
  response variable in a linear model.
}
\usage{
box.cox.var(y)
}

\arguments{
  \item{y}{response variable.}
}
\details{
  The constructed variable is defined as \eqn{y[\log (y/\widetilde{y})-1]}{y[log(y/y') -1]}, where 
  \eqn{\widetilde{y}}{y'} is the geometric mean of \code{y}. 
  
  The constructed variable is meant to be
  added to the right-hand-side of the linear model. The t-test for the
  coefficient of the constructed variable is an approximate score test for whether a
  transformation is required. 
  
  If \eqn{b} is the coefficient of the constructed variable,
  then an estimate of the normalizing power transformation based on the score statistic
  is \eqn{1-b}{1 - b}. An added-variable plot for the constructed
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
  
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage. 
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[MASS]{boxcox}}, \code{\link{box.cox}}, \code{\link{box.cox.powers}}, 
  \code{\link{box.cox.axis}}, \code{\link{av.plots}}}

\examples{
mod <- lm(interlocks + 1 ~ assets, data=Ornstein)
mod.aux <- update(mod, . ~ . + box.cox.var(interlocks + 1))
summary(mod.aux)
## Call:
## lm(formula = interlocks + 1 ~ assets + box.cox.var(interlocks + 
##     1), data = Ornstein)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -23.189  -6.701   0.541   6.773  12.051 
## 
## Coefficients:
##                               Estimate Std. Error t value Pr(>|t|)
## (Intercept)                  1.461e+01  5.426e-01  26.920   <2e-16
## assets                      -7.142e-05  5.119e-05  -1.395    0.164
## box.cox.var(interlocks + 1)  7.427e-01  4.136e-02  17.956   <2e-16
## 
## Residual standard error: 7.247 on 245 degrees of freedom
## Multiple R-Squared: 0.7986,     Adjusted R-squared: 0.797 
## F-statistic: 485.7 on 2 and 245 degrees of freedom,     p-value:     0 
av.plots(mod.aux, "box.cox.var(interlocks + 1)")
}
\keyword{manip}
\keyword{regression}
