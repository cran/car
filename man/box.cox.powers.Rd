\name{box.cox.powers}
\alias{box.cox.powers}
\alias{print.box.cox.powers}
\alias{summary.box.cox.powers}

\title{Multivariate Unconditional Box-Cox Transformations}
\description{
  Estimates multivariate unconditional power transformations
  to multinormality
  by the method of maximum likelihood. The univariate case is
  obtained when only one variable is specified.
}
\usage{
box.cox.powers(X, start=NULL, hypotheses=NULL, ...)

\method{print}{box.cox.powers}(x, digits=4, ...)

\method{summary}{box.cox.powers}(object, digits=4, ...)
}

\arguments{
  \item{X}{a numeric matrix of variables (or a vector for one variable)
    to be transformed.}
  \item{start}{start values for the power transformation parameters;
    if \code{NULL} (the default), univariate Box-Cox transformations will
    be computed and used as the start values.}
  \item{hypotheses}{if non-\code{NULL}, a list of hypotheses to be tested;
    each hypothesis should be a vector of values giving the power for each 
    column of \code{X}. Note that the hypotheses that all powers are 1 and
    that all powers are 0 (log) are always tested.}
  \item{...}{optional arguments to be passed to the \code{optim} function.}
  \item{digits}{number of places to round result.}
  \item{x, object}{\code{box.cox.powers} object.}
}
\details{
  Note that this is \emph{unconditional} Box-Cox. That is, there is
  no regression model, and there are no predictors. The object is to
  make the distribution of the variable(s) as (multi)normal as possible.
  For Box-Cox regression, see the \code{boxcox} function in the 
  \code{MASS} package.
  
  The function estimates the Box-Cox powers, 
  \eqn{x_{j}^{\prime }=(x_{j}^{\lambda _{j}}-1)/\lambda _{j}}{x' = (x^p - 1)/p} 
  for \eqn{\lambda _{j} \neq 0}{p != 0} and \eqn{x_{j}^{\prime }=\log x_{j}}{x' = log(x)} 
  for \eqn{\lambda _{j}=0}{p = 0}. Subsequently using ordinary power
  transformations (i.e., \eqn{x^p} for \eqn{p \neq 0}{p != 0}) 
  achieves the same result.
}
\value{
  returns an object of class \code{box.cox.powers}, which may be printed
  or summarized. the \code{print} and \code{summary} methods are now identical; I've 
  retained the latter for backwards compatibility.
}

\references{
  Box, G. E. P. and Cox, D. R. (1964)
  An analysis of transformations.
  \emph{JRSS B} \bold{26} 211--246.
  
  Cook, R. D. and Weisberg, S. (1999)
  \emph{Applied Regression, Including Computing and Graphics.} Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[MASS:boxcox]{boxcox}}, \code{\link{box.cox}}, \code{\link{box.cox.var}}, 
    \code{\link{box.cox.axis}}}

\examples{
attach(Prestige)
box.cox.powers(cbind(income, education))
## Box-Cox Transformations to Multinormality  
## 
##           Est.Power Std.Err. Wald(Power=0) Wald(Power=1) 
## income       0.2617   0.1014         2.580        -7.280 
## education    0.4242   0.4033         1.052        -1.428 
## 
## L.R. test, all powers = 0:  7.694   df = 2   p = 0.0213 
## L.R. test, all powers = 1:  48.8727   df = 2   p = 0  
plot(income, education)
plot(box.cox(income, .26), box.cox(education, .42))

box.cox.powers(income)
## Box-Cox Transformation to Normality 
## 
##  Est.Power Std.Err. Wald(Power=0) Wald(Power=1)
##     0.1793   0.1108         1.618        -7.406
## 
## L.R. test, power = 0:  2.7103   df = 1   p = 0.0997
## L.R. test, power = 1:  47.261   df = 1   p = 0 
qq.plot(income)
qq.plot(income^.18)
}

\keyword{multivariate}
\keyword{models}
