\name{skewPower}
\alias{skewPower}
\alias{contour.skewpowerTransform}

\title{Skew Power Transformations}
\description{
Transform the elements of a vector using, the skew power transformations.
}
\usage{
skewPower(U, lambda, jacobian.adjusted=FALSE, gamma)

\method{contour}{skewpowerTransform}(x, ksds = 4, levels = c(0.5, 0.95, 0.99, 0.999), 
                 main = "Skew Power Log-likelihood", ...)
}
\arguments{
  \item{U}{A vector, matrix or data.frame of values to be transformed}
  \item{lambda}{The one-dimensional transformation parameter, usually in 
  the range from \eqn{-2} to \eqn{2}, or if \code{U} is a matrix or data frame, a vector of length
  \code{ncol(U)} of transformation parameters}
  \item{jacobian.adjusted}{If \code{TRUE}, the transformation is normalized to have
  Jacobian equal to one.  The default is \code{FALSE}.}
  \item{gamma}{The skew power family has two parameters, the power lambda, and the start gamma, which can be any positive value.}
  \item{x}{
  An object created by a call to \code{\link{powerTransform}} using the \code{family="skewPower"} option.
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
  \item{\dots}{
Additional arguments passed to the generic contour function.
}
}

\details{
  The \code{skewPower} family is a modification of the Box-Cox power family, see \code{\link{bcPower}} that allows the untransformed response to include negative values.  Ignoring the Jacobian, correction, the Box-Cox family is defined to be \eqn{(U^{\lambda}-1)/\lambda}{(U^(lambda)-1)/lambda} 
  for \eqn{\lambda \neq 0}{lambda not equal to zero}, and 
  \eqn{\log(U)}{log(U)} if \eqn{\lambda =0}{lambda = 0}.  The \code{skewPower family} is the Box-Cox transformation of \eqn{z = .5 * (U + (U^2 + \gamma^2)^{1/2})}. Zero or negative values for \eqn{U}{\code{U}} are permitted if \eqn{\gamma}{gamma} is positive. 

If \code{jacobian.adjusted} is \code{TRUE}, then the transformed values are scaled to have Jacobian equal to 1, as suggested by Box and Cox (1964).  This simplifies computing the maximum likelihood-like estimates of the transformation parameters.  
  
Missing values are permitted, and return \code{NA} where ever \code{U} is equal to \code{NA}.

The contour method allows drawing a contour plot for the two-dimensional log-likelihood for the skew power parameters with a univariate response.
}

\value{
  Returns a vector or matrix of transformed values.  The contour method is used for the side-effect of drawing a contour plot.
}

\references{
Fox, J. and Weisberg, S. (2011) 
\emph{An R Companion to Applied Regression}, Second Edition, Sage.

Hawkins, D. and Weisberg, S. (2015)
Combining the Box-Cox Power and Genralized Log Transformations to Accomodate Negative Responses,
submitted for publication.
  
Weisberg, S. (2014)  \emph{Applied Linear Regression}, Fourth Edition, Wiley
Wiley, Chapter 7.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{\code{\link{powerTransform}}, \code{\link{bcPower}}, \code{\link{boxCox}}}
\examples{
# Univariate Example from Hawkins and Weisberg (2015)
m1 <- lm(I1L1 ~ pool, LoBD)
p1 <- powerTransform(m1, family="skewPower")
summary(p1)
# summary prints estimate, se and conf.ints for both parameters
# helper functions
coef(p1)
vcov(p1)
p1$value # maximum value of log-likelihood surface

# tests are for lambda, maximizing over gamma (profile log-likelihoods)
testTransform(p1, lambda=0.5)

# Contour plot of the log-likelihood
contour(p1, main="", levels=c(.5, .95, .99))

# the boxCox function can provide profile log-likelihoods for each of the two parameters:
boxCox(m1, family="skewPower", param="lambda", lambda=seq(0.25, 1.1, length=100))
boxCox(m1, family="skewPower", param="gamma", gamma=seq(3, 80, length=100))

#  Create new variate corresponding to the transformation, then update the regression
LoBD$z <- skewPower(LoBD$I1L1, lambda=p1$lambda, gamma=p1$gamma, jacobian.adjusted=FALSE)
m1.t <- update(m1, z ~ . -1)
#  The example in Hawkins and Weisberg(2015) requires computing the inverse of the skew power
# transformation and certain summaries of the updated regression
mu0 <- mean(coef(m1.t)[1:4])
sd <- sigmaHat(m1.t)
LoB <- mu0 + 1.645*sd
LoD <- LoB + 1.645*sd
skewPowerInverse <- function(z, lambda, gamma){
  q <- if(abs(lambda) <= 1.e-9) {2 * log(z)} else {2*(z*lambda + 1) ^(1/lambda)}
  (q^2 - gamma^2)/(2 * q)
}
y.untransformed <- round(skewPowerInverse(c(LoB, LoD), p1$lambda, p1$gamma), 2)

# Multivariate Response
p3 <- powerTransform(update(m1, as.matrix(cbind(LoBD$I1L2, LoBD$I1L1)) ~ .), family="skewPower")
summary(p3)
# Construct test for the same transformations for each response.  Get new data file by stacking
d2 <- data.frame(assay = c(LoBD$I1L2, LoBD$I1L1), pool=factor(c(LoBD$pool, LoBD$pool)),
                 set=factor(rep(c("A", "B"), c(84, 84))))
# the following allows for different within set means, but common transformation parameters
p5 <- powerTransform(lm(assay ~ pool*set, d2), family="skewPower") # allow grps to differ by set.
c(df=(df <- length(coef(p3)) - length(coef(p5))),
  test = (test <- -.5 * (p5$value - p3$value)),
  pvalue = 1 - pchisq(test, df))

}

\keyword{regression}
