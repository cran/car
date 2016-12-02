\name{skewPower}
\alias{skewPower}
\alias{contour.skewpowerTransform}
\alias{contour.skewpowerTransformlmer}

\title{Skew Power Transformations}
\description{
Transform the elements of a vector using, the skew power transformations.
}
\usage{
skewPower(U, lambda, jacobian.adjusted=FALSE, gamma)

\method{contour}{skewpowerTransform}(x, ksds = 4, levels = c(0.5, 0.95, 0.99, 0.999), 
                 main = "Skew Power Log-likelihood", ...)
                 
\method{contour}{skewpowerTransformlmer}(x, ksds = 4, levels = c(0.5, 0.95, 0.99, 0.999), 
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

Hawkins, D. and Weisberg, S. (2016)
Combining the Box-Cox Power and Generalized Log Transformations to Accomodate Negative Responses
In Linear and Mixed-Effects Linear Models,
submitted for publication.
  
Weisberg, S. (2014)  \emph{Applied Linear Regression}, Fourth Edition, Wiley
Wiley, Chapter 7.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{\code{\link{powerTransform}}, \code{\link{bcPower}}, \code{\link{boxCox}}}
\examples{
# Univariate Example, from Hawkins and Weisberg (2016)
m1 <- lm(I1L1 ~ pool, LoBD)
p1 <- powerTransform(m1, family="skewPower")
summary(p1)
# summary prints estimate, se and conf.ints for both parameters
# helper functions
c(p1$lambda, p1$gamma, LogLik=p1$llik)
vcov(p1) # Estimated covarinace from inverse Hessian
# tests are for lambda, maximizing over gamma (profile log-likelihoods
testTransform(p1, lambda=0.5)
# Contour plot of the log-likelihood
contour(p1, main="", levels=c(.5, .95, .99))
# the boxCox function can provide profile log-likelihoods for each of the two parameters:
boxCox(m1, family="skewPower", param="lambda", lambda=seq(0.25, 1.1, length=100))
boxCox(m1, family="skewPower", param="gamma", gamma=seq(3, 80, length=100))

# Fit with fixed gamma, in this case fixed at the estimate from p1
p2 <- powerTransform(m1, family="skewPower", gamma=p1$gamma)
# summary gives different tests because gamma is fixed rather than maximized over
summary(p2)

# Multivariate Response
p3 <- powerTransform(update(m1, as.matrix(cbind(LoBD$I1L2, LoBD$I1L1)) ~ .), 
     family="skewPower")
summary(p3)

# gamma fixed
p4 <- powerTransform(update(m1, as.matrix(cbind(LoBD$I1L2, LoBD$I1L1)) ~ .), 
     family="skewPower", gamma=p3$gamma)
summary(p4)

# mixed models fit with lmer - requires lmer and nloptr packages
\dontrun{
data <- reshape(LoBD[1:20,], varying=names(LoBD)[-1], direction="long", v.names="y")
names(data) <- c("pool", "assay", "y", "id")
data$assay <- factor(data$assay)
require(lme4)
m2 <- lmer(y ~ pool + (1|assay), data)
f1 <- powerTransform(m2, family="skewPower")
summary(f1)
f2 <- powerTransform(m2, family="skewPower", gamma= 10)
summary(f2)
}
}

\keyword{regression}
