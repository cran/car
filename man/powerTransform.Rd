\name{powerTransform}
\alias{powerTransform}
\alias{powerTransform.default}
\alias{powerTransform.lm}
\alias{powerTransform.formula}
\alias{powerTransform.lmerMod}

\title{Finding Univariate or Multivariate Power Transformations}
\description{
  \code{powerTransform} uses the maximum likelihood-like approach of Box and Cox (1964) to select a transformatiion of a univariate or multivariate response for normality, linearity and/or constant variance.  Available families of transformations are the default Box-Cox power family and two additioal families that are modifications of the Box-Cox family that allow for (a few) negative responses. The summary method automatically computes two or three likelihood ratio type tests concerning the transformation powers.
}
\usage{
powerTransform(object, ...)

\S3method{powerTransform}{default}(object, family="bcPower", ...)

\S3method{powerTransform}{lm}(object, family="bcPower", ...)

\S3method{powerTransform}{formula}(object, data, subset, weights, na.action,
    family="bcPower", ...)

\S3method{powerTransform}{lmerMod}(object, family="bcPower", ...)
}

\arguments{
\item{object}{This can either be an object of class \code{lm} or \code{lmerMod}, a formula, or a
   matrix or vector; see below.
}
\item{family}{The quoted name of a family of transformations.  The available options are \code{"bcPower"}
   for the default for the Box-Cox power family; \code{"bcnPower"} for a two-parameter modification of
   the Box-Cox family that allows negative responses (Hawkins and Weisberg (2017)), and the
   \code{"yjPower"} family (Yeo and Johnson(2000)), another
   modifiation of the Box-Cox family that allows a few negative values.  All three families are documented
   at \code{\link{bcPower}}.
}
  \item{data}{A data frame or environment, as in \sQuote{\link{lm}}.}
  \item{subset}{Case indices to be used, as in \sQuote{\link{lm}}.}
  \item{weights}{Weights as in \sQuote{\link{lm}}.}
  \item{na.action}{Missing value action, as in \sQuote{lm}.}
  \item{...}{Additional arguments that used in the interative algorithm; defaults are generally adequate.  For use with the \code{bcnPower} family, a convergence criterion can be set with \code{conv=.0001} the default, and a minimum positive value of the location parameter can be set, with default \code{gamma.min=.1}.}
}

\details{This function implements the Box and Cox (1964) method of selecting a power transformation of a variable toward normality, and its generalization by Velilla (1993) to a multivariate response.  Cook and Weisberg (1999) and Weisberg (2014) suggest the usefulness of transforming a set of predictors \code{z1, z2, z3} for multivariate normality.  It also includes two additional families that allow for negative values.

If the \code{object} argument is of class \sQuote{lm} or \sQuote{lmerMod}, the Box-Cox procedure is applied to the conditional distribution of the response given the predictors.  For \sQuote{lm} objects, the respose may be multivariate, and each column will have its own transformation.  With \sQuote{lmerMod} the response must be univariate.

The \code{object} argument may also be a formula.  For example, \code{z ~ x1 + x2 + x3} will estimate a transformation for the response \code{z} from a family after fitting a linear model with the given formula.  \code{cbind(y1, y2, y3) ~ 1} specifies transformations
to multivariate normality with no predictors.  A vector value for \code{object}, for example
\code{powerTransform(ais$LBM)}, is equivalent to\code{powerTransform(LBM ~ 1, ais)}.  Similarly, \code{powerTransform(cbind(ais$LBM, ais$SSF))}, where the first argument is a matrix rather than a formula is equivalent to specification of a mulitvariate linear model \code{powerTransform(cbind(LBM, SSF) ~ 1, ais)}.

Three families of power transformations are available. The default Box-Cox power family (\code{family="bcPower"}) of power transformations effectively replaces a vector by that vector raised to a power, generally in the range from -3 to 3.  For powers close to zero, the log-transformtion is suggested.  In practical situations, after estimating a power using the \code{powerTransform} function, a variable would be replaced by a simple power transformation of it, for example, if \eqn{\lambda\approx 0.5}{lamba is about 0.5}, then the correspoding variable would be replaced by its square root; if \eqn{\lambda}{\lambda} is close enough to zero, the the variable would be replaced by its natural logarithm.  The Box-Cox family requires the responses to be strictly positive.

The \code{family="bcnPower"}, or Box-Cox with negatives, family proposed by Hawkins and Weisberg (2017)  allows for (a few) non-positive values, while allowing for the transformed data to be interpreted similarly to the interpretation of Box-Cox transformed values.  This family is the Box-Cox transformation of \eqn{z = .5 * (y + (y^2 + \gamma^2)^{1/2})} that depends on a location parameter \eqn{\gamma}.  The quantity \eqn{z} is positive for all values of \eqn{y}.  If \eqn{\gamma = 0} and \eqn{y} is strictly positive, then the Box-Cox and the bcnPower transformations are identical.  When fitting the Box-Cox with negatives family, \code{lambda} is restricted to the range [-3, 3], and \code{gamma} is restricted to the range from \code{gamma.min=.1} to the largest positive value of the variable, since values outside these ranges are unreasonable in practice.  The value of \code{gamma.min} can be changed with an argument to \code{powerTransform}.

The final family \code{family="yjPower"} uses the Yeo-Johnson transformation, which is the Box-Cox transformation of \eqn{U+1} for nonnegative values, and of \eqn{|U|+1} with parameter \eqn{2-\lambda}{2-lambda} for \eqn{U} negative and thus it provides a family for fitting when (a few) observations are negative.  Because of the unusual constraints on the powers for positive and negative data, this transformation is not used very often, as results are difficult to interpret.  In practical problems, a variable would be replaced by its Yeo-Johnson transformation computed using the \code{\link{yjPower}} function.

The function \code{\link{testTransform}} is used to obtain likelihood ratio tests for any specified value for the transformation parameter(s).

Computations maximize the likelihood-like functions described by Box and Cox (1964) and by Velilla (1993).  For univariate responses, the computations are very stable and problems are unlikely, although for \sQuote{lmer} models computations may be very slow because the model is refit many times. For multivariate responses with the \code{bcnPower} family, the computing algorithm may fail.  In this case we recommend adding the argument \code{itmax = 1} to the call to \code{powerTransform}.  This will return the starting value estimates of the transformation parameters, fitting a d-dimensional response as if all the d responses were independent.

}
\value{
An object of class \code{powerTransform} or class \code{bcnPowerTransform} if \code{family="bcnPower"} that
inherits from \code{powerTransform} is returned, including the components listed below.

A \code{summary} method presents estimated values for the transformation power \code{lambda} and for the \sQuote{bcnPower} family the location parameter \code{gamma} as well.  Standard errors and Wald 95\% confidence intervals based on the standard errors are computed from the inverse of the sample Hessian matrix evaluted at the estimates.  The interval estimates for the \code{gamma} parameters will generally be very wide, reflecting little information available about the location parameter.  Likelihood ratio type tests are also provided.  For the \sQuote{bcnPower} family these are based on the profile loglikelihood for \code{lambda} alone; that is, we treat \code{gamma} as a nusiance parameter and average over it.

The components of the returned object includes
\item{lambda}{Estimated transformation parameter}
\item{roundlam}{Convenient rounded values for the estimates.  These rounded
values will  usually be the desired transformations.}
\item{gamma}{Estimated location parameters for \code{bcnPower}, \code{NULL} otherwise}
\item{invHess}{Estimated covariance matrix of the estimated parameters}
\item{llik}{Value of the log-likelihood at the estimates}

The \code{summary} method for \code{powerTransform} returns an array with columns labeled "Est Power" for the value of \code{lambda} that maximizes the likelihood; "Rounded Pwr" for \code{roundlam},  and columns "Wald Lwr Bnd" and "Wald Ur Bnd" for a 95 percent Wald normal theory confidence interval for \code{lambda} computed as the estimate plus or minus 1.96 times the standard error.
}
\references{Box, G. E. P. and Cox, D. R. (1964) An analysis of transformations.
\emph{Journal
of the Royal Statisistical Society, Series B}. 26 211-46.

Cook, R. D. and Weisberg, S. (1999)  \emph{Applied Regression Including
Computing
and Graphics}.  Wiley.

Fox, J. and Weisberg, S. (2019)
\emph{An R Companion to Applied Regression}, Third Edition, Sage.

Hawkins, D. and Weisberg, S. (2017)
Combining the Box-Cox Power and Generalized Log Transformations to Accomodate Nonpositive  Responses In Linear and Mixed-Effects Linear Models \emph{South African Statistics Journal}, 51, 317-328.

Velilla, S. (1993)  A note on the multivariate Box-Cox transformation to
normality.  \emph{Statistics and Probability Letters}, 17, 259-263.

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.

Yeo, I. and Johnson, R. (2000) A new family of
power transformations to improve normality or symmetry.
\emph{Biometrika}, 87, 954-959.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{\code{\link{testTransform}}, \code{\link{bcPower}}, \code{\link{bcnPower}}, \code{\link{transform}}, \code{\link{optim}}, \code{\link{boxCox}}. }
\examples{
# Box Cox Method, univariate
summary(p1 <- powerTransform(cycles ~ len + amp + load, Wool))
# fit linear model with transformed response:
coef(p1, round=TRUE)
summary(m1 <- lm(bcPower(cycles, p1$roundlam) ~ len + amp + load, Wool))

# Multivariate Box Cox uses Highway1 data
summary(powerTransform(cbind(len, adt, trks, sigs1) ~ 1, Highway1))

# Multivariate transformation to normality within levels of 'htype'
summary(a3 <- powerTransform(cbind(len, adt, trks, sigs1) ~ htype, Highway1))

# test lambda = (0 0 0 -1)
testTransform(a3, c(0, 0, 0, -1))

# save the rounded transformed values, plot them with a separate
# color for each highway type
transformedY <- bcPower(with(Highway1, cbind(len, adt, trks, sigs1)),
                        coef(a3, round=TRUE))
\dontrun{scatterplotMatrix( ~ transformedY|htype, Highway1) }

# With negative responses, use the bcnPower family
m2 <- lm(I1L1 ~ pool, LoBD)
summary(p2 <- powerTransform(m2, family="bcnPower"))
testTransform(p2, .5)
summary(powerTransform(update(m2, cbind(LoBD$I1L2, LoBD$I1L1) ~ .), family="bcnPower"))

\dontrun{ 
  # tajes a few seconds:
  # multivariate bcnPower, with 8 responses
  summary(powerTransform(update(m2, as.matrix(LoBD[, -1]) ~ .), family="bcnPower"))
  # multivariate bcnPower, fit with one iteration using starting values as estimates
  summary(powerTransform(update(m2, as.matrix(LoBD[, -1]) ~ .), family="bcnPower", itmax=1))
}

# mixed effects model
\dontrun{
  # uses the lme4 package
  data <- reshape(LoBD[1:20, ], varying=names(LoBD)[-1], direction="long", v.names="y")
  names(data) <- c("pool", "assay", "y", "id")
  data$assay <- factor(data$assay)
  require(lme4)
  m2 <- lmer(y ~ pool + (1|assay), data)
  summary(l2 <- powerTransform(m2, family="bcnPower", verbose=TRUE))
}
}
\keyword{ regression}% at least one, from doc/KEYWORDS
