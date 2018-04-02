%-------------------------------------------------------------------------------
% Revision history:
% checked in 2008-12-29 by J. Fox (corresponds to version 1.2-10 of car)
% 2009-09-16: added argument singular.ok to lm method. J. Fox
% 2012-04-04: weighted lm now allowed. John
% 2015-07-13: removed URL causing note in R CMD check. John
% 2015-08-05: fixed typo. John
%-------------------------------------------------------------------------------

\name{hccm}
\alias{hccm}
\alias{hccm.lm}
\alias{hccm.default}

\title{Heteroscedasticity-Corrected Covariance Matrices}
\description{
  Calculates heteroscedasticity-corrected covariance matrices 
  linear models fit by least squares or weighted least squares. 
  These are also called \dQuote{White-corrected} or \dQuote{White-Huber}
  covariance matrices.
}
\usage{
hccm(model, ...)

\method{hccm}{lm}(model, type=c("hc3", "hc0", "hc1", "hc2", "hc4"), 
	singular.ok=TRUE, ...)

\method{hccm}{default}(model, ...)
}

\arguments{
  \item{model}{a unweighted or weighted linear model, produced by \code{lm}.}
  \item{type}{one of \code{"hc0"}, \code{"hc1"}, \code{"hc2"}, \code{"hc3"}, or \code{"hc4"}; the
    first of these gives the classic White correction. The \code{"hc1"}, \code{"hc2"}, and \code{"hc3"}
    corrections are described in Long and Ervin (2000); \code{"hc4"} is described in Cribari-Neto (2004).}
  \item{singular.ok}{if \code{FALSE} (the default is \code{TRUE}), a model with aliased coefficients
    produces an error; otherwise, the aliased coefficients are ignored in the coefficient covariance
    matrix that's returned.}
  \item{...}{arguments to pass to \code{hccm.lm}.}
}

\details{
  The original White-corrected coefficient covariance matrix (\code{"hc0"}) for an unweighted model is
  \deqn{V(b)=(X^{\prime }X)^{-1}X^{\prime }diag(e_{i}^{2})X(X^{\prime }X)^{-1}}{V(b) = inv(X'X) X' diag(e^2) X inv(X'X)}
  where \eqn{e_{i}^{2}}{e^2} are the squared residuals, and \eqn{X} is the model
  matrix. The other methods represent adjustments to this formula. If there are weights, these are incorporated in the
  corrected covariance matrix.
  
  The function \code{hccm.default} simply catches non-\code{lm} objects.
}
\value{
  The heteroscedasticity-corrected covariance matrix for the model.
}
\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.

  Cribari-Neto, F. (2004)
  Asymptotic inference under heteroskedasticity of unknown form.
  \emph{Computational Statistics and Data Analysis} \bold{45}, 215--233.
  
  Long, J. S. and Ervin, L. H. (2000)
  Using heteroscedasity consistent standard errors in the linear regression model. 
  \emph{The American Statistician} \bold{54}, 217--224.
  
  White, H. (1980)
  A heteroskedastic consistent covariance matrix estimator and a direct test of heteroskedasticity.
  \emph{Econometrica} \bold{48}, 817--838.
  }
\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
options(digits=4)
mod<-lm(interlocks~assets+nation, data=Ornstein)
vcov(mod)
##             (Intercept)     assets  nationOTH   nationUK   nationUS
## (Intercept)   1.079e+00 -1.588e-05 -1.037e+00 -1.057e+00 -1.032e+00
## assets       -1.588e-05  1.642e-09  1.155e-05  1.362e-05  1.109e-05
## nationOTH    -1.037e+00  1.155e-05  7.019e+00  1.021e+00  1.003e+00
## nationUK     -1.057e+00  1.362e-05  1.021e+00  7.405e+00  1.017e+00
## nationUS     -1.032e+00  1.109e-05  1.003e+00  1.017e+00  2.128e+00
hccm(mod)             
##             (Intercept)     assets  nationOTH   nationUK   nationUS
## (Intercept)   1.664e+00 -3.957e-05 -1.569e+00 -1.611e+00 -1.572e+00
## assets       -3.957e-05  6.752e-09  2.275e-05  3.051e-05  2.231e-05
## nationOTH    -1.569e+00  2.275e-05  8.209e+00  1.539e+00  1.520e+00
## nationUK     -1.611e+00  3.051e-05  1.539e+00  4.476e+00  1.543e+00
## nationUS     -1.572e+00  2.231e-05  1.520e+00  1.543e+00  1.946e+00
}

\keyword{regression}
