\name{hccm}
\alias{hccm}
\alias{hccm.lm}
\alias{hccm.default}

\title{Heteroscedasticity-Corrected Covariance Matrices}
\description{
  Calculates heteroscedasticity-corrected covariance matrices for
  unweighted linear models. These are also called ``White-corrected''
  covariance matrices.
}
\usage{
hccm(model, ...)

hccm.lm(model, type=c("hc3", "hc0", "hc1", "hc2"))

hccm.default(model, ...)
}

\arguments{
  \item{model}{an unweighted linear model, produced by \code{lm}.}
  \item{type}{one of \code{"hc0"}, \code{"hc1"}, \code{"hc2"}, or \code{"hc3"}; the
    first of these gives the classic White correction. The others are described in
    Long and Ervin (2000).
 }
}
\details{
  The classical White-corrected coefficient covariance matrix (\code{"hc0"}) is
  \deqn{V(b)=(X^{\prime }X)^{-1}X^{\prime }diag(e_{i}^{2})X(X^{\prime }X)^{-1}}{V(b) = inv(X'X) X' diag(e^2) X inv(X'X)}
  where \eqn{e_{i}^{2}}{e^2} are the squared residuals, and \eqn{X} is the model
  matrix. The other methods represent adjustments to this formula.
  
  The function \code{hccm.default} simply catches non-\code{lm} objects.
}
\value{
  The heteroscedasticity-corrected covariance matrix for the model.
}
\references{
  Long, J. S. and Ervin, L. H. (2000) 
  Using heteroscedasity consistent standard errors in the linear regression model. 
  \emph{The American Statistician} \bold{54}, 217--224.
  
  White, H. (1980)
  A heterskedastic consistent covariance matrix estimator and a direct test of heteroskedasticity.
  \emph{Econometrica} \bold{48}, 817--838.
  }
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{ncv.test}}, \code{\link{spread.level.plot}} }

\examples{
options(digits=4)
data(Ornstein)
mod<-lm(interlocks~assets+nation, data=Ornstein)
Var(mod)
##             (Intercept)     assets  nationOTH   nationUK   nationUS
## (Intercept)   1.079e+00 -1.588e-05 -1.037e+00 -1.057e+00 -1.032e+00
## assets       -1.588e-05  1.642e-09  1.155e-05  1.362e-05  1.109e-05
## nationOTH    -1.037e+00  1.155e-05  7.019e+00  1.021e+00  1.003e+00
## nationUK     -1.057e+00  1.362e-05  1.021e+00  7.405e+00  1.017e+00
## nationUS     -1.032e+00  1.109e-05  1.003e+00  1.017e+00  2.128e+00
hccm(mod)             
##             (Intercept)     assets  nationOTH   nationUK   nationUS
## (Intercept)   1.4808521 -1.640e-05 -1.445e+00 -1.460e+00 -1.455e+00
## assets       -0.0000164  2.425e-09  1.079e-05  1.324e-05  1.161e-05
## nationOTH    -1.4451983  1.079e-05  6.743e+00  1.432e+00  1.430e+00
## nationUK     -1.4600207  1.324e-05  1.432e+00  3.722e+00  1.440e+00
## nationUS     -1.4550976  1.161e-05  1.430e+00  1.440e+00  1.840e+00
}

\keyword{regression}
