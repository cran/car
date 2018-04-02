\name{vif}
\alias{vif}
\alias{vif.default}
\alias{vif.merMod}

\title{Variance Inflation Factors}
\description{
  Calculates variance-inflation and generalized variance-inflation factors
  for linear, generalized linear, and other models.
}
\usage{
vif(mod, ...)

\method{vif}{default}(mod, ...)

\method{vif}{merMod}(mod, ...)
}
 
\arguments{
  \item{mod}{for the default method, an object that responds to 
  \code{\link{coef}}, \code{\link{vcov}}, and
  \code{\link{model.matrix}}, such as an \code{lm} or \code{glm} object.}
  \item{\dots}{not used.}
}
\details{
  If all terms in an unweighted linear model have 1 df, then the usual variance-inflation
  factors are calculated.
  
  If any terms in an unweighted linear model have more than 1 df, then generalized variance-inflation factors
  (Fox and Monette, 1992) are calculated. These are interpretable as the inflation
  in size of the confidence ellipse or ellipsoid for the coefficients of the term in
  comparison with what would be obtained for orthogonal data. 
  
  The generalized vifs
  are invariant with respect to the coding of the terms in the model (as long as
  the subspace of the columns of the model matrix pertaining to each term is
  invariant). To adjust for the dimension of the confidence ellipsoid, the function
  also prints \eqn{GVIF^{1/(2\times df)}}{GVIF^[1/(2*df)]} where \eqn{df} is the degrees of freedom
  associated with the term.
  
  Through a further generalization, the implementation here is applicable as well to other sorts of models,
  in particular weighted linear models,  generalized linear models, and mixed-effects models.
}
\value{
  A vector of vifs, or a matrix containing one row for each term in the model, and
  columns for the GVIF, df, and \eqn{GVIF^{1/(2\times df)}}{GVIF^[1/(2*df)]}.
}
\references{ 
  Fox, J. and Monette, G. (1992)
  Generalized collinearity diagnostics.
  \emph{JASA}, \bold{87}, 178--183.

  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2018) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}
\author{Henric Nilsson and John Fox \email{jfox@mcmaster.ca}}

\examples{
vif(lm(prestige ~ income + education, data=Duncan))
vif(lm(prestige ~ income + education + type, data=Duncan))

}

\keyword{regression}
