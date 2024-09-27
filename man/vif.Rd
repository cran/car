\name{vif}
\alias{vif}
\alias{vif.default}
\alias{vif.lm}
\alias{vif.merMod}
\alias{vif.polr}
\alias{vif.svyolr}

\title{Variance Inflation Factors}
\description{
  Calculates variance-inflation and generalized variance-inflation factors
  (VIFs and GVIFs) for linear, generalized linear, and other regression models.
}
\usage{
vif(mod, ...)

\method{vif}{default}(mod, ...)

\method{vif}{lm}(mod, type=c("terms", "predictor"), ...)

\method{vif}{merMod}(mod, ...)

\method{vif}{polr}(mod, ...)

\method{vif}{svyolr}(mod, ...)
}
 
\arguments{
  \item{mod}{for the default method, an object that responds to 
  \code{\link{coef}}, \code{\link{vcov}}, and
  \code{\link{model.matrix}}, such as a \code{glm} object.}
  \item{type}{for unweighted \code{lm} objects only, how to handle models that contain interactions: see Details below.}
  \item{\dots}{not used.}
}
\details{
  If all terms in an unweighted linear model have 1 df, then the usual variance-inflation
  factors are calculated.
  
  If any terms in an unweighted linear model have more than 1 df, then generalized variance-inflation factors
  (Fox and Monette, 1992) are calculated. These are interpretable as the inflation
  in size of the confidence ellipse or ellipsoid for the coefficients of the term in
  comparison with what would be obtained for orthogonal data. 
  
  The generalized VIFs
  are invariant with respect to the coding of the terms in the model (as long as
  the subspace of the columns of the model matrix pertaining to each term is
  invariant). To adjust for the dimension of the confidence ellipsoid, the function
  also prints \eqn{GVIF^{1/(2\times df)}}{GVIF^[1/(2*df)]} where \eqn{df} is the degrees of freedom
  associated with the term.
  
  Through a further generalization, the implementation here is applicable as well to other sorts of models,
  in particular weighted linear models,  generalized linear models, and mixed-effects models: Rather than computing the GVIFs from the correlation matrix of the columns of the model matrix (disregarding the regression constant), \code{vif()} uses the correlation matrix of the coefficients. Using either correlation matrix produces the same result for a linear model fit by least-squares (as pointed out by Henric Nilsson; see Fox, 2016, Exercise 13.5), but using the correlation matrix of the coefficients works with other regression models, such as GLMs, while using the correlations among the columns of the model matrix does not.
  
  Two methods of computing GVIFs are provided for unweighted linear models: 
  \itemize{
    \item Setting \code{type="terms"} (the default) behaves like the default method, and computes the GVIF for each term in the model, ignoring relations of marginality among the terms in models with interactions. GVIFs computed in this manner aren't generally sensible.
    
    \item Setting \code{type="predictor"} focuses in turn on each predictor in the model, combining the main effect for that predictor with the main effects of the predictors with which the focal predictor interacts and the interactions; e.g., in the model with formula \code{y ~ a*b + b*c}, the GVIF for the predictor \code{a} also includes the \code{b} main effect and the \code{a:b} interaction regressors; the GVIF for the predictor \code{c} includes the \code{b} main effect and the  \code{b:c} interaction; and the GVIF for the predictor \code{b} includes the \code{a} and \code{c} main effects and the \code{a:b} and \code{a:c} interactions (i.e., the whole model), and is thus necessarily 1. These predictor GVIFs should be regarded as experimental.
    }
  
  Specific methods are provided for ordinal regression model objects produced by \code{\link[MASS:polr]{polr}} in the \pkg{MASS} package and \code{\link[survey:svyolr]{svyolr}} in the \pkg{survey} package, which are "intercept-less"; VIFs or GVIFs for linear and similar regression models without intercepts are generally not sensible.
}
\value{
  A vector of VIFs, or a matrix containing one row for each term, and
  columns for the GVIF, df, and \eqn{GVIF^{1/(2\times df)}}{GVIF^[1/(2*df)]}, the last 
  of which is intended to be comparable across terms of different dimension.
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
\author{John Fox \email{jfox@mcmaster.ca} and Henric Nilsson}

\examples{
vif(lm(prestige ~ income + education, data=Duncan))
vif(lm(prestige ~ income + education + type, data=Duncan))
vif(lm(prestige ~ (income + education)*type, data=Duncan),
    type="terms") # not recommended
vif(lm(prestige ~ (income + education)*type, data=Duncan),
    type="predictor")
}

\keyword{regression}
