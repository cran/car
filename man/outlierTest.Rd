\name{outlierTest}
\alias{outlierTest}
\alias{outlierTest.lm}
\alias{outlierTest.glm}
\alias{print.outlierTest}

\title{Bonferroni Outlier Test}

\description{
  Reports the Bonferroni p-values for Studentized residuals in
  linear and generalized linear models, based on a t-test for
  linear models and normal-distribution test for generalized linear models.
}

\usage{
outlierTest(model, ...)

\method{outlierTest}{lm}(model, cutoff=0.05, n.max=10, order=TRUE, 
	labels=names(rstudent), ...)

\method{print}{outlierTest}(x, digits=5, ...)
}

\arguments{
  \item{model}{an \code{lm} or \code{glm} model object.}
  \item{cutoff}{observations with Bonferonni p-values exceeding
  	\code{cutoff} are not reported, unless no observations are
  	nominated, in which case the one with the largest Studentized
  	residual is reported.}
  \item{n.max}{maximum number of observations to report (default, \code{10}).}
  \item{order}{report Studenized residuals in descending order of magnitude? 
  	(default, \code{TRUE}).}
  \item{labels}{an optional vector of observation names.}
  \item{...}{arguments passed down to methods functions.}
  \item{x}{\code{outlierTest} object.}
  \item{digits}{number of digits for reported p-values.}
}

\details{
  For a linear model, p-values reported use the t distribution with degrees of
  freedom one less than the residual df for the model. For a generalized
  linear model, p-values are based on the standard-normal distribution. The Bonferroni
  adjustment multiplies the usual two-sided p-value by the number of 
  observations. The \code{lm} method works for \code{glm} objects. To show all
  of the observations set \code{cutoff=Inf} and \code{n.max=Inf}.
}

\value{
  an object of class \code{outlierTest}, which is normally just
  printed.
}

\references{
  Cook, R. D. and Weisberg, S. (1982)
  \emph{Residuals and Influence in Regression.} Chapman and Hall.
  
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.
  
  Williams, D. A. (1987)
  Generalized linear model diagnostics using the deviance and single
  case deletions. \emph{Applied Statistics} \bold{36}, 181--191.
}

\author{John Fox \email{jfox@mcmaster.ca} and Sanford Weisberg}

\examples{
outlierTest(lm(prestige ~ income + education, data=Duncan))
}

\keyword{regression}
\keyword{htest}
