\name{levene.test}
\alias{levene.test}
\alias{levene.test.formula}
\alias{levene.test.lm}
\alias{levene.test.default}

\title{Levene's Test}
\description{
 Computes Levene's test for homeogeneity of variance across groups.
}
\usage{
levene.test(y, ...)
\method{levene.test}{formula}(y, ...)
\method{levene.test}{lm}(y, ...)
\method{levene.test}{default}(y, group, ...)
}

\arguments{
  \item{y}{response variable for the default method, or \code{lm} or 
  \code{formula} object. If \code{y} is a linear-model object or a formula,
  the variables on the right-hand-side of the model must all be factors and
  must be completely crossed.}
  \item{group}{factor defining groups.}
  \item{...}{arguments to be passed down, e.g., \code{data} for the
  \code{formula} and \code{lm} methods.}
}

\details{
  The version of Levene's test implemented here is based on absolute
  deviations from the group medians. Other, less robust, versions, not implemented,
  are based on absolute deviations from group means or trimmed means.
}

\value{
    returns an object meant to be printed showing the results of the test.
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}
     
\author{John Fox \email{jfox@mcmaster.ca}; generic version 
    contributed by Derek Ogle}
\note{adapted from a response posted by Brian Ripley to the R-help email list.}

\examples{
with(Moore, levene.test(conformity, fcategory))
## Levene's Test for Homogeneity of Variance
##       Df F value Pr(>F)
## group  2   0.046 0.9551
##       42  

with(Moore, levene.test(conformity, interaction(fcategory, partner.status)))
## Levene's Test for Homogeneity of Variance
##       Df F value Pr(>F)
## group  5  1.4694 0.2219
##       39 

levene.test(conformity ~ fcategory*partner.status, data=Moore)
levene.test(lm(conformity ~ fcategory*partner.status, data=Moore))
}        

\keyword{htest}
