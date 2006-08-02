\name{levene.test}
\alias{levene.test}

\title{Levene's Test}
\description{
 Computes Levene's test for homeogeneity of variance across groups.
}
\usage{
levene.test(y, group)
}

\arguments{
  \item{y}{response variable.}
  \item{group}{factor defining groups.}
}

\value{
    returns an object meant to be printed showing the results of the test.
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}
     
\author{John Fox \email{jfox@mcmaster.ca}}
\note{adapted from a response posted by Brian Ripley to the R-help email list.}

\examples{
attach(Moore)
levene.test(conformity, fcategory)
## Levene's Test for Homogeneity of Variance
##       Df F value Pr(>F)
## group  2   0.046 0.9551
##       42  

levene.test(conformity, interaction(fcategory, partner.status))
## Levene's Test for Homogeneity of Variance
##       Df F value Pr(>F)
## group  5  1.4694 0.2219
##       39 
}        

\keyword{htest}
