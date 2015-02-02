\name{AMSsurvey}
\alias{AMSsurvey}
\docType{data}
\title{
American Math Society Survey Data
}
\description{
Counts of new PhDs in the mathematical sciences for 2008-09 and 2011-12
categorized by type of institution, gender, and US citizenship status.
}
\usage{AMSsurvey}
\format{
  A data frame with 24 observations on the following 5 variables.
  \describe{
    \item{type}{a factor with levels \code{I(Pu)} for group I
public universities, \code{I(Pr)} for group I private universities, \code{II}
and \code{III} for groups II and III, \code{IV} for statistics and
biostatistics programs, and \code{Va} for applied mathemeatics programs.}
    \item{sex}{a factor with levels \code{Female}, \code{Male} of the
recipient}
    \item{citizen}{a factor with levels \code{Non-US}, \code{US} giving
citizenship status}
    \item{count}{The number of individuals of each type in 2008-09}
    \item{count11}{The number of individuals of each type in 2011-12}
  }
}
\details{
These data are produced yearly by the American Math Society.
}
\source{
\url{http://www.ams.org/employment/surveyreports.html}
Supplementary Table 4 in the 2008-09 data.
}
\references{
 Fox, J. and Weisberg, S. (2011) 
 \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
Phipps, Polly, Maxwell, James W. and Rose, Colleen
(2009), \emph{2009 Annual Survey of the Mathematical Sciences},
57, {250--259}, Supplementary Table 4,
\url{http://www.ams.org/employment/2009Survey-First-Report-Supp-Table4.pdf}
}

\keyword{datasets}
