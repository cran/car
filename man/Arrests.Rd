\name{Arrests}
\alias{Arrests}
\docType{data}

\title{Arrests for Marijuana Possession}

\description{
Data on police treatment of individuals arrested in
Toronto for simple possession of small quantities of marijuana. The data are part of a larger
data set featured in a series of articles in the Toronto Star newspaper.
}

\usage{data(Arrests)}

\format{
  A data frame with 5226 observations on the following 8 variables.
  \describe{
    \item{released}{Whether or not the arrestee was released with a summons; 
    a factor with levels:
        \code{No};
        \code{Yes}. }
    \item{colour}{The arrestee's race; a factor with levels:
        \code{Black}; 
        \code{White}. }
    \item{year}{1997 through 2002; a numeric vector.}
    \item{age}{in years; a numeric vector.}
    \item{sex}{a factor with levels:
        \code{Female};
        \code{Male}. }
    \item{employed}{a factor with levels:
        \code{No};
        \code{Yes}. }
    \item{citizen}{a factor with levels:
        \code{No};
        \code{Yes}. }
    \item{checks}{Number of police data bases (of previous
    arrests, previous convictions, parole status, etc. -- 6 in all) on which
    the arrestee's name appeared; a numeric vector}
    }
}

\source{
  Personal communication from Michael Friendly, York University.
}

\keyword{datasets}
