\name{Highway1}
\alias{Highway1}
\docType{data}
\title{Highway Accidents}
\description{
The data comes from a unpublished master's paper by Carl Hoffstedt. They relate the automobile accident rate, 
in accidents per million vehicle miles to several potential terms.  The data include 39 sections of large 
highways in the state of Minnesota in 1973.  The goal of this analysis was to understand the impact 
of design variables, \code{Acpts}, \code{Slim}, \code{Sig}, and \code{Shld} that are under the control of 
the highway department, on accidents.
}
\format{This data frame contains the following columns:
 \describe{
    \item{rate}{
 1973 accident rate per million vehicle miles
    }
    \item{len}{
 length of the \code{Highway1} segment in miles
    }
    \item{ADT}{
 average daily traffic count in thousands
    }
    \item{trks}{
 truck volume as a percent of the total volume
    }
    \item{sigs1}{
 (number of signalized interchanges per mile times len + 1)/len, the number
of signals per mile of roadway, adjusted to have no zero values.
    }
    \item{slim}{
 speed limit in 1973
    }
    \item{shld}{
 width in feet of outer shoulder on the roadway
    }
    \item{lane}{
 total number of lanes of traffic
    }
    \item{acpt}{
 number of access points per mile
    }
    \item{itg}{
 number of freeway-type interchanges per mile
    }
    \item{lwid}{
 lane width, in feet
    }
    \item{hwy}{An indicator of the type of roadway or the source of funding
    for the road, either MC, FAI, PA, or  MA
    }
 }
}

\source{Carl Hoffstedt.  This differs from the dataset \code{highway} in the
\code{alr3} package only by transformation of some of the columns.}

\references{
Fox, J. and Weisberg, S. (2011) \emph{An R Companion to Applied
Regression}, Second Edition, Sage.

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition,  Wiley, Section 7.2.
}

\usage{
Highway1
}
\keyword{datasets}
