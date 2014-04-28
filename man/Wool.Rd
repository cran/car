\name{Wool}
\alias{Wool}
\docType{data}
\title{Wool data}
\description{
This is a three-factor experiment with each factor at three levels, for
a total of 27 runs.  Samples of worsted yarn were with different levels of
the three factors were given a cyclic load until the sample failed.  The
goal is to understand how cycles to failure depends on the factors.
}
\format{This data frame contains the following columns:
 \describe{
    \item{len}{
 length of specimen (250, 300, 350 mm)
    }
    \item{amp}{
 amplitude of loading cycle (8, 9, 10 min)
    }
    \item{load}{
 load (40, 45, 50g)
    }
    \item{cycles}{
 number of cycles until failure
    }
 }
}

\source{
Box, G. E. P. and Cox, D. R. (1964).  An analysis of transformations (with
discussion).  \emph{J. Royal Statist. Soc.}, B26, 211-46.}
\references{
Fox, J. and Weisberg, S. (2011) 
\emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition,  Wiley, Section 6.3.
}
\usage{
Wool
}
\keyword{datasets}
