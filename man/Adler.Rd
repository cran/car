\name{Adler}
\alias{Adler}
\docType{data}

\title{Experimenter Expectations} 

\description{
The \code{Adler} data frame has 97 rows and 3 columns.  
                                                        
The ``experimenters'' were the actual subjects of the study. 
They collected ratings of the apparent successfulness of
people in pictures who were pre-selected for their
average appearance. The experimenters were told prior to collecting 
data that the pictures were either high or low in their 
appearance of success, and were instructed to get good data, 
scientific data, or were given no such instruction. Each 
experimenter collected ratings from 18 randomly assigned
respondents; a few subjects were deleted at random to produce
an unbalanced design.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{instruction}{
     a factor with levels:
     \code{GOOD}, good data;
     \code{NONE}, no stress;
     \code{SCIENTIFIC}, scientific data.
     }
    \item{expectation}{
     a factor with levels:
     \code{HIGH}, expect high ratings;
     \code{LOW}, expect low ratings.
     }
    \item{rating}{
     The average rating obtained.  
     }
   }
 }

\source{
 Adler, N. E. (1973)
 Impact of prior sets given experimenters and subjects on the experimenter
 expectancy effect.
 \emph{Sociometry} \bold{36}, 113--126.
}

\references{
  Erickson, B. H., and Nosanchuk, T. A. (1977)
  \emph{Understanding Data.} McGraw-Hill Ryerson.
}

\usage{
Adler
}

\keyword{datasets}
