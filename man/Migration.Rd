\name{Migration}
\alias{Migration}
\docType{data}

\title{Canadian Interprovincial Migration Data}

\description{
The \code{Migration} data frame has 90 rows and 8 columns.

}

\format{
  This data frame contains the following columns:
  \describe{
      \item{source}{
        Province of origin (source).
        A factor with levels:
        \code{ALTA}, Alberta;
        \code{BC}, British Columbia;
        \code{MAN}, Manitoba;
        \code{NB}, New Brunswick;
        \code{NFLD}, New Foundland;
        \code{NS}, Nova Scotia;
        \code{ONT}, Ontario;
        \code{PEI}, Prince Edward Island;
        \code{QUE}, Quebec;
        \code{SASK}, Saskatchewan.
        }
      \item{destination}{
        Province of destination (1971 residence).    
        A factor with levels:
        \code{ALTA}, Alberta;
        \code{BC}, British Columbia;
        \code{MAN}, Manitoba;
        \code{NB}, New Brunswick;
        \code{NFLD}, New Foundland;
        \code{NS}, Nova Scotia;
        \code{ONT}, Ontario;
        \code{PEI}, Prince Edward Island;
        \code{QUE}, Quebec;
        \code{SASK}, Saskatchewan.
        }
      \item{migrants}{
        Number of migrants (from source to destination) in the
        period 1966--1971.
        }
      \item{distance}{
        Distance (between principal cities of provinces):
        NFLD, St. John;
        PEI, Charlottetown;
        NS, Halifax;
        NB, Fredricton;
        QUE, Montreal;
        ONT, Toronto;
        MAN, Winnipeg;
        SASK, Regina;
        ALTA, Edmonton;
        BC, Vancouver.
        }
      \item{pops66}{
        1966 population of source province.
        }
      \item{pops71}{
        1971 population of source province.
        }
      \item{popd66}{
        1966 population of destination province.
        }
      \item{popd71}{
        1971 population of destination province.
        }
      }
    }

\details{
There is one record in the data file for each migration stream. You can average the 1966 
and 1971 population figures for each of the source and destination 
provinces.
}

\source{
 Canada (1962)
 \emph{Map}. Department of Mines and Technical Surveys.
 
 Canada (1971)
 \emph{Census of Canada}.
  Statistics Canada,
  Vol. 1, Part 2 [Table 32]. 
  
  Canada (1972)
  \emph{Canada Year Book}. Statistics Canada [p. 1369].
}

\references{
 Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
}

\usage{
Migration
}
\keyword{datasets}
