\name{Ornstein}
\alias{Ornstein}
\docType{data}
\title{Interlocking Directorates Among Major Canadian Firms}

\description{
The \code{Ornstein} data frame has 248 rows and 4 columns.
The observations are the 248 largest Canadian firms with
publicly available information in the mid-1970s. The names
of the firms were not available.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{assets}{
      Assets in millions of dollars.
      }
    \item{sector}{
      Industrial sector.
      A factor with levels:
      \code{AGR}, agriculture, food, light industry;
      \code{BNK}, banking; 
      \code{CON}, construction; 
      \code{FIN}, other financial;
      \code{HLD}, holding companies; 
      \code{MAN}, heavy manufacturing; 
      \code{MER}, merchandizing;
      \code{MIN}, mining, metals, etc.;
      \code{TRN}, transport;
      \code{WOD}, wood and paper.
      }
    \item{nation}{
      Nation of control.
      A factor with levels:
      \code{CAN}, Canada;
      \code{OTH}, other foreign;
      \code{UK}, Britain;
      \code{US}, United States.
      }
    \item{interlocks}{
      Number of interlocking director and executive positions shared with
      other major firms.
      }
    }
  }

\source{
 Ornstein, M. (1976) 
 The boards and executives of the largest Canadian corporations. 
 \emph{Canadian Journal of Sociology} \bold{1}, 411--437.
 
 Personal communication from M. Ornstein, Department of Sociology,
 York University.
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}


\usage{
Ornstein
}

\keyword{datasets}
