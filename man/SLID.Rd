\name{SLID}
\alias{SLID}
\docType{data}
\title{Survey of Labour and Income Dynamics}
\description{
The \code{SLID} data frame has 7425 rows and 5 columns.
The data are from the 1994 wave of the Canadian Survey of Labour and Income Dynamics,
for the province of Ontario.
There are missing data, particularly for wages.
}
\format{
This data frame contains the following columns:
\describe{
\item{wages}{
Composite hourly wage rate from all jobs.
}
\item{education}{
Number of years of schooling.
}
\item{age}{
in years.
}
\item{sex}{
A factor with levels:
\code{Female},
\code{Male}.
}
\item{language}{
A factor with levels:
\code{English},
\code{French},
\code{Other}.
}
}
}


\source{
 The data are taken from the public-use dataset made available by
 Statistics Canada, and prepared by the Institute for Social Research,
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
SLID
}
\keyword{datasets}
