\name{Chile}
\alias{Chile}
\docType{data}
\title{Voting Intentions in the 1988 Chilean Plebiscite}
\description{
The \code{Chile} data frame has 2700 rows and 8 columns.
The data are from a national survey conducted in April and May of 1988
by FLACSO/Chile. There are some missing data.
}
\format{
This data frame contains the following columns:
\describe{
\item{region}{
A factor with levels:
\code{C}, Central;
\code{M}, Metropolitan Santiago area;
\code{N}, North;
\code{S}, South;
\code{SA}, city of Santiago.
}
\item{population}{
Population size of respondent's community.
}
\item{sex}{
A factor with levels:
\code{F}, female;
\code{M}, male.
}
\item{age}{
in years.
}
\item{education}{
A factor with levels (note: out of order):
\code{P}, Primary;
\code{PS}, Post-secondary;
\code{S}, Secondary.
}
\item{income}{
Monthly income, in Pesos.
}
\item{statusquo}{
Scale of support for the status-quo.
}
\item{vote}{
a factor with levels:
\code{A}, will abstain; 
\code{N}, will vote no (against Pinochet);
\code{U}, undecided;
\code{Y}, will vote yes (for Pinochet).
}
}
}

\source{
 Personal communication from FLACSO/Chile.
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}

\usage{
Chile
}
\keyword{datasets}
