\name{Export}
\alias{Export}
\title{
Export a data frame to disk in one of many formats
}
\description{
Uses the \code{export} function in the \pkg{rio} package to export a file to disk.  This function adds an argument for converting row.names to a column in the resulting file.
}
\usage{
Export(x, file, format, ..., keep.row.names)
}
\arguments{
  \item{x}{
A data frame or matrix to be written to a file.
}
  \item{file}{
A character string naming a file. If the file name has an extension, such as \code{.xlsx}, the extention is used to infer the type of file to be exported.  See \code{\link[rio]{export}} for the file types supported.
}
  \item{format}{see \code{\link[rio]{export}}.
}
  \item{\dots}{
Additional arguments; see \code{\link[rio]{export}}.
}
  \item{keep.row.names}{
If set to \code{TRUE}, then the data frame's row.names are appended to the left of the data frame with the name "id".  If set to quoted character string, the row.names are added using the character string as its name.  If set to \code{FALSE} row.names are lost.
}
}
\details{
This is a convenience function in the \pkg{car} package for exporting (writing) a data frame to a file in a wide variety of formats including csv, Microsoft Excel.  It optionally allows converting the row.names for the data frame to a column before writing.  It then calls \code{\link[rio]{export}} in the \code{rio} package.  That function in turn uses many other packages and functions for writing the function to a file.
}
\value{
The name of the output file as a character string (invisibly).
}

\author{Sanford Weisberg \email{sandy@umn.edu}}

\references{
Chung-hong Chan, Geoffrey CH Chan, Thomas J. Leeper, and Jason Becker (2017). rio: A
  Swiss-army knife for data file I/O. R package version 0.5.0.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[rio]{export}}, \code{\link{Import}}
}
\examples{
Export(Duncan, "Duncan.csv", keep.row.names="occupation")
Duncan2 <- Import("Duncan.csv") # Automatically restores row.names
identical(Duncan, Duncan2)
# cleanup
unlink("Duncan.csv")
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
\keyword{ connections }% __ONLY ONE__ keyword per line
