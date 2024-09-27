\name{pointLabel}
\alias{pointLabel}
\title{ Label placement for points to avoid overlaps  }
\description{
  Use optimization routines to find good locations for point labels
  without overlaps. 
}
\usage{
pointLabel(x, y = NULL, labels = seq(along = x), cex = 1,
           method = c("SANN", "GA"),
           allowSmallOverlap = FALSE,
           trace = FALSE,
           doPlot = TRUE,
           ...)
}
\arguments{
  \item{x, y}{ as with \code{plot.default}, these provide the x and y coordinates for
          the point labels.  Any reasonable way of defining the coordinates is
          acceptable.  See the function \code{xy.coords} for details. }
  \item{labels}{ as with \code{text}, a character vector or expression specifying the text to be
          written.  An attempt is made to coerce other language objects
          (names and calls) to expressions, and vectors and other
          classed objects to character vectors by \code{as.character}. }
  \item{cex}{ numeric character expansion factor as with \code{text}. }
  \item{method}{ the optimization method, either \dQuote{SANN} for simulated
          annealing (the default) or \dQuote{GA} for a genetic algorithm. }
  \item{allowSmallOverlap}{ logical; if \code{TRUE}, labels are allowed
          a small overlap. The overlap allowed is 2\% of the diagonal
          distance of the plot area.}
  \item{trace}{ logical; if \code{TRUE}, status updates are given as the optimization algorithms
          progress.}
  \item{doPlot}{ logical; if \code{TRUE}, the labels are plotted on the
          existing graph with \code{text}.}
  \item{\ldots}{ arguments passed along to \code{text} to specify
          labeling parameters such as \code{col}. }

}
\details{
  
  Eight positions are candidates for label placement, either
  horizontally, vertically, or diagonally offset from the points. The
  default position for labels is the top right diagonal relative to the
  point (considered the preferred label position). 

  With the default settings, simulating annealing solves faster than the
  genetic algorithm. It is an open question as to which settles into a
  global optimum the best (both algorithms have parameters that may be
  tweaked).

  The label positioning problem is NP-hard (nondeterministic
  polynomial-time hard). Placement becomes difficult and slows
  considerably with large numbers of points. This function places all
  labels, whether overlaps occur or not. Some placement algorithms
  remove labels that overlap.
  
  Note that only \code{cex} is used to calculate string width and
  height (using \code{strwidth} and \code{strheight}), so passing a
  different font may corrupt the label dimensions. You could get around
  this by adjusting the font parameters with \code{par} prior to running
  this function.
  
}
\value{
  An \code{xy} list giving the \code{x} and \code{y} positions of the
  label as would be placed by \code{text(xy, labels)}.
}

\section{Note}{This function was moved from the \pkg{maptools} package in anticipation of the retirement of that package, and with the permission of the function author.}

\references{

\url{https://en.wikipedia.org/wiki/Automatic_label_placement}

Old URL: \code{https://i11www.iti.uni-karlsruhe.de/map-labeling/bibliography/}

\url{http://www.eecs.harvard.edu/~shieber/Projects/Carto/carto.html}

The genetic algorithm code was adapted from the python code at 

\url{https://meta.wikimedia.org/wiki/Map_generator}.

The simulated annealing code follows the algorithm and guidelines in:

Jon Christensen, Joe Marks, and Stuart Shieber. Placing text labels on
maps and diagrams. In Paul Heckbert, editor, Graphics Gems IV, pages
497-504. Academic Press, Boston, MA, 1994.
\url{http://www.eecs.harvard.edu/~shieber/Biblio/Papers/jc.label.pdf}

}
\author{ Tom Short, EPRI, \email{tshort@epri.com} }
\seealso{ \code{\link{text}}, \code{\link[plotrix]{thigmophobe.labels}}
  in package \pkg{plotrix} }
\examples{
n <- 50
x <- rnorm(n)*10
y <- rnorm(n)*10
plot(x, y, col = "red", pch = 20)
pointLabel(x, y, as.character(round(x,5)), offset = 0, cex = .7)

plot(x, y, col = "red", pch = 20)
pointLabel(x, y, expression(over(alpha, beta[123])), offset = 0, cex = .8)

}
\keyword{aplot}
