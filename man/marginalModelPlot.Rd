\name{mmps}
\alias{mmps}
\alias{mmp}
\alias{mmp.lm}
\alias{mmp.glm}
\alias{mmp.default}
\alias{marginalModelPlot}
\alias{marginalModelPlots}

\title{Marginal Model Plotting}
\description{
For a regression object, plots the response on the vertical axis versus
a linear combination \eqn{u} of terms in the mean function on the horizontal
axis.  Added to the plot are a \code{loess} smooth for the graph, along with
a loess smooth from the plot of the fitted values on \eqn{u}.  \code{mmps} is an alias
for \code{marginalModelPlots}, and \code{mmp} is an alias for \code{marginalModelPlot}.
}
\usage{
marginalModelPlots(...)

mmps(model, terms= ~ ., fitted=TRUE, layout=NULL, ask,
        main, AsIs=FALSE, ...)

marginalModelPlot(...)

\method{mmp}{lm}(model, variable, mean = TRUE, sd = FALSE,
    xlab = deparse(substitute(variable)), degree = 1, span = 2/3, key=TRUE,
    col.line = palette()[c(4, 2)], 
    ...)
    
\method{mmp}{default}(model, variable, mean = TRUE, sd = FALSE, xlab =
                 deparse(substitute(variable)), degree = 1, span = 2/3,
                 key = TRUE, col.line = palette()[c(4,2)],   
                 labels, id.method = "y", 
                 id.n=if(id.method[1]=="identify") Inf else 0, 
                 id.cex = 1, id.col=palette()[1], grid=TRUE, ...)
   
\method{mmp}{glm}(model, variable, mean = TRUE, sd = FALSE,
    xlab = deparse(substitute(variable)), degree = 1, span = 2/3, key=TRUE,
    col.line = palette()[c(4, 2)],
    labels, id.method="y", 
    id.n=if(id.method[1]=="identify") Inf else 0,
    id.cex=1, id.col=palette()[1], grid=TRUE, ...)
}

\arguments{
\item{model}{A regression object, usually of class either \code{lm} or \code{glm},
   for which there is a \code{predict} method defined. }
\item{terms}{A one-sided formula.  A marginal model plot will be drawn for
  each variable on the right-side of this formula that is not a factor.  The
  default is \code{~ .}, which specifies that all the terms in 
  \code{formula(object)} will be used.  See examples below.}
\item{fitted}{If the default \code{TRUE}, then a marginal model plot in the direction
  of the fitted values or linear predictor of a generalized linear model will
  be drawn.}
\item{layout}{A reasonable layout for the plots in the window is
  determined by the program.  If you don't like the default you can set your
  own layout:  \code{c(2, 3)} means two rows and three columns.}
\item{ask}{If \code{TRUE}, ask before clearing the graph window to draw more plots.}
        \item{main}{
Main title for the array of plots.  Use \code{main=""} to suppress the title;
if missing, a title will be supplied.
}
  \item{AsIs}{
If \code{FALSE}, the default, terms that use the \dQuote{as-is} function \code{I}
are skipped; if \code{TRUE} they are included.
}
\item{\dots}{
Additional arguments passed from \code{mmps} to \code{mmp} and
then to \code{plot}.  Users should generally use \code{mmps}, or equivalently 
\code{marginalModelPlots}.
}
\item{variable}{ The quantity to be plotted on the horizontal axis.  The
  default is the predicted values \code{predict(object)}.  Can be any other
  vector of length equal to the number of observations in the object. Thus the
  \code{mmp} function can be used to get a marginal model plot versus any
  predictor or term while the \code{mmps} function can be used only to get
  marginal model plots for the first-order terms in the formula.  In 
  particular, terms defined by a spline basis are skipped by \code{mmps}, but
  you can use \code{mmp} to get the plot for the variable used to define
  the splines.}
\item{mean}{ If \code{TRUE}, compare mean smooths}
\item{sd}{ If \code{TRUE}, compare sd smooths.  For a binomial regression with all
sample sizes equal to one, this argument is ignored as the SD bounds don't
make any sense. }
\item{xlab}{ label for horizontal axis }
\item{degree}{Degree of the local polynomial, passed to \code{loess}.  The
  usual default for \code{loess} is 2, but the default here is 1.}
\item{span}{Span, the smoothing parameter for \code{loess}.  }
\item{key}{If \code{TRUE}, include a key at the top of the plot, if \code{FALSE} omit the
key}
  \item{id.method,labels,id.n,id.cex,id.col}{Arguments for labelling
points.  The default \code{id.n=0} suppresses labelling, and setting this
argument greater than zero will include labelling.  See
\code{\link{showLabels}} for these arguments.}
\item{col.line}{ colors for data and model smooth, respectively.  Using
the default palette, these are blue and red. }
\item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
}   
\details{
\code{mmp} and \code{marginalModelPlot} draw one marginal model plot against
whatever is specified as the horizontal axis.  
\code{mmps} and \code{marginalModelPlots} draws marginal model plots
versus each of the terms in the \code{terms} argument and versus fitted values.  
\code{mmps} skips factors, interactions, and spline terms.  For polynomial
terms in one variable, plots are against the first-order term.
}

\value{
Used for its side effect of producing plots. 
}

\seealso{\code{\link{loess}}, \code{\link{plot}} }

\references{ 

Fox, J. and Weisberg, S. (2011)
\emph{An R Companion to Applied Regression}, Second Edition.  Sage.


Weisberg, S. (2005) \emph{Applied
Linear Regression}, Third Edition, Wiley, Chapter 8.
}

\author{Sanford Weisberg, \email{sandy@stat.umn.edu}}

\examples{
c1 <- lm(infant.mortality ~ gdp, UN)
mmps(c1)
# include SD lines
p1 <- lm(prestige ~ income + education, Prestige)
mmps(p1, sd=TRUE)
# logisitic regression example
m1 <- glm(lfp ~ ., family=binomial, data=Mroz)
mmps(m1)

} 
\keyword{hplot }% at least one, from doc/KEYWORDS
\keyword{regression}
