\name{showLabels}
\Rdversion{1.1}
\alias{showLabels}
\title{
Utility Functions to Identify and Mark Extreme Points in a 2D Plot.  
}
\description{
This function is called by several graphical functions in the \code{car}
package to mark extreme points in a 2D plot.  Although the user is unlikely
to call this function directly, the documentation below applies to all
these other functions.
}
\usage{
showLabels(x, y, labels=NULL, id.method="identify",  
  id.n = length(x), id.cex=1, id.col=palette()[1], ...)
     
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Plotted horizontal coordinates.
}
  \item{y}{
Plotted vertical coordinates.
}
  \item{labels}{
Plotting labels.  If \code{NULL}, case numbers will be used.  If labels are long, the
\code{\link{substr}} or \code{\link{abbreviate}} function can be used to shorten them.
}
  \item{id.method}{
How points are to be identified.  See Details below.
}
  \item{id.n}{
Number of points to be identified.  If set to zero, no points are identified.  
}
  \item{id.cex}{
Controls the size of the plotted labels.  The default is \code{1}.
}
  \item{id.col}{
Controls the color of the plotted labels. 
}
  \item{...}{additional arguments passed to \code{identify} or to \code{text}.}
}
\details{
The argument \code{id.method} determine how the points 
to be identified are selected.  For the default value of \code{id.method="identify"},
the \code{\link{identify}} function is used to identify points 
interactively using the mouse.  Up to \code{id.n} points can be identified, 
so if  \code{id.n=0}, which is the default in many functions in the \code{car} 
package, then no point identification is done.

Automatic point identification can be done depending on the value of the
argument \code{id.method}.
\itemize{
   \item \code{id.method = "x"} select points according to their value of \code{abs(x - mean(x))}
   \item \code{id.method = "y"} select points according to their value of \code{abs(y - mean(y))}
   \item \code{id.method = "mahal"} Treat \code{(x, y)} as if it were a bivariate sample, and
select cases according to their Mahalanobis distance from  \code{(mean(x), mean(y))}
  \item \code{id.method} can be a vector of the same length as \code{x} consisting of 
values to determine the points to be labeled.  For example, for a linear model
\code{m}, setting \code{id.method=cooks.distance(m), id.n=4} will label the 
points corresponding to the four largest values of Cook's distance, or
\code{id.method = abs(residuals(m, type="pearson")), id.n=2} would label
the two observations corresponding to the largest absolute Pearson residuals.
  \item \code{id.method} can be a vector of case numbers or case-labels, in which case
  those cases will be labeled, as long as \code{id.n} is greater than zero.
}

With \code{showLabels}, the \code{id.method} argument can be list, so, for
example \code{id.method=list("x", "y")} would label according to the horizontal 
and vertical axes variables.

Finally, if the axes in the graph are logged, the function uses logged-variables
where appropriate.
}
\value{
A utility function used for its side-effect of drawing labels on a plot. 
Although intended for use with other functions in the \code{car} package, 
this function can be used directly.
}
\references{
  Fox, J. and Weisberg, S. (2011) \emph{An R Companion to Applied Regression},
    Second Edition, Sage.
  
  Weisberg, S. (2005) \emph{Applied Linear Regression}, Third Edition, Wiley.
}
\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}}

\seealso{\code{\link{avPlots}}, \code{\link{residualPlots}},
\code{\link{crPlots}}, \code{\link{leveragePlots}}
}

\examples{
plot(income ~ education, Prestige)
with(Prestige, showLabels(education, income, 
     labels = rownames(Prestige), id.method=list("x", "y"), id.n=3))
m <- lm(income ~ education, Prestige)
plot(income ~ education, Prestige)
abline(m)
with(Prestige, showLabels(education, income, 
     labels=rownames(Prestige), id.method=abs(residuals(m)), id.n=4))
}

\keyword{ utilities }

