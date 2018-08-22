\name{showLabels}
\Rdversion{1.1}
\alias{showLabels}
\title{
 Functions to Identify and Mark Extreme Points in a 2D Plot.
}
\description{
This function is called by several graphical functions in the \pkg{car}
package to mark extreme points in a 2D plot.  Although the user is unlikely
to call this function directly, the documentation below applies to all
these other functions.
}
\usage{
showLabels(x, y, labels=NULL, method="identify",
  n = length(x), cex=1, col=carPalette()[1], location=c("lr", "ab", "avoid"), ...)

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
Plotting labels.  When called from within a \pkg{car} plotting function, the labels are automatically obtained from the row names in the data frame used to create the modeling object. If \code{labels=NULL}, case numbers will be used.  If labels are long, the
\code{\link{substr}} or \code{\link{abbreviate}} functions can be used to shorten them. Users may supply their own labels as a character vector of length equal to the number of plotted points.  For use with \pkg{car} plotting functions, the number of plotted points is equal to the number of rows of data that have neither missing values nor are excluded using the `subset' argument.  When called directly, the length of labels shoud equal the length of x.
}
  \item{method}{
How points are to be identified.  See Details below.
}
  \item{n}{
Number of points to be identified.  If set to 0, no points are identified.
}
  \item{cex}{
Controls the size of the plotted labels.  The default is \code{1}.
}
  \item{col}{
Controls the color of the plotted labels. The default is the first element returned by \code{carPalette()}.
}
  \item{location}{Where should the label be drawn?  The default is \code{"lr"} to draw the label to the left of the point for points in the right-half of the graph and to the right for points in the left-half.  The other option is \code{"ab"} for above the point for points below the middle of the graph and above the point below the middle. Finally, \code{"avoid"} tries to avoid over-plotting labels.
}
  \item{...}{not used.}
}
\details{
The argument \code{method} determine how the points
to be identified are selected.  For the default value of \code{method="identify"},
the \code{\link{identify}} function is used to identify points
interactively using the mouse.  Up to \code{n} points can be identified,
so if  \code{n=0}, which is the default in many functions in the \pkg{car}
package, then no point identification is done.

Automatic point identification can be done depending on the value of the
argument \code{method}.
\itemize{
   \item \code{method = "x"} select points according to their value of \code{abs(x - mean(x))}
   \item \code{method = "y"} select points according to their value of \code{abs(y - mean(y))}
   \item \code{method = "r"} select points according to their value of \code{abs(y)}, as may be
appropriate in residual plots, or others with a meaningful origin at 0
   \item \code{method = "mahal"} Treat \code{(x, y)} as if it were a bivariate sample, and
select cases according to their Mahalanobis distance from  \code{(mean(x), mean(y))}
  \item \code{method} can be a vector of the same length as \code{x} consisting of
values to determine the points to be labeled.  For example, for a linear model
\code{m}, setting \code{method=cooks.distance(m)} will label the
points corresponding to the largest values of Cook's distance, or
\code{method = which(abs(residuals(m, type="pearson")) > 2} would label
all observations with Pearson residuals greater than 2 in absolute value.
Warning:  If missing data are present, points may be incorrectly labelled.
  \item \code{method} can be a vector of case numbers or case-labels, in which case
  those cases will be labeled.  Warning:  If missing data are present, a list of
  case numbers may identify the wrong points.  A list of case labels, however,
  will work correctly with missing values.
  \item \code{method = "none"} causes no point labels to be shown.
}

With \code{showLabels}, the \code{method} argument can be a list, so, for
example \code{method=list("x", "y")} would label according to the horizontal
and vertical axes variables.

Finally, if the axes in the graph are logged, the function uses logged-variables
where appropriate.
}
\value{
A function primarily used for its side-effect of drawing
point labels on a plot.
Returns invisibly the labels of the selected points, or NULL if no
points are selected. Although intended for use with other functions in
the \pkg{car} package, this function can be used directly.
}
\references{
  Fox, J. and Weisberg, S. (2019) \emph{An R Companion to Applied Regression},
    Third Edition, Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}}

\seealso{\code{\link{avPlots}}, \code{\link{residualPlots}},
\code{\link{crPlots}}, \code{\link{leveragePlots}}
}

\examples{
plot(income ~ education, Prestige)
with(Prestige, showLabels(education, income,
     labels = rownames(Prestige), method=list("x", "y"), n=3))
m <- lm(income ~ education, Prestige)
plot(income ~ education, Prestige)
abline(m)
with(Prestige, showLabels(education, income,
     labels=rownames(Prestige), method=abs(residuals(m)), n=4))
}

\keyword{ utilities }

