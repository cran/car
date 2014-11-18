\name{mcPlots}
\alias{mcPlots}
\alias{mcPlot}
\alias{mcPlot.lm}


\title{Draw Linear Model Marginal and Conditional Plots in Parallel or Overlaid}
\description{
  These functions take a fitted \code{lm} object and draw two plots for each regressor, one showing the marginal relationship between the response and the regressor, centered but otherwise ignoring all other regressors, and the other an added-variable plot, showing the conditional relationship of Y and X given all other regressors.  The added variable plot by default is drawn using the same xlim and ylim as the centered marginal plot to emphasize that conditioning removes variation in both the regressor and the response.
}

\usage{
mcPlots(model, terms=~., layout=NULL, ask, overlaid=TRUE, ...)

mcPlot(model, ...) 

\method{mcPlot}{lm}(model, variable,
                      id.method = list(abs(residuals(model, type="pearson")), "x"),
                      labels, 
                      id.n = if(id.method[1]=="identify") Inf else 0,
                      id.cex=1, id.col=palette()[1],
                      col.marginal="blue", col.conditional="red", col.arrows="gray",
                      pch = c(16, 1), lwd = 2, grid=TRUE,   ###removed arg main
                      ellipse=FALSE, ellipse.args=list(levels=0.5), 
                      overlaid=TRUE, new=TRUE, ...)
}
\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.
}
  \item{terms}{ A one-sided formula that specifies a subset of the predictors. 
One added-variable plot is drawn for each term.   For example, the
specification \code{terms = ~.-X3} would plot against all terms
except for \code{X3}. If this argument is a quoted name of one of the terms, the 
added-variable plot is drawn for that term only.
}
\item{variable}{A quoted string giving the name of a regressor in the model matrix for the horizontal
  axis
}
  \item{layout}{ 
If set to a value like \code{c(1, 2)} or \code{c(6, 2)}, the layout
of the graph will have this many rows and columns.  If not set, the program
will select an appropriate layout of up to c(4,2); with more than 4 regressors, graphs
will be drawing one-pair-at-a-time, with user input required between pairs of graphs
}
\item{ask}{If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE} don't ask.
}
\item{\dots}{\code{mcPlots} passes these arguments to \code{mcPlot}.  
  \code{mcPlot} passes arguments to \code{plot}. 
}
  \item{id.method,labels,id.n,id.cex,id.col}{Arguments for the labelling of
    points.  The default is \code{id.n=0} for labeling no points.  See
    \code{\link{showLabels}} for details of these arguments.
}
  \item{overlaid}{If TRUE, overlay the marginal and conditional plots; otherwise
  plot them side-by-side}
  \item{col.marginal, col.conditional}{colors for points, lines,
  ellipses in the marginal and conditional plots, respectively}
  \item{col.arrows}{color for the arrows with \code{overlaid=TRUE}}
  \item{pch}{Plotting character for marginal and conditional plots, respectively.}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{grid}{If \code{TRUE}, the default, a light-gray background grid is put on the graph.}
  \item{ellipse}{If \code{TRUE}, plot a concentration ellipse; default is \code{FALSE}.}
  \item{ellipse.args}{Arguments to pass to the \code{link{dataEllipse}} function, in the form of a list
    with named elements; e.g., \code{ellipse.args=list(robust=TRUE))} will cause the ellipse to be plotted using
    a robust covariance-matrix.}
  \item{new}{if \code{TRUE}, the default, the plot window is reset when \code{overlaid=FALSE} using \code{par{mfrow=c(1, 2)}}.  If \code{FALSE}, the layout of the plot window is not reset.  
  Users will ordinarily ignore this argument.}
}  

\details{
With an \code{lm} object, suppose the response is Y, X is a focal regressor of interest, and Z is all the remaining regressors.  This function produces two graphs.  The marginal plot is simply of Y versus X, with each variable centered around its mean.  The conditional plot is the added-variable plot of e(Y|Z) versus e(X|Z) where e(a|b) means the ordinary residuals from the regression of a on b.  If \code{overlaid=TRUE}, these two plots are overlaid in one graph, with the points in different colors.  In addition, each point in the marginal plot is joined to its value in the conditional plot by an arrow.  Least squares regression lines fit to the marginal and conditional graphs are also shown; data ellipsoids can also be added.  If \code{overlaid=FALSE}, then the two graphs are shown in side-by-side plots, and the arrows are omitted.
  
These graphs are primarily for teaching, as the marginal plot shows the relationship between Y and X ignoring Z, while the conditional is the relationship between Y and X given X.  By keeping the scales the same in both graphs the effect of conditioning on both X and Y can be visualized.

These graphs are closely related to the ARES plots proposed by Cook and Weisberg (1989).
This plot would benefit from animation.
}

\value{
  These functions are used for their side effect id producing
  plots.
}

\references{
  Cook, R. D. and Weisberg, S. (1989)
  \emph{Regression diagnostics with dynamic graphics,} Technometrics, 31, 277.

  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}}

\seealso{\code{\link{avPlots}}, \code{\link{residualPlots}}, \code{\link{crPlots}}, \code{\link{ceresPlots}}, 
\code{\link{dataEllipse}}
}

\examples{
m1 <- lm(partic ~ tfr + menwage + womwage + debt + parttime, data = Bfox)
mcPlot(m1, "womwage")
mcPlot(m1, "womwage", overlaid=FALSE, ellipse=TRUE)
}

\keyword{hplot}
\keyword{regression}
