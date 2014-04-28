\name{ceresPlots}
\alias{ceresPlots}
\alias{ceresPlot}
\alias{ceresPlot.lm}
\alias{ceresPlot.glm}

\title{Ceres Plots}
\description{
  These functions draw Ceres plots for linear and generalized
  linear models.
}
\usage{
ceresPlots(model, terms = ~., layout = NULL, ask, main, 
    ...)

ceresPlot(model, ...)

\method{ceresPlot}{lm}(model, variable, 
  id.method = list(abs(residuals(model, type="pearson")), "x"),
  labels, 
  id.n = if(id.method[1]=="identify") Inf else 0,
  id.cex=1, id.col=palette()[1],
  line=TRUE, smoother=loessLine, smoother.args=list(),
    smooth, span,
	col=palette()[1], col.lines=palette()[-1],
  xlab, ylab, pch=1, lwd=2,  
  grid=TRUE, ...)

\method{ceresPlot}{glm}(model, ...)
}

\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.}
  \item{terms}{ A one-sided formula that specifies a subset of the predictors. 
One component-plus-residual plot is drawn for each term.  The default
\code{~.} is to plot against all numeric predictors.  For example, the
specification \code{terms = ~ . - X3} would plot against all predictors
except for \code{X3}.  Factors and nonstandard predictors such as B-splines are 
skipped.  If this argument is a quoted name of one of the predictors, the 
component-plus-residual plot is drawn for that predictor only.
}
  \item{layout}{ 
If set to a value like \code{c(1, 1)} or \code{c(4, 3)}, the layout
of the graph will have this many rows and columns.  If not set, the program
will select an appropriate layout.  If the number of graphs exceed nine, you 
must select the layout yourself, or you will get a maximum of nine per page.
If \code{layout=NA}, the function does not set the layout and the user can
use the \code{par} function to control the layout, for example to have
plots from two models in the same graphics window.
}
  \item{ask}{If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE}, the default, don't ask.
This is relevant only if not all the graphs can be drawn in one window.}
  \item{main}{Overall title for any array of cerers plots; if missing a default
is provided.}
  \item{\dots}{\code{ceresPlots} passes these arguments to \code{ceresPlot}.  
\code{ceresPlot} passes them to \code{plot}. 
}  
  \item{variable}{A quoted string giving the name of a variable for the horizontal
axis}
  \item{id.method,labels,id.n,id.cex,id.col}{Arguments for the labelling of
    points.  The default is \code{id.n=0} for labeling no points.  See
    \code{\link{showLabels}} for details of these arguments.
}
  \item{line}{\code{TRUE} to plot least-squares line.
}
  \item{smoother}{Function to add a nonparametric smooth.
}
  \item{smoother.args}{see \code{\link{ScatterplotSmoothers}} for available
  smooethers and arguments.
}
  \item{smooth, span}{these arguments are included for backwards compatility: 
    if \code{smooth=TRUE} then \code{smoother} is set to \code{loessLine},
    and if \code{span} is specified, it is added to \code{smoother.args}.}
  \item{col}{color for points; the default is the first entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}).
}
  \item{col.lines}{a list of at least two colors.  The first color is used for the
    ls line and the second color is used for the fitted lowess line.  To use
    the same color for both, use, for example, \code{col.lines=c("red", "red")}
}
  \item{xlab,ylab}{labels for the x and y axes, respectively.  If not set 
    appropriate labels are created by the function.}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).
}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).
}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph
}
}
\details{
  Ceres plots are a generalization of component+residual (partial
  residual) plots that are less prone to leakage of nonlinearity
  among the predictors.
  
  The function intended for direct use is \code{ceresPlots}. 
  
  The model cannot contain interactions, but can contain factors. 
  Factors may be present in the model, but Ceres plots cannot be drawn
  for them.
}
\value{
  \code{NULL}. These functions are used for their side effect: producing
  plots.
}
\references{ 
  Cook, R. D. and Weisberg, S. (1999)
  \emph{Applied Regression, Including Computing and Graphics.} Wiley.

  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Weisberg, S. (2014)  \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{crPlots}}, \code{\link{avPlots}}, \code{\link{showLabels}}}

\examples{
ceresPlots(lm(prestige~income+education+type, data=Prestige), terms= ~ . - type)
}

\keyword{hplot}
\keyword{regression}
