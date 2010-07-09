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
    id.method = list(abs(residuals(model, type = "pearson")), "x"), 
    labels, 
    id.n = if(id.method[1]=="identify") Inf else 0,
    id.cex = 1, id.col=palette()[1], 
    line = TRUE, smooth = TRUE, span = 0.5, iter, las = par("las"), 
    col = palette()[2], pch = 1, lwd = 2, main = "Ceres Plot", 
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
  \item{layout}{ If set to a value like \code{c(1,1)} or \code{c(4,3)}, the layout
of the graph will have this many rows and columns.  If not set, the program will
select an appropriate layout.  If the number of graphs exceed nine, you must 
select the layout yourself, or you will get a maximum of nine per page.}
  \item{ask}{If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE}, the default, don't ask.
This is relevant only if not all the graphs can be drawn in one window.}
  \item{main}{The title of the plot; if missing, one will be supplied.
}
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
  \item{smooth}{\code{TRUE} to plot nonparametric-regression (lowess) line.
}
  \item{span}{span for lowess smoother.
}
  \item{iter}{number of robustness iterations for nonparametric-regression smooth;
    defaults to 3 for a linear model and to 0 for a non-Gaussian glm.
}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).
}
  \item{col}{color for points and lines; the default is the \emph{second} entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}).
}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).
}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).
  }
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
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
  
  Weisberg, S. (2005)  \emph{Applied Linear Regression}, Third Edition. Wiley.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{crPlots}}, \code{\link{avPlots}}, \code{\link{showLabels}}}

\examples{
ceresPlots(lm(prestige~income+education+type, data=Prestige), terms= ~ . - type)
}

\keyword{hplot}
\keyword{regression}
