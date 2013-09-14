\name{influencePlot}
\alias{influencePlot}
\alias{influencePlot.lm}
\alias{influence.plot}

\title{Regression Influence Plot}
\description{
  This function creates a \dQuote{bubble} plot of Studentized residuals by hat values, with the areas of the
  circles representing the observations proportional to Cook's distances. Vertical reference lines are drawn
  at twice and three times the average hat value, horizontal reference lines at -2, 0, and 2 on the
  Studentized-residual scale.
}

\usage{
influencePlot(model, ...)

\method{influencePlot}{lm}(model, scale=10, 
	xlab="Hat-Values", ylab="Studentized Residuals", 
    labels, id.method = "noteworthy",
    id.n = if(id.method[1]=="identify") Inf else 0, 
    id.cex=1, id.col=palette()[1], ...)
}

\arguments{
  \item{model}{a linear or generalized-linear model.
  }
  \item{scale}{a factor to adjust the size of the circles.
  }
  \item{xlab, ylab}{axis labels.
  }
  \item{labels, id.method, id.n, id.cex, id.col}{settings for labelling
points; see \code{link{showLabels}} for details. To omit point labelling, set
\code{id.n=0}, the default.  The default \code{id.method="noteworthy"} is used
only in this function and indicates setting labels for points with large
Studentized residuals, hat-values or Cook's distances.  Set 
\code{id.method="identify"} for interactive point identification.
}
  \item{\dots}{arguments to pass to the \code{plot} and \code{points} functions.
  }
}

\value{
  If points are identified, returns a data frame with the hat values, 
Studentized residuals and Cook's distance of the identified points.  If
no points are identified, nothing is returned.  This function is primarily
used for its side-effect of drawing a plot.
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}, minor changes by S. Weisberg
\email{sandy@umn.edu}}


\seealso{\code{\link[stats]{cooks.distance}}, \code{\link[stats]{rstudent}}, 
\code{\link[stats]{hatvalues}}, \code{\link{showLabels}}}

\examples{
influencePlot(lm(prestige ~ income + education, data=Duncan))
}

\keyword{regression}
