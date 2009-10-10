\name{Transformation Axes}
\alias{power.axis}
\alias{box.cox.axis}
\alias{prob.axis}

\title{Axes for Transformed Variables}
\description{
  These functions produce axes for the original scale of 
  transformed variables. Typically these would appear as additional
  axes to the right or
  at the top of the plot, but if the plot is produced with 
  \code{axes=FALSE}, then these functions could be used for axes below or to
  the left of the plot as well.
}
\usage{
power.axis(power, base=exp(1), side=c("right", "above", "left", "below"), 
  at, grid=FALSE, grid.col=gray(0.5), grid.lty=3, 
  axis.title="Untransformed Data", cex=1, las=par("las"))

box.cox.axis(power, side=c("right", "above", "left", "below"), 
  at, grid=FALSE, grid.col=gray(0.5), grid.lty=3, 
  axis.title="Untransformed Data", cex=1, las=par("las"))

prob.axis(at, side=c("right", "above", "left", "below"), grid=FALSE, grid.lty=3, 
  grid.col=gray(0.5), axis.title="Probability", interval=0.1, cex=1, las=par("las"))
}

\arguments{
  \item{power}{power for Box-Cox or power transformation.}
  \item{side}{side at which the axis is to be drawn; numeric
   codes are also permitted: \code{side = 1} for the bottom of the plot,
   \code{side=2} for the left side, 
   \code{side = 3} for the top, \code{side = 4} for the right side.}
  \item{at}{numeric vector giving location of tick marks on
    original scale; if missing, the function will try to pick
    nice locations for the ticks.}
  \item{grid}{if \code{TRUE} grid lines for the axis will be drawn.}
  \item{grid.col}{color of grid lines.}
  \item{grid.lty}{line type for grid lines.}
  \item{axis.title}{title for axis.}
  \item{cex}{relative character expansion for axis label.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{base}{base of log transformation for \code{power.axis}
    when \code{power = 0}.}
  \item{interval}{desired interval between tick marks on the probability
    scale.}
}

\details{
  The transformations corresponding to the three functions are as follows:
  \describe{
    \item{\code{power.axis}:}{ \eqn{x^{\prime }=x^{p}}{x' = x^p} for \eqn{p\neq 0}{p != 0}
      and \eqn{x^{\prime }=\log x}{x' = log x} for \eqn{p=0}{p = 0}.}
    \item{\code{box.cox.axis}:}{\eqn{x^{\prime }=(x^{\lambda }-1)/\lambda}{x' = (x^p - 1)/p} 
      for \eqn{\lambda \neq 0}{x != 0} and \eqn{x^{\prime }=\log x}{x' = log(x)} 
      for \eqn{\lambda =0}{p = 0}.}
    \item{\code{prob.axis}:}{logit \eqn{=\log [p/(1-p)]}{= log[p/(1 - p)]}.}
  }
  
  These functions will try to place tick marks at reasonable locations, but
  producing a good-looking graph sometimes requires some fiddling with the 
  \code{at} argument.
}
\value{
  These functions are used for their side effects: to draw axes.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{ \code{\link{box.cox}}, \code{\link{logit}}}

\examples{
UN<-na.omit(UN)
attach(UN)
par(mar=c(5, 4, 4, 4)+.1)

plot(log(gdp, 10), log(infant.mortality, 10))
power.axis(0, base=10, side="above", 
  at=c(50,200,500,2000,5000,20000),grid=TRUE, axis.title="GDP per capita")
power.axis(0, base=10, side="right",
  at=c(5,10,20,50,100), grid=TRUE, axis.title="infant mortality rate per 1000")

plot(box.cox(gdp, 0), box.cox(infant.mortality, 0))
box.cox.axis(0, side="above", 
  grid=TRUE, axis.title="GDP per capita")
box.cox.axis(0, side="right",
  grid=TRUE, axis.title="infant mortality rate per 1000")

qq.plot(logit(infant.mortality/1000))
prob.axis()

qq.plot(logit(infant.mortality/1000))
prob.axis(c(.005, .01, .02, .04, .08, .16))
}
\keyword{aplot}
