\name{TransformationAxes}
\alias{basicPowerAxis}
\alias{bcPowerAxis}
\alias{bcnPowerAxis}
\alias{yjPowerAxis}
\alias{probabilityAxis}

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
basicPowerAxis(power, base=exp(1), 
    side=c("right", "above", "left", "below"), 
    at, start=0, lead.digits=1, n.ticks, grid=FALSE, grid.col=gray(0.50), 
    grid.lty=2,
    axis.title="Untransformed Data", cex=1, las=par("las"))

bcPowerAxis(power, side=c("right", "above", "left", "below"), 
    at, start=0, lead.digits=1, n.ticks, grid=FALSE, grid.col=gray(0.50), 
    grid.lty=2,
    axis.title="Untransformed Data", cex=1, las=par("las"))
    
bcnPowerAxis(power, shift, side=c("right", "above", "left", "below"), 
    at, start=0, lead.digits=1, n.ticks, grid=FALSE, grid.col=gray(0.50), 
    grid.lty=2,
    axis.title="Untransformed Data", cex=1, las=par("las"))
    
yjPowerAxis(power, side=c("right", "above", "left", "below"), 
	at, lead.digits=1, n.ticks, grid=FALSE, grid.col=gray(0.50), 
  grid.lty=2,
	axis.title="Untransformed Data", cex=1, las=par("las"))

probabilityAxis(scale=c("logit", "probit"), 
	side=c("right", "above", "left", "below"),
	at, lead.digits=1, grid=FALSE, grid.lty=2, grid.col=gray(0.50),
    axis.title = "Probability", interval = 0.1, cex = 1, las=par("las"))
}

\arguments{
  \item{power}{power for Box-Cox, Box-Cox with negatives, Yeo-Johnson, or simple power transformation.}
  \item{shift}{the shift (gamma) parameter for the Box-Cox with negatives family.}
  \item{scale}{transformation used for probabilities, \code{"logit"}
    (the default) or \code{"probit"}.}
  \item{side}{side at which the axis is to be drawn; numeric
   codes are also permitted: \code{side = 1} for the bottom of the plot,
   \code{side=2} for the left side, 
   \code{side = 3} for the top, \code{side = 4} for the right side.}
  \item{at}{numeric vector giving location of tick marks on
    original scale; if missing, the function will try to pick
    nice locations for the ticks.}
  \item{start}{if a \emph{start} was added to a variable (e.g., to make all
  	data values positive), it can now be subtracted from the tick labels.}
  \item{lead.digits}{number of leading digits for determining `nice' numbers 
	for tick labels (default is \code{1}.}
  \item{n.ticks}{number of tick marks; if missing, same as corresponding
  	transformed axis.}
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
    \item{\code{basicPowerAxis}:}{Simple power transformation, 
      \eqn{x^{\prime }=x^{p}}{x' = x^p} for \eqn{p\neq 0}{p != 0}
      and \eqn{x^{\prime }=\log x}{x' = log x} for \eqn{p=0}{p = 0}.}
    \item{\code{bcPowerAxis}:}{Box-Cox power transformation,
      \eqn{x^{\prime }=(x^{\lambda }-1)/\lambda}{x' = (x^p - 1)/p} 
      for \eqn{\lambda \neq 0}{x != 0} and \eqn{x^{\prime }=\log x}{x' = log(x)} 
      for \eqn{\lambda =0}{p = 0}.}
    \item{\code{bcnPowerAxis}:}{Box-Cox with negatives power transformation, the Box-Cox power transformation of \eqn{z = .5 * (y + (y^2 + \gamma^2)^{1/2})}, where \eqn{\gamma}{gamma} is strictly positive if \eqn{y}{y} includes negative values and non-negative otherwise.  The value of \eqn{z}{z} is always positive.}
    \item{\code{yjPowerAxis}:}{Yeo-Johnson power transformation, 
      for non-negative \eqn{x}{x}, the Box-Cox transformation of
      \eqn{x + 1}{x + 1}; for negative \eqn{x}{x}, the Box-Cox transformation of
      \eqn{|x| + 1}{|x| + 1} with power \eqn{2 - p}{2 - p}.}
    \item{\code{probabilityAxis}:}{logit or probit transformation,
      logit \eqn{=\log [p/(1-p)]}{= log[p/(1 - p)]}, or 
      probit \eqn{=\Phi^{-1}(p)}{= Phi^-1(p)}, where \eqn{\Phi^{-1}}{Phi^-1} is the
      standard-normal quantile function.}
  }
  
  These functions will try to place tick marks at reasonable locations, but
  producing a good-looking graph sometimes requires some fiddling with the 
  \code{at} argument.
}
\value{
  These functions are used for their side effects: to draw axes.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\seealso{\code{\link{basicPower}}, \code{\link{bcPower}}, \code{\link{yjPower}}, 
	\code{\link{logit}}.}

\examples{
UN <- na.omit(UN)
par(mar=c(5, 4, 4, 4) + 0.1) # leave space on right

with(UN, plot(log(ppgdp, 10), log(infantMortality, 10)))
basicPowerAxis(0, base=10, side="above", 
  at=c(50, 200, 500, 2000, 5000, 20000), grid=TRUE, 
  axis.title="GDP per capita")
basicPowerAxis(0, base=10, side="right",
  at=c(5, 10, 20, 50, 100), grid=TRUE, 
  axis.title="infant mortality rate per 1000")

with(UN, plot(bcPower(ppgdp, 0), bcPower(infantMortality, 0)))
bcPowerAxis(0, side="above", 
  grid=TRUE, axis.title="GDP per capita")
bcPowerAxis(0, side="right",
  grid=TRUE, axis.title="infant mortality rate per 1000")

with(UN, qqPlot(logit(infantMortality/1000)))
probabilityAxis()

with(UN, qqPlot(qnorm(infantMortality/1000)))
probabilityAxis(at=c(.005, .01, .02, .04, .08, .16), scale="probit")

qqPlot(bcnPower(Ornstein$interlocks, lambda=1/3, gamma=0.1))
bcnPowerAxis(1/3, 0.1, at=c(o=0, 5, 10, 20, 40, 80))
}
\keyword{aplot}
