\name{scatter3d}
\alias{scatter3d}
\alias{scatter3d.formula}
\alias{scatter3d.default}
\alias{Identify3d}
\title{Three-Dimensional Scatterplots and Point Identification}

\description{
  The \code{scatter3d} function uses the \code{rgl} package to draw 3D scatterplots
  with various regression surfaces. The function \code{Identify3d}
  allows you to label points interactively with the mouse:
  Press the right mouse button (on a two-button mouse) or the centre button (on a
  three-button mouse), drag a
  rectangle around the points to be identified, and release the button.
  Repeat this procedure for each point or
  set of \dQuote{nearby} points to be identified. To exit from point-identification mode,
  click the right (or centre) button in an empty region of the plot.
}

\usage{
scatter3d(x, ...)

\method{scatter3d}{formula}(formula, data, subset, radius, xlab, ylab, zlab, id=FALSE, ...)

\method{scatter3d}{default}(x, y, z,
	  xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
	  zlab=deparse(substitute(z)), axis.scales=TRUE, axis.ticks=FALSE,
	  revolutions=0,
	  bg.col=c("white", "black"),
	  axis.col=if (bg.col == "white") c("darkmagenta", "black", "darkcyan")
		            else c("darkmagenta", "white", "darkcyan"),
	  surface.col=carPalette()[-1], surface.alpha=0.5,
	  neg.res.col="magenta", pos.res.col="cyan",
	  square.col=if (bg.col == "white") "black" else "gray",
	  point.col="yellow", text.col=axis.col,
	  grid.col=if (bg.col == "white") "black" else "gray",
	  fogtype=c("exp2", "linear", "exp", "none"),
	  residuals=(length(fit) == 1),
	  surface=TRUE, fill=TRUE,
	  grid=TRUE, grid.lines=26, df.smooth=NULL, df.additive=NULL,
	  sphere.size=1, radius=1, threshold=0.01, speed=1, fov=60,
	  fit="linear", groups=NULL, parallel=TRUE,
	  ellipsoid=FALSE, level=0.5, ellipsoid.alpha=0.1, id=FALSE,
	  model.summary=FALSE, ...)

Identify3d(x, y, z, axis.scales=TRUE, groups = NULL, labels = 1:length(x),
	col = c("blue", "green", "orange", "magenta", "cyan", "red", "yellow", "gray"),
	offset = ((100/length(x))^(1/3)) * 0.02)
}

\note{
You have to install the \code{rgl} package to produce 3D plots.  On a Macintosh (but not on Windows or Linux), you may also need to install the X11 windowing system.  Go to \url{https://www.xquartz.org/} and click on the link for XQuartz. Double-click on the downloaded disk-image file, and then double-click on \code{XQuartz.pkg} to start the installer. You may take all of the defaults in the installation. After XQuartz is installed, you should restart your Macintosh.
}

\arguments{
  \item{formula}{``model'' formula, of the form \code{y ~ x + z} or
    to plot by groups \code{y ~ x + z | g}, where \code{g} evaluates to a factor
    or other variable dividing the data into groups.}
  \item{data}{data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{x}{variable for horizontal axis.}
  \item{y}{variable for vertical axis (response).}
  \item{z}{variable for out-of-screen axis.}
  \item{xlab, ylab, zlab}{axis labels.}
  \item{axis.scales}{if \code{TRUE}, label the values of the ends of the axes.
    \emph{Note:} For \code{Identify3d} to work properly, the value of this argument must
    be the same as in \code{scatter3d}.}
  \item{axis.ticks}{if \code{TRUE}, print interior axis-``tick'' labels; the default is \code{FALSE}.
    (The code for this option was provided by David Winsemius.)}
  \item{revolutions}{number of full revolutions of the display.}
  \item{bg.col}{background colour; one of \code{"white"}, \code{"black"}.}
  \item{axis.col}{colours for axes; if \code{axis.scales} is \code{FALSE}, then
    the second colour is used for all three axes.}
  \item{surface.col}{vector of colours for regression planes,
    used in the order specified by \code{fit}; for multi-group plots, the colours are used for the regression
    surfaces and points in the several groups.}
  \item{surface.alpha}{transparency of regression surfaces, from \code{0.0} (fully transparent) to \code{1.0} (opaque); default is \code{0.5}.}
  \item{neg.res.col, pos.res.col}{colours for lines representing negative
    and positive residuals.}
  \item{square.col}{colour to use to plot squared residuals.}
  \item{point.col}{colour of points.}
  \item{text.col}{colour of axis labels.}
  \item{grid.col}{colour of grid lines on the regression surface(s).}
  \item{fogtype}{type of fog effect; one of \code{"exp2"}, \code{"linear"},
    \code{"exp"}, \code{"none".}}
  \item{residuals}{plot residuals if \code{TRUE}; if  \code{residuals="squares"},
    then the squared residuals are shown as squares (using code adapted from Richard
    Heiberger). Residuals are available only when there is one surface plotted.}
  \item{surface}{plot surface(s) (\code{TRUE} or \code{FALSE}).}
  \item{fill}{fill the plotted surface(s) with colour (\code{TRUE} or \code{FALSE}).}
  \item{grid}{plot grid lines on the regression surface(s) (\code{TRUE} or \code{FALSE}).}
  \item{grid.lines}{number of lines (default, 26) forming the grid, in each of the x and z directions.}
  \item{df.smooth}{degrees of freedom for the two-dimensional smooth regression surface;
    if \code{NULL} (the default), the \code{\link[mgcv]{gam}} function will select the degrees of freedom
    for a smoothing spline by generalized cross-validation; if a positive number, a fixed
    regression spline will be fit with the specified degrees of freedom.}
  \item{df.additive}{degrees of freedom for each explanatory variable in an additive regression;
    if \code{NULL} (the default), the \code{gam} function will select degrees of freedom
    for the smoothing splines by generalized cross-validation; if a positive number
    or a vector of two positive numbers, fixed
    regression splines will be fit with the specified degrees of freedom for each term.}
  \item{sphere.size}{general size of spheres representing points;
    the actual size is dependent on the number of observations.}
  \item{radius}{relative radii of the spheres representing the points. This is normally a vector of the
    same length as the variables giving the coordinates of the points, and for the \code{formula} method,
    that must be the case or the argument may be omitted, in which case spheres are the same size;
    for the \code{default} method, the default for the argument, \code{1}, produces spheres
    all of the same size. The radii are scaled so that their median is 1.}
  \item{threshold}{if the actual size of the spheres is less than the threshold,
    points are plotted instead.}
  \item{speed}{relative speed of revolution of the plot.}
  \item{fov}{field of view (in degrees); controls degree of perspective.}
  \item{fit}{one or more of \code{"linear"}, \code{"quadratic"}, \code{"smooth"},
    \code{"additive"}; to display fitted surface(s); partial matching is supported --
    e.g., \code{c("lin", "quad")}.}
  \item{groups}{if \code{NULL} (the default), no groups are defined; if a factor, a different surface
    or set of surfaces is plotted for each level of the factor; in this event, the colours in
    \code{surface.col} are used successively for the points, surfaces, and residuals corresponding to
    each level of the factor.}
  \item{parallel}{when plotting surfaces by \code{groups}, should the surfaces be constrained to be
    parallel? A logical value, with default \code{TRUE}.}
  \item{ellipsoid}{plot concentration ellipsoid(s) (\code{TRUE} or \code{FALSE}).}
  \item{level}{expected proportion of bivariate-normal observations included in the
    concentration ellipsoid(s); default is 0.5.}
  \item{ellipsoid.alpha}{transparency of ellipsoids, from \code{0.0} (fully transparent) to \code{1.0} (opaque); default is \code{0.1}.}
  \item{id}{\code{FALSE}, \code{TRUE}, or a list controlling point identification, 
    similar to \code{\link{showLabels}} for 2D plots (see Details).}
  \item{model.summary}{print summary or summaries of the model(s) fit
    (\code{TRUE} or \code{FALSE}). \code{scatter3d} rescales the three variables
    internally to fit in the unit cube; this rescaling will affect regression
    coefficients.}
  \item{labels}{text labels for the points, one for each point;  defaults to
    the observation indices.}
  \item{col}{colours for the point labels, given by group. There must be at
    least as many colours as groups; if there are no groups, the first colour is used. Normally, the colours
    would correspond to the \code{surface.col} argument to \code{scatter3d}.}
  \item{offset}{vertical displacement for point labels (to avoid overplotting the points).}
  \item{\dots}{arguments to be passed down.}
}

\value{
  \code{scatter3d} does not return a useful value; it is used for its side-effect of
  creating a 3D scatterplot. \code{Identify3d} returns the labels of the
  identified points.
}

\details{
The \code{id} argument to \code{scatter3d} can be \code{FALSE}, \code{TRUE} (in which case 2
points will be identified according to their Mahalanobis distances from the center of the data),
or a list containing any or all of the following elements:
\describe{
   \item{method}{if \code{"mahal"} (the default), relatively extreme points are identified automatically
  	according to their Mahalanobis distances from the centroid (point of means);
  	if \code{"identify"}, points are identified interactively by right-clicking and dragging a box around them; right-click
  	in an empty area to exit from interactive-point-identification mode; if \code{"xz"}, identify extreme points in the
  	predictor plane; if \code{"y"}, identify unusual values of the response; if \code{"xyz"} identify
  	unusual values of an variable; if \code{"none"}, no point identification.
  	See \code{\link{showLabels}} for more information.}
  \item{n}{Number of relatively extreme points to identify automatically (default, \code{2},
    unless \code{method="identify"}, in which case identification continues until the user
    exits).}
  \item{labels}{text labels for the points, one for each point; in the \code{default} method defaults to
    the observation indices, in the \code{formula} method to the row names of the data.}
  \item{offset}{vertical displacement for point labels (to avoid overplotting the points).}
}
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\seealso{\code{\link[rgl]{rgl-package}}, \code{\link[mgcv]{gam}}}

\examples{
    if(interactive() && require(rgl) && require(mgcv)){
scatter3d(prestige ~ income + education, data=Duncan, id=list(n=3))
Sys.sleep(5) # wait 5 seconds
scatter3d(prestige ~ income + education | type, data=Duncan)
Sys.sleep(5)
scatter3d(prestige ~ income + education | type, surface=FALSE,
	ellipsoid=TRUE, revolutions=3, data=Duncan)
scatter3d(prestige ~ income + education, fit=c("linear", "additive"),
	data=Prestige)
Sys.sleep(5)
scatter3d(prestige ~ income + education | type,
    radius=(1 + women)^(1/3), data=Prestige)
	}
	\dontrun{
# drag right mouse button to identify points, click right button in open area to exit
scatter3d(prestige ~ income + education, data=Duncan, id=list(method="identify"))
scatter3d(prestige ~ income + education | type, data=Duncan, id=list(method="identify"))
    }
}

\keyword{hplot}
