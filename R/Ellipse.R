# Ellipses (J. Fox and G. Monette)

# last modified 24 January 2010 by J. Fox
# added grid lines, 25 May 2010 by S. Weisberg

ellipse <- function(center, shape, radius, log="", center.pch=19, center.cex=1.5, segments=51, add=TRUE, 
	xlab="", ylab="", las=par("las"), col=palette()[2], lwd=2, lty=1, 
  grid=TRUE, ...) {
	logged <- function(axis=c("x", "y")){
		axis <- match.arg(axis)
		0 != length(grep(axis, log))
	}
	if (! (is.vector(center) && 2==length(center))) stop("center must be a vector of length 2")
	if (! (is.matrix(shape) && all(2==dim(shape)))) stop("shape must be a 2 by 2 matrix")
	angles <- (0:segments)*2*pi/segments 
	unit.circle <- cbind(cos(angles), sin(angles)) 
	ellipse <- t(center + radius*t(unit.circle %*% chol(shape))) 
	colnames(ellipse) <- c("x", "y")
	if (logged("x")) ellipse[, "x"] <- exp(ellipse[, "x"])
	if (logged("y")) ellipse[, "y"] <- exp(ellipse[, "y"])
	if (add) lines(ellipse, col=col, lwd=lwd, lty=lty, ...) 
	else {plot(ellipse, type="n", xlab = xlab, ylab = ylab, las=las ) 
        if(grid){
             grid(lty=1, equilogs=FALSE)
             box()}
        lines(ellipse, col=col, lwd=lwd, lty=lty, ... )} 	
	if (center.pch) points(center[1], center[2], pch=center.pch, cex=center.cex, col=col)
  }

dataEllipse <- function(x, y, log="", levels=c(0.5, 0.95), center.pch=19, 
  center.cex=1.5,
	plot.points=TRUE, add=!plot.points, segments=51, robust=FALSE, 
	xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), 
	las=par("las"), col=palette()[2], pch=1, lwd=2, lty=1, grid=TRUE, ...) {
	if (length(col) == 1) col <- rep(col, 2)
	if(missing(y)){
		if (is.matrix(x) && ncol(x) == 2) {
			if (missing(xlab)) xlab <- colnames(x)[1]
			if (missing(ylab)) ylab <- colnames(x)[2]
			y <- x[,2]
			x <- x[,1]
		}
		else stop("x and y must be vectors, or x must be a 2 column matrix")
	}
	else if(!(is.vector(x) && is.vector(y) && length(x) == length(y)))
		stop("x and y must be vectors of the same length")
	if (plot.points && !add) {
      plot(x, y, type="n", xlab=xlab, ylab=ylab, las=las, ...) 
	    if(grid){
        grid(lty=1, equilogs=FALSE)
        box()}
      points(x, y, col=col[1], pch=pch, ...)}
	if (plot.points && add)  points(x, y, col=col[1], pch=pch, ...)
	
	dfn <- 2
	dfd <- length(x) - 1
	if (robust) {
		v <- cov.trob(cbind(x, y))
		shape <- v$cov
		center <- v$center
	}
	else {
		shape <- var(cbind(x, y))
		center <- c(mean(x), mean(y))
	}
	for (level in levels) {
		radius <- sqrt(dfn * qf(level, dfn, dfd ))
		ellipse(center, shape, radius, log=log,
			center.pch=center.pch, center.cex=center.cex, segments=segments, 
			col=col[2], lty=lty, lwd=lwd, ...)
	}
}

confidenceEllipse <- function (model, ...) {
	UseMethod("confidenceEllipse")
}

confidenceEllipse.lm <- function(model, which.coef, levels=0.95, Scheffe=FALSE, 
	center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
	las=par("las"), col=palette()[2], lwd=2, lty=1,  ...){
	which.coef <- if(length(coefficients(model)) == 2) c(1, 2)
		else{
			if (missing(which.coef)){
				if (has.intercept(model)) c(2,3) else c(1, 2)
			} else which.coef
		}
	coef <- coefficients(model)[which.coef]
	xlab <- if (missing(xlab)) paste(names(coef)[1], "coefficient")
	ylab <- if (missing(ylab)) paste(names(coef)[2], "coefficient")
	dfn <- if (Scheffe) sum(df.terms(model)) else 2
	dfd <- df.residual(model)
	shape <- vcov(model)[which.coef, which.coef]
	for (level in rev(sort(levels))){
		radius <- sqrt(dfn*qf(level, dfn, dfd))
		add<-!level == max(levels)
		ellipse(coef, shape, radius, add=add, xlab=xlab, ylab=ylab,
			center.pch=center.pch, center.cex=center.cex, segments=segments, 
			col=col, lwd=lwd, lty=lty, las=las, ...)
	}
}


confidenceEllipse.glm <- function(model, which.coef, levels=0.95, Scheffe=FALSE, 
	center.pch=19, center.cex=1.5, segments=51, xlab, ylab,
	las=par("las"), col=palette()[2], lwd=2, lty=1, ...){
	which.coef <- if(length(coefficients(model)) == 2) c(1, 2)
		else{
			if (missing(which.coef)){
				if (has.intercept(model)) c(2, 3) else c(1, 2)
			} else which.coef
		}
	coef <- coefficients(model)[which.coef]
	xlab <- if (missing(xlab)) paste(names(coef)[1], "coefficient")
	ylab <- if (missing(ylab)) paste(names(coef)[2], "coefficient")
	df <- if (Scheffe) sum(df.terms(model)) else 2
	sumry <- summary(model, corr = FALSE)
	shape <- vcov(model)[which.coef, which.coef]
	for (level in rev(sort(levels))){
		radius <- sqrt(qchisq(level, df))
		add <- !level==max(levels)
		ellipse(coef, shape, radius, add=add, xlab=xlab, ylab=ylab,
			center.pch=center.pch, center.cex=center.cex, segments=segments,
			col=col, lwd=lwd, lty=lty, las=las, ...)
	}
}
