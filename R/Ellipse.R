# Ellipses (orignally by J. Fox and G. Monette)

# added grid lines, 25 May 2010 by S. Weisberg
# arguments more consistent with other functions; ... passes args to plot, 5 Sept 2010 by J. Fox
# confidenceEllipse.lm and .glm can add to current plot, applying patch from Rafael Laboissiere, 17 Oct 2010 by J. Fox
# added fill and fill.alpha arguments for translucent fills (suggested by Michael Friendly), 14 Nov 2010 by J. Fox
# modified 2 May 2011 by Michael Friendly
#   - allow pivot=TRUE (with warning)
#   - barf on non-symmetric shape
#   - return coordinates of ellipse invisibly
# dataEllipse() and confidenceEllipse() invisibly return coordinates,  3 May 2011 by J. Fox
# Modified 5 May 2011 by Michael Friendly
#   - dataEllipse now honors add=FALSE, plot.points=FALSE
# Modified 16 May 2011 by Michaell Friendly
#   - corrected bug introduced in dataEllipse via allowing pivot=TRUE 
# Modified 7 Aug 2011 by J. Fox: added draw argument
# Modified 28 Nov 2011  by J. Fox (suggested by Michael Friendly):
#   - corrected bug in xlab, ylab in confidenceEllipse()
#   - added dfn argument to .lm and .glm methods for confidenceEllipse()
# Modified 14&16 Dec 2011 by J. Fox (suggested by Michael Friendly) to add weights argument to dataEllipse().

ellipse <- function(center, shape, radius, log="", center.pch=19, center.cex=1.5, segments=51, draw=TRUE, add=draw, 
		xlab="", ylab="", col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3,
		grid=TRUE, ...) {
	trans.colors <- function(col, alpha=0.5, names=NULL) {
		# this function by Michael Friendly
		nc <- length(col)
		na <- length(alpha)
		# make lengths conform, filling out to the longest
		if (nc != na) {
			col <- rep(col, length.out=max(nc,na))
			alpha <- rep(alpha, length.out=max(nc,na))
		}
		clr <-rbind(col2rgb(col)/255, alpha=alpha)
		col <- rgb(clr[1,], clr[2,], clr[3,], clr[4,], names=names)
		col
	}
	logged <- function(axis=c("x", "y")){
		axis <- match.arg(axis)
		0 != length(grep(axis, log))
	}
	if (! (is.vector(center) && 2==length(center))) stop("center must be a vector of length 2")
	if (! (is.matrix(shape) && all(2==dim(shape)))) stop("shape must be a 2 by 2 matrix")
	if (max(abs(shape - t(shape)))/max(abs(shape)) > 1e-10) stop("shape must be a symmetric matrix")
	angles <- (0:segments)*2*pi/segments 
	unit.circle <- cbind(cos(angles), sin(angles)) 
#	ellipse <- t(center + radius*t(unit.circle %*% chol(shape,pivot=TRUE))) 
	Q <- chol(shape, pivot=TRUE)
	order <- order(attr(Q, "pivot"))
	ellipse <- t( center + radius*t( unit.circle %*% Q[,order]))
	colnames(ellipse) <- c("x", "y")
	if (logged("x")) ellipse[, "x"] <- exp(ellipse[, "x"])
	if (logged("y")) ellipse[, "y"] <- exp(ellipse[, "y"])
	fill.col <- trans.colors(col, fill.alpha)
	if (draw) {
		if (add) {
			lines(ellipse, col=col, lwd=lwd, ...) 
			if (fill) polygon(ellipse, col=fill.col, border=NA)
		}
		else {
			plot(ellipse, type="n", xlab = xlab, ylab = ylab, ...) 
			if(grid){
				grid(lty=1, equilogs=FALSE)
				box()}
			lines(ellipse, col=col, lwd=lwd, ... )
			if (fill) polygon(ellipse, col=fill.col, border=NA)
		} 	
		if (center.pch) points(center[1], center[2], pch=center.pch, cex=center.cex, col=col)
	}
	invisible(ellipse)
}

dataEllipse <- function(x, y, weights, log="", levels=c(0.5, 0.95), center.pch=19, 
		center.cex=1.5, draw=TRUE,
		plot.points=draw, add=!plot.points, segments=51, robust=FALSE, 
		xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), 
		col=palette()[1:2], lwd=2, fill=FALSE, fill.alpha=0.3, grid=TRUE, ...) {
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
	if (missing(weights)) weights <- rep(1, length(x))
	if (length(weights) != length(x)) stop("weights must be of the same length as x and y")
	if(draw) {
		if (!add) {
			plot(x, y, type="n", xlab=xlab, ylab=ylab,  ...) 
			if(grid){
				grid(lty=1, equilogs=FALSE)
				box()}
		}
		if (plot.points)  points(x, y, col=col[1], ...)
	}
	dfn <- 2
	dfd <- length(x) - 1
	if (robust) {
		use <- weights > 0
		v <- cov.trob(cbind(x[use], y[use]), wt=weights[use])
		shape <- v$cov
		center <- v$center
	}
	else {
		v <- cov.wt(cbind(x, y), wt=weights)
		shape <- v$cov
		center <- v$center
	}
	result <- vector("list", length=length(levels))
	names(result) <- levels
	for (i in seq(along=levels)) {
		level <- levels[i]
		radius <- sqrt(dfn * qf(level, dfn, dfd ))
		result[[i]] <- ellipse(center, shape, radius, log=log,
				center.pch=center.pch, center.cex=center.cex, segments=segments, 
				col=col[2], lwd=lwd, fill=fill, fill.alpha=fill.alpha, draw=draw, ...)
	}
	invisible(if (length(levels) == 1) result[[1]] else result)
}

confidenceEllipse <- function (model, ...) {
	UseMethod("confidenceEllipse")
}

confidenceEllipse.lm <- function(model, which.coef, levels=0.95, Scheffe=FALSE, dfn,
		center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
		col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, draw=TRUE, add=!draw, ...){
	which.coef <- if(length(coefficients(model)) == 2) c(1, 2)
			else{
				if (missing(which.coef)){
					if (has.intercept(model)) c(2,3) else c(1, 2)
				} else which.coef
			}
	coef <- coefficients(model)[which.coef]
	if (missing(xlab)) xlab <- paste(names(coef)[1], "coefficient")
	if (missing(ylab)) ylab <-  paste(names(coef)[2], "coefficient")
	if (missing(dfn)) dfn <- if (Scheffe) sum(df.terms(model)) else 2
	dfd <- df.residual(model)
	shape <- vcov(model)[which.coef, which.coef]
	levels <- rev(sort(levels))
	result <- vector("list", length=length(levels))
	names(result) <- levels
	for (i in seq(along=levels)){
		level <- levels[i]
		radius <- sqrt(dfn*qf(level, dfn, dfd))
		add.plot <- !level==max(levels) | add
		result[[i]] <- ellipse(coef, shape, radius, add=add.plot, xlab=xlab, ylab=ylab,
				center.pch=center.pch, center.cex=center.cex, segments=segments, 
				col=col, lwd=lwd, fill=fill, fill.alpha=fill.alpha, draw=draw, ...)
	}
	invisible(if (length(levels) == 1) result[[1]] else result)
}


confidenceEllipse.glm <- function(model, which.coef, levels=0.95, Scheffe=FALSE, dfn,
		center.pch=19, center.cex=1.5, segments=51, xlab, ylab,
		col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, draw=TRUE, add=!draw, ...){
	which.coef <- if(length(coefficients(model)) == 2) c(1, 2)
			else{
				if (missing(which.coef)){
					if (has.intercept(model)) c(2, 3) else c(1, 2)
				} else which.coef
			}
	coef <- coefficients(model)[which.coef]
	xlab <- if (missing(xlab)) paste(names(coef)[1], "coefficient")
	ylab <- if (missing(ylab)) paste(names(coef)[2], "coefficient")
	df <- if (!missing(dfn)) dfn
			else if (Scheffe) sum(df.terms(model)) else 2
	sumry <- summary(model, corr = FALSE)
	shape <- vcov(model)[which.coef, which.coef]
	levels <- rev(sort(levels))
	result <- vector("list", length=length(levels))
	names(result) <- levels
	for (i in seq(along=levels)){
		level <- levels[i]
		radius <- sqrt(qchisq(level, df))
		add.plot <- !level==max(levels) | add
		result[[i]] <- ellipse(coef, shape, radius, add=add.plot, xlab=xlab, ylab=ylab,
				center.pch=center.pch, center.cex=center.cex, segments=segments,
				col=col, lwd=lwd, fill=fill, fill.alpha=fill.alpha, draw=draw, ...)
	}
	invisible(if (length(levels) == 1) result[[1]] else result)
}
