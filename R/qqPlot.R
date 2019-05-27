# Quantile-comparison plots (J. Fox)

# last modified 30 September 2009 by J. Fox
# November 2009 by S. Weisberg -- changed to use showLabels for point identification
# 14 April 2010: set id.n = 0. J. Fox
# 1 June 2010: set reps=100 in qqPlot.lm. J. Fox
# 28 June 2010: fixed labeling bug S. Weisberg
# 11 March 2011: moved up ... argument. J. Fox
# 23 May 2012: line="none" now honored in qqPlot.default. J. Fox
# 2 May 2013: qqPlot.lm() now works with "aov" objects (fixing problem reported by Thomas Burk). J. Fox
# 2015-12-12: allow vectorized col, pch, and cex arguments (suggestion of Emmanuel Curis)
# 2017-02-12: consolidated id argument. J. Fox
# 2017-02-16: replace rlm() with MASS::rlm(). J. Fox
# 2017-03-25: fixed handling of indexes so that unsorted indexes are reported. J. Fox
# 2017-06-27: added formula method and plotting by groups. J. Fox
# 2017-10-26: fix qqPlot.lm() so that it doesn't return names identifical to indices. J. Fox
# 2017-11-30: substitute carPalette() for palette(). J. Fox
# 2018-03-23: properly return point IDs when method="identify"
# 2018-11-04: fixed handling of NAs when plotting by groups in qqPlot.default(). J. Fox
# 2019-04-09: respect order of factor levels when plotting by groups (problem and fix by Vilmantas Gegzna). J. Fox

qqp <- function(...) qqPlot(...)

qqPlot<-function(x, ...) {
	UseMethod("qqPlot")
}

qqPlot.default <- function(x, distribution="norm", groups, layout,
        ylim=range(x, na.rm=TRUE), ylab=deparse(substitute(x)),
		xlab=paste(distribution, "quantiles"), glab=deparse(substitute(groups)),
		main=NULL, las=par("las"),
		envelope=.95,
		col=carPalette()[1], col.lines=carPalette()[2], lwd=2, pch=1, cex=par("cex"),
		line=c("quartiles", "robust", "none"), id=TRUE, grid=TRUE, ...){
    if (!missing(groups)){
        if (isTRUE(id)) id <- list(n=2)
        if (is.null(id$labels)) id$labels <- seq(along=x)
        grps <- levels(as.factor(groups))
        if (missing(layout)) layout <- mfrow(length(grps), max.plots=12)
        if (prod(layout) < length(grps)) stop("layout cannot accomodate ", length(grps), " plots")
        oldpar <- par(mfrow=layout)
        on.exit(par(oldpar))
        for (group in grps){
            id.gr <- id
            if (!isFALSE(id)) id.gr$labels <- id$labels[groups == group]
            qqPlot.default(x[groups == group], distribution=distribution, ylim=ylim, ylab=ylab,
                           xlab=xlab, main=paste(glab, "=", group), las=las, envelope=envelope, col=col,
                           col.lines=col.lines, pch=pch, cex=cex, line=line, id=id.gr, grid=grid, ...)
        }
        return(invisible(NULL))
    }
    id <- applyDefaults(id, defaults=list(method="y", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
    if (isFALSE(id)){
        id.n <- 0
        id.method <- "none"
        labels <- id.cex <- id.col <- id.location <- NULL
    }
    else{
        labels <- id$labels
        if (is.null(labels)) labels <- if(!is.null(names(x))) names(x) else seq(along=x)
        id.method <- id$method
        id.n <- if ("identify" %in% id.method) Inf else id$n
        id.cex <- id$cex
        id.col <- id$col
        id.location <- id$location
    }
	line <- match.arg(line)
	index <- seq(along=x)
	good <- !is.na(x)
	ord <- order(x[good])
	if (length(col) == length(x)) col <- col[good][ord]
	if (length(pch) == length(x)) pch <- pch[good][ord]
	if (length(cex) == length(x)) cex <- cex[good][ord]
	ord.x <- x[good][ord]
	ord.lab <- labels[good][ord]
	q.function <- eval(parse(text=paste("q", distribution, sep="")))
	d.function <- eval(parse(text=paste("d", distribution, sep="")))
	n <- length(ord.x)
	P <- ppoints(n)
	z <- q.function(P, ...)
	plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=main, las=las, ylim=ylim)
	if(grid){
		grid(lty=1, equilogs=FALSE)
		box()}
	points(z, ord.x, col=col, pch=pch, cex=cex)
	if (line == "quartiles" || line == "none"){
		Q.x <- quantile(ord.x, c(.25,.75))
		Q.z <- q.function(c(.25,.75), ...)
		b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
		a <- Q.x[1] - b*Q.z[1]
		if (line == "quartiles") abline(a, b, col=col.lines, lwd=lwd)
	}
	if (line=="robust") {
		coef <- coef(MASS::rlm(ord.x ~ z))
		a <- coef[1]
		b <- coef[2]
		abline(a, b, col=col.lines, lwd=lwd)
	}
	conf <- if (envelope == FALSE) .95 else envelope
	zz <- qnorm(1 - (1 - conf)/2)
	SE <- (b/d.function(z, ...))*sqrt(P*(1 - P)/n)
	fit.value <- a + b*z
	upper <- fit.value + zz*SE
	lower <- fit.value - zz*SE
	if (envelope != FALSE) {
		lines(z, upper, lty=2, lwd=lwd, col=col.lines)
		lines(z, lower, lty=2, lwd=lwd, col=col.lines)
	}
	extreme <- showLabels(z, ord.x, labels=ord.lab,
			method = id.method, n = id.n, cex=id.cex, col=id.col, location=id.location)
	if (is.numeric(extreme)){
  	nms <- names(extreme)
  	extreme <- index[good][ord][extreme]
  	if (!all(as.character(extreme) == nms)) names(extreme) <- nms
	}
	if (length(extreme) > 0) extreme else invisible(NULL)
}

qqPlot.formula <- function (formula, data, subset, id=TRUE, ylab, glab, ...) {
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent()))))
        m$data <- as.data.frame(data)
    m$formula <- m$... <- m$id <- m$ylab <- NULL
    m[[1]] <- as.name("model.frame")
    if (missing(ylab)) ylab <- as.character(formula[[2]])
    if (length(formula) == 2){
        groups <- FALSE
    }
    else if (length(formula) == 3){
        groups <- TRUE
        if(missing(glab)) glab <- as.character(formula[[3]])
    }
    m$formula <-formula
    if (missing(data)){
        x <- na.omit(eval(m, parent.frame()))
    }
    else{
        x <- eval(m, parent.frame())
    }
    if (!isFALSE(id)){
        if (isTRUE(id)){
            id <- list(n=2)
        }
        if (is.null(id$labels)){
            id$labels <- rownames(x)
        }
    }
    if (!groups && ncol(x) > 1) stop("more than one variable specified")
    if (groups && ncol(x) != 2) stop("formula must be of the form variable ~ groups")
    if (!groups) qqPlot(x[, 1], id=id, ylab=ylab, ...)
    else qqPlot(x[, 1], groups=x[, 2], id=id, ylab=ylab, glab=glab, ...)
}


qqPlot.lm <- function(x, xlab=paste(distribution, "Quantiles"),
		ylab=paste("Studentized Residuals(", deparse(substitute(x)), ")", sep=""), main=NULL,
		distribution=c("t", "norm"), line=c("robust", "quartiles", "none"), las=par("las"),
		simulate=TRUE, envelope=.95,  reps=100,
		col=carPalette()[1], col.lines=carPalette()[2], lwd=2, pch=1, cex=par("cex"),
		id=TRUE, grid=TRUE, ...){
    distribution <- match.arg(distribution)
    force(xlab)
    force(ylab)
    x <- update(x, na.action="na.exclude")
    id <- applyDefaults(id, defaults=list(method="y", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
    if (isFALSE(id)){
        id.n <- 0
        id.method <- "none"
        labels <- id.cex <- id.col <- id.location <- NULL
    }
    else{
        labels <- id$labels
        if (is.null(labels)) labels <- names(residuals(x))
        if (length(labels) != length(residuals(x))) 
            stop("the number of labels, ", length(labels), ", differs from the number of cases, ", length(residuals(x)))
        id.method <- id$method
        id.n <- if ("identify" %in% id.method) Inf else id$n
        id.cex <- id$cex
        id.col <- id$col
        id.location <- id$location
    }
	result <- NULL
	line <- match.arg(line)
	rstudent <- rstudent(x)
	index <- seq(along=rstudent)
	sumry <- summary.lm(x)
	res.df <- sumry$df[2]
	if(!simulate){
		result <- qqPlot(rstudent, distribution=if (distribution == "t") "t" else "norm", df=res.df-1, line=line,
				main=main, xlab=xlab, ylab=ylab, las=las, envelope=envelope,
				col=col, col.lines=col.lines, lwd=lwd, pch=pch, cex=cex,
				id=list(method=id.method, n=id.n, cex=id.cex,
				col=id.col, location="lr", labels=labels), ...)
	}
	else {
	    good <- !is.na(rstudent)
	    rstudent <- rstudent[good]
	    labels <- labels[good]
		n <- length(rstudent)
		ord <- order(rstudent)
		ord.x <- rstudent[ord]
		ord.lab <- labels[ord]
		P <- ppoints(n)
		z <- if (distribution == 't') qt(P, df=res.df-1) else qnorm(P)
		plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=main, las=las)
		if(grid) grid(lty=1, equilogs=FALSE)
		points(z, ord.x, pch=pch, col=col, cex=cex)
		yhat <- na.omit(fitted.values(x))
		S <- sumry$sigma
		Y <- matrix(yhat, n, reps) + matrix(rnorm(n*reps, sd=S), n, reps)
		X <- model.matrix(x)
		rstud <- apply(rstudent(lm(Y ~ X - 1)), 2, sort)
		lower <- apply(rstud, 1, quantile, prob=(1 - envelope)/2)
		upper <- apply(rstud, 1, quantile, prob=(1 + envelope)/2)
		lines(z, upper, lty=2, lwd=lwd, col=col.lines)
		lines(z, lower, lty=2, lwd=lwd, col=col.lines)
		if (line == "quartiles"){
			Q.x <- quantile(rstudent, c(.25,.75))
			Q.z <- if (distribution == 't') qt(c(.25,.75), df=res.df - 1) else qnorm(c(.25,.75))
			b <- (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
			a <- Q.x[1] - b*Q.z[1]
			abline(a, b, col=col.lines, lwd=lwd)
		}
		if (line=="robust"){
			coef <- coefficients(MASS::rlm(ord.x~z))
			a <- coef[1]
			b <- coef[2]
			abline(a, b, col=col.lines, lwd=lwd)
		}
		result <- showLabels(z, ord.x,labels=ord.lab,
				method = id.method, n = id.n, cex=id.cex, col=id.col, location=id.location)
		nms <- names(result)
		result <- index[good][ord][result]
		names(result) <- nms
		result
	}
	if (all(as.character(result) == names(result))) names(result) <- NULL
	if (length(result) == 0) invisible(result) else if (is.numeric(result)) sort(result) else result
}

qqPlot.glm <- function(x, ...){
	stop("QQ plot for studentized residuals not available for glm")
}
