# fancy scatterplot matrices (J. Fox)

# 2010-09-04: J. Fox: changed color choice
# 2010-09-16: fixed point color when col is length 1
# 2011-03-08: J. Fox: changed col argument

scatterplotMatrix <- function(x, ...){
	UseMethod("scatterplotMatrix")
}

scatterplotMatrix.formula <- function (x, data=NULL, subset, labels, ...) {
	m <- match.call(expand.dots = FALSE)
	if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
		m$data <- as.data.frame(data)
	m$labels <- m$formula <- m$... <- NULL
	m[[1]] <- as.name("model.frame")
	if (!inherits(x, "formula") | length(x) != 2) 
		stop("invalid formula")
	rhs <- x[[2]]
	if ("|" != deparse(rhs[[1]])){
		groups <- FALSE
	}
	else{
		groups <- TRUE
		x<-as.character(c(x))
		x<-as.formula(sub("\\|", "+", x))   
	}
	m$formula <-x
	if (missing(data)){ 
		X <- na.omit(eval(m, parent.frame()))
		if (missing(labels)) labels <- gsub("X", "", row.names(X))
	}
	else{
		if (!missing(labels)) row.names(data) <- labels
		X <- eval(m, parent.frame())
		labels <- row.names(X)
	}
	if (!groups) scatterplotMatrix(X, labels=labels, ...)
	else{
		ncol<-ncol(X)
		scatterplotMatrix.default(X[, -ncol], groups=X[, ncol], labels=labels, ...)
	}
}

scatterplotMatrix.default <- function(x, var.labels=colnames(x), 
	diagonal=c("density", "boxplot", "histogram", "oned", "qqplot", "none"), adjust=1, nclass,
	plot.points=TRUE, smooth=TRUE, spread=smooth && !by.groups, span=.5, loess.threshold=5, reg.line=lm, 
	transform=FALSE, family=c("bcPower", "yjPower"),
	ellipse=FALSE, levels=c(.5, .95), robust=TRUE,
	groups=NULL, by.groups=FALSE, 
	labels, id.method="mahal", id.n=0, id.cex=1, id.col=palette()[1],
	col=if (n.groups == 1) palette()[3:1] else rep(palette(), length=n.groups),
	pch=1:n.groups, lwd=1, lwd.smooth=lwd, lwd.spread=lwd, lty=1, lty.smooth=lty, lty.spread=2,
	cex=par("cex"), cex.axis=par("cex.axis"), cex.labels=NULL, 
	cex.main=par("cex.main"), 
	legend.plot=length(levels(groups)) > 1, row1attop=TRUE, ...){
	spread # force evaluation
	if (id.method == "identify") stop("interactive point identification not permitted")
	family <- match.arg(family)
	err <- ""
	lowess.line <- function(x, y, col, span) {
		warn <- options(warn=-1)
		valid <- complete.cases(x, y)
		x <- x[valid]
		y <- y[valid]
		ord <- order(x)
		x <- x[ord]
		y <- y[ord]
#		if (length(unique(x)) < lowess.threshold || length(unique(y)) < lowess.threshold) return()
		if (length(unique(y)) < loess.threshold) return()
		if (!spread){
			fit <- try(loess.smooth(x, y, span=span), silent=TRUE)
			if (class(fit) == "try-error"){
				err <<- c(err, "smooth")
				options(warn)
				return()
			}
			lines(fit$x, fit$y, lty=lty.smooth, lwd=lwd.smooth, col=col)
		}
		else{
			fit <- try(loess(y ~ x, degree=1, family="symmetric", span=span), silent=TRUE)
			if (class(fit) == "try-error"){
				err <<- c(err, "smooth")
				options(warn)
				return()
			}
			res <- residuals(fit)
			pos <- res > 0
			pos.fit <- try(loess(res^2 ~ x, span=span, degree=0, family="gaussian", subset=pos), silent=TRUE)
			neg.fit <- try(loess(res^2 ~ x, span=span, degree=0, family="gaussian", subset=!pos), silent=TRUE)
			if (class(pos.fit) == "try-error" || class(neg.fit) == "try-error"){
				err <<- c(err, "spread")
				options(warn)
				return()
			}
			lines(x, fitted(fit), lty=lty.smooth, lwd=lwd.smooth, col=col)
			y.pos <- fitted(fit)[pos] + sqrt(fitted(pos.fit))
			lines(x[pos], y.pos, lty=lty.spread, lwd=lwd.spread, col=col)
			y.neg <- fitted(fit)[!pos] - sqrt(fitted(neg.fit))
			lines(x[!pos], y.neg, lty=lty.spread, lwd=lwd.spread, col=col)
		}
		options(warn)
	}
	if (missing(labels)){
		labels <- rownames(x)
		if (is.null(labels)) labels <- as.character(seq(length.out=nrow(x)))
	}
	if (!(missing(groups))){
		x <- na.omit(data.frame(groups, labels, x, stringsAsFactors=FALSE))
#		groups <- as.factor(as.character(x[, 1]))
		if (!is.factor(groups)) groups <- as.factor(as.character(x[,1]))
		labels <- x[, 2]
		x <- x[, -(1:2)]
	}
	else {
		x <- na.omit(data.frame(labels, x, stringsAsFactors=FALSE))
		labels <- x[, 1]
		x <- x[, -1]
	}
	if (missing(nclass)) nclass <- "FD"
	reg <- function(x, y, col){
		mod<-reg.line(y ~ x)
		y.hat <- fitted.values(mod)
		x <- model.matrix(mod)[,2]
		min <- which.min(x)
		max <- which.max(x)
		lines(c(x[min], x[max]), c(y.hat[min], y.hat[max]), lty=lty, lwd=lwd, col=col)
	}
	legendPlot <- function(){
		usr <- par("usr")
		legend("bottomleft", bg="white",
			legend=levels(groups), pch=pch, col=col[1:n.groups],
			cex=cex)
	}	
	do.legend <- legend.plot	
# The following panel function adapted from Richard Heiberger
	panel.density <- function(x, ...){
		dens.x <- density(x, adjust = adjust)
		lines(dens.x$x, min(x) + dens.x$y * diff(range(x))/diff(range(dens.x$y)))
		rug(x)
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	panel.histogram <- function(x, ...){
		par(new=TRUE)
		hist(x, main="", axes=FALSE, breaks=nclass, col=col[1])
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	panel.boxplot <- function(x, ...){
		par(new=TRUE)
		boxplot(x, axes=FALSE, main="", col=col[1])
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	# The following panel function adapted from Richard Heiberger
	panel.oned <- function(x, ...) {
		range <- range(x)
		delta <- diff(range)/50
		y <- mean(range)
		segments(x - delta, x, x + delta, x, col = col[3])
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	panel.qqplot <- function(x, ...){
		par(new=TRUE)
		qqnorm(x, axes=FALSE, xlab="", ylab="", main="", col=col[3])
		qqline(x, col=col[1])
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	panel.blank <- function(x, ...){
		if (do.legend) legendPlot()
		do.legend <<- FALSE
	}
	which.fn <- match(match.arg(diagonal),
		c("density", "boxplot", "histogram", "oned", "qqplot", "none"))
	diag <- list(panel.density, panel.boxplot, panel.histogram, panel.oned, panel.qqplot, panel.blank)[[which.fn]]
	groups <- as.factor(if(missing(groups)) rep(1, length(x[, 1])) else groups)
	n.groups <- length(levels(groups))
	if (n.groups > length(col)) stop("number of groups exceeds number of available colors")
	if (length(col) == 1) col <- rep(col, 3)
	if (transform != FALSE | length(transform) == ncol(x)){
		if (transform == TRUE & length(transform) == 1){
			transform <- if (by.groups) coef(powerTransform(as.matrix(x) ~ groups, family=family), round=TRUE)
				else coef(powerTransform(x, family=family), round=TRUE)
		}
		for (i in 1:ncol(x)){
			x[, i] <- if (family == "bcPower") 
					bcPower(x[, i], transform[i])
				else yjPower(x[, i], transform[i])
			var.labels[i] <- paste(var.labels[i], "^(", round(transform[i],2), ")", sep="")
		}
	}
	labs <- labels
	pairs(x, labels=var.labels, 
		cex.axis=cex.axis, cex.main=cex.main, cex.labels=cex.labels, cex=cex,
		diag.panel=diag, row1attop = row1attop,
		panel=function(x, y, ...){ 
			for (i in 1:n.groups){
				subs <- groups == levels(groups)[i]
				if (plot.points) points(x[subs], y[subs], pch=pch[i], col=col[if (n.groups == 1) 3 else i], cex=cex)
				if (by.groups){
					if (smooth) lowess.line(x[subs], y[subs], col=col[i], span)
					if (is.function(reg.line)) reg(x[subs], y[subs], col=col[i])
					if (ellipse) dataEllipse(x[subs], y[subs], plot.points=FALSE, 
							levels=levels, col=col[i], robust=robust, lwd=1)
					showLabels(x[subs], y[subs], labs[subs], id.method=id.method, 
					    id.n=id.n, id.col=col[i], id.cex=id.cex)
					#if (id.method != "none") 
					#	showLabelsScatter(x[subs], y[subs], labs[subs], id.var=id.var, id.method=id.method,
					#		id.n=id.n, id.col=col[i], id.cex=id.cex)
				}
			}
			if (!by.groups){
				if (is.function(reg.line)) abline(reg.line(y ~ x), lty=lty, lwd=lwd, col=col[1])
				if (smooth) lowess.line(x, y, col=col[2], span)
				if (ellipse) dataEllipse(x, y, plot.points=FALSE, levels=levels, col=col[1],
						robust=robust, lwd=1)
				showLabels(x, y, labs, id.method=id.method, 
					    id.n=id.n, id.col=id.col, id.cex=id.cex)
				#if (id.method != "none") showLabelsScatter(x, y, labs, id.var=id.var, id.method=id.method, 
				#		id.n=id.n, id.col=col[1], id.cex=id.cex)
			}
		}, ...
	)
	if ("smooth" %in% err) warning("some smooths failed")
	if ("spread" %in% err) warning("some spreads failed")
}

spm <- function(x, ...){
	scatterplotMatrix(x, ...)
} 



