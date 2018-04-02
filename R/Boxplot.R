# checked in 26 December 2009 by J. Fox
# 2012-12-12: Fixed Boxplot.default() so that it works properly when g is numeric. J. Fox
# 2013-04-10: handles at argument properly, contribution of Steve Ellison. J. Fox
# 2013-08-19: removed loading of stats package. J. Fox
# 2016-09-30: added list, data.frame, and matrix methods, suggestion of Michael Friendly. J. Fox
# 2016-10-01: tweaked data.frame and list methods. J. Fox
# 2017-01-11: consolidate id argument
# 2017-10-03: add col and cex to id argument

Boxplot <- function(y, ...){
  arg.list <- list(...)
  if (!is.null(arg.list$horizontal) && isTRUE(arg.list$horizontal))
    stop("Boxplot does not support horizontal=TRUE")
	UseMethod("Boxplot")
}

Boxplot.default <- function(y, g, id=TRUE, xlab, ylab, ...){
    if (isFALSE(id)) {
        id.method="none"
        labels <- NULL
        id.n <- 0
        id.cex <- NULL
        id.col <- NULL
    }
    else{
        id <- applyDefaults(id, defaults=list(method="y", n=10, location="lr", cex=1,  col=carPalette()[1]), type="id")
        id.method <- match.arg(id$method, c("y", "identify", "none"))
        id.n <- id$n
        id.location <- id$location
        labels <- if (is.null(id$labels)) seq(along = y)
            else id$labels
        id.cex <- id$cex
        id.col <- id$col
    }
    if (missing(ylab)) 
        ylab <- deparse(substitute(y))
    pars <- list(...)
    if (missing(g)) {
        valid <- complete.cases(y, labels)
        y <- y[valid]
        labels <- labels[valid]
        b <- boxplot(y, ylab = ylab, ...)
        if (id.method == "none" | id.n == 0) 
            return(invisible(NULL))
        else if (id.method == "identify") {
            res <- identify(rep(1, length(y)), y, labels, cex=id.cex, col=id.col)
            return(if (length(res) == 0) invisible(NULL) else labels[res])
        }
        else if (length(b$out) > 0) {
            sel <- y %in% b$out
            yy <- y[sel]
            labs <- labels[sel]
            which.low <- yy < b$stats[1, 1]
            y.low <- yy[which.low]
            labs.low <- labs[which.low]
            if (length(y.low) > id.n) {
                ord.low <- order(y.low)[1:id.n]
                y.low <- y.low[ord.low]
                labs.low <- labs.low[ord.low]
            }
            which.high <- yy > b$stats[5, 1]
            y.high <- yy[which.high]
            labs.high <- labs[which.high]
            if (length(y.high) > id.n) {
                ord.high <- order(y.high, decreasing = TRUE)[1:id.n]
                y.high <- y.high[ord.high]
                labs.high <- labs.high[ord.high]
            }
            labs <- c(labs.low, labs.high)
            at <- if(!is.null(pars$at)) pars$at else 1		
            if (id.location == "lr") text(at, c(y.low, y.high), labs, pos = 2, xpd=TRUE, cex=id.cex, col=id.col)
            else maptools::pointLabel(c(at, at), c(y.low, y.high, y.low, y.high), 
                                      c(paste0(" ", labs, " "), rep(" ", length(labs))), 
                                      xpd=TRUE, col=id.col, cex=id.cex)
            return(if (length(labs) == 0) invisible(NULL) else labs)
        }
        else return(invisible(NULL))
    }
    else {
        if (missing(xlab)) 
            xlab = deparse(substitute(g))
        valid <- complete.cases(y, labels, g)
        y <- y[valid]
        labels <- labels[valid]
        g <- g[valid]
        b <- boxplot(split(y, g), ylab = ylab, xlab = xlab, ...)
        levels <- if (is.factor(g)) 
            levels(g)
        else sort(unique(g))
        gg <- as.numeric(g)
        if (id.method == "none" | id.n == 0) 
            return(invisible(NULL))
        else if (id.method == "identify") {
            res <- identify(gg, y, labels)
            return(if (length(res) == 0) invisible(NULL) else labels[res])
        }
        else {
            midx <- mean(par("usr")[1:2])
            identified <- character(0)
            if (length(b$out) > 0) {
                groups <- unique(b$group)
                for (group in groups) {
                    grp <- g == levels[group]
                    yy <- y[grp]
                    labs <- labels[grp]
                    sel <- yy %in% b$out[b$group == group]
                    yy <- yy[sel]
                    glabs <- labs[sel]
                    which.low <- yy < b$stats[1, group]
                    y.low <- yy[which.low]
                    labs.low <- glabs[which.low]
                    if (length(y.low) > id.n) {
                        ord.low <- order(y.low)[1:id.n]
                        y.low <- y.low[ord.low]
                        labs.low <- labs.low[ord.low]
                    }
                    which.high <- yy > b$stats[5, group]
                    y.high <- yy[which.high]
                    labs.high <- glabs[which.high]
                    if (length(y.high) > id.n) {
                        ord.high <- order(y.high, decreasing = TRUE)[1:id.n]
                        y.high <- y.high[ord.high]
                        labs.high <- labs.high[ord.high]
                    }
                    pos <- if (group < midx) 
                        4
                    else 2
                    at <- if(!is.null(pars$at)) pars$at[group] else group
                    labs <- c(labs.low, labs.high)
                    if (id.location == "lr") text(at, c(y.low, y.high), labs, pos = pos, xpd=TRUE, col=id.col, cex=id.cex)
                    else maptools::pointLabel(c(at, at), c(y.low, y.high, y.low, y.high), 
                                                   c(paste0(" ", labs, " "), rep(" ", length(labs))), 
                                                   xpd=TRUE, col=id.col, cex=id.cex)
                    identified <- c(identified, c(labs.low, labs.high))
                }
            }
            return(if (length(identified) == 0) invisible(NULL) else identified)
        }
    }
}


Boxplot.formula <- function(formula, data=NULL, subset, na.action=NULL, id=TRUE, xlab, ylab, ...){
	# much of this function adapted from graphics:boxplot.formula
    if (isFALSE(id)) {
        id.method="none"
        id.n <- 0
        labels <- NULL
        id.location <- "y"
        id.col <- NULL
        id.cex <- NULL
    }
    else{
        id <- applyDefaults(id, defaults=list(method="y", n=10, location="lr", cex=1,  col=carPalette()[1]), type="id")
        id.method <- match.arg(id$method, c("y", "identify", "none"))
        id.n <- id$n
        id.location <- id$location
        labels <- id$labels
        id.col <- id$col
        id.cex <- id$cex
    }
	m <- match.call(expand.dots = FALSE)
	if (is.matrix(eval(m$data, parent.frame()))) 
		m$data <- as.data.frame(data)
	m$xlab <- m$ylab <- m$id <- m$... <- NULL
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")
	mf <- eval(m, parent.frame())
	mf$"(labels.)" <- if (is.null(labels)) rownames(mf) else labels
	lab.var <- which(names(mf) == "(labels.)")
	if (length(formula) == 3){
		response <- attr(attr(mf, "terms"), "response")
		if (missing(ylab)) ylab <- names(mf)[response]
		if (missing(xlab)) xlab <- names(mf)[-c(response, lab.var)]
        x <- mf[, -c(response, lab.var)]
		if (is.data.frame(x)) x <- do.call("interaction", as.list(x))
        if (length(xlab) > 1) xlab <- paste(xlab, collapse="*")
		Boxplot(mf[[response]], x, id=list(method=id.method, labels=mf[[lab.var]], 
		                                   n=id.n, location=id.location, col=id.col, cex=id.cex), 
            xlab=xlab, ylab=ylab, ...)
	}
	else if (length(formula) == 2){
		if (missing(ylab)) ylab <- names(mf)[-lab.var]
		Boxplot(mf[, -lab.var], id=list(method=id.method, 
		                                labels=mf[[lab.var]], n=id.n, location=id.location, col=id.col, cex=id.cex), ylab=ylab, ...)
	}
	else stop("improper Boxplot formula")   
}

Boxplot.list <- function(y, xlab="", ylab="", ...){
  if (is.null(names(y))) names(y) <- 1:length(y)
  g <- factor(rep(names(y), sapply(y, length)), levels=names(y))
  y <- do.call(c, y)
  Boxplot(y, g, xlab=xlab, ylab=ylab, ...)
}

Boxplot.data.frame <-  function(y, id=TRUE, ...){
    if (isFALSE(id)) {
        id.method="none"
        id.n <- 0
        labels <- NULL
        id.location <- "y"
        id.col <- NULL
        id.cex <- NULL
    }
    else{
        id <- applyDefaults(id, defaults=list(method="y", n=10, location="lr", labels=rownames(y), cex=1,  col=carPalette()[1]), type="id")
        id.method <- match.arg(id$method, c("y", "identify", "none"))
        id.n <- id$n
        id.location <- id$location
        labels <- rep(id$labels, ncol(y))
        id.col <- id$col
        id.cex <- id$cex
    }
  Boxplot(as.list(y), id=list(method=id.method, n=id.n, location=id.location, labels=labels, cex=id.cex, col=id.col), ...)
}

Boxplot.matrix <- function(y, ...){
  Boxplot(as.data.frame(y), ...)
}
