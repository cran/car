# Leverage plots (J. Fox)

# last modified 9 October 2009 by J. Fox
# modified 25 November for layout and marking points only
# changed 'vars' to 'terms' 16 March 2010 SW
# 14 April 2010: set id.n = 0. J. Fox
# 15 August 2010 S. Weisberg, added col.lines and col arguments
# 5 Sept 2010 J. Fox, pass ... down to plot() and points() etc.
# 16 June 2011 allow layout=NA, in which case the layout is not set in this
#  function, so it is the responsibility of the user
# 25 April 2016:  checks na.action for compatibility with Rcmdr.  SW
# 2017-02-12: consolidate id argument. J. Fox
# 2017-11-09: made consistent with vcov() in R 2.5.0. J Fox
# 2017-11-30: substitute carPalette() for palette(). J. Fox
# 2019-11-14: change class(x) == "y" to inherits(x, "y")

# these functions to be rewritten; simply renamed for now

leveragePlots <- function(model, terms= ~ ., layout=NULL, ask, 
		main, ...){
# Added for compatibility with Rcmdr
  if(inherits(model$na.action, "exclude")) model <- update(model, na.action=na.omit)
# End addition
	terms <- if(is.character(terms)) paste("~",terms) else terms
	vform <- update(formula(model),terms)
	terms.model <- attr(attr(model.frame(model), "terms"), "term.labels")
	terms.vform <- attr(terms(vform), "term.labels")
	good <- terms.model[match(terms.vform, terms.model)]
	nt <- length(good)
	if (nt == 0) stop("No plots specified")
	if (missing(main)) main <- if (nt == 1) "Leverage Plot" else "Leverage Plots"
  nr <- 0  
  if (nt > 1 & (is.null(layout) || is.numeric(layout))) {
    if(is.null(layout)){
         layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2), 
                             c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))
    }
    ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
            oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
    on.exit(par(op))
    }
	for (term in good) leveragePlot(model, term, main="", ...)
	mtext(side=3,outer=TRUE,main, cex=1.2)
	invisible(0)
}

leveragePlot <- function (model, ...) {
	UseMethod("leveragePlot")
}

leveragePlot.lm <- function(model, term.name,
		id=TRUE, col=carPalette()[1], col.lines=carPalette()[2], lwd=2, 
		xlab, ylab, main="Leverage Plot", grid=TRUE, ...){
    id <- applyDefaults(id, defaults=list(method=list(abs(residuals(model, type="pearson")), "x"), n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
    if (isFALSE(id)){
        id.n <- 0
        id.method <- "none"
        labels <- id.cex <- id.col <- id.location <- NULL
    }
    else{
        labels <- id$labels
        if (is.null(labels)) labels <- names(na.omit(residuals(model)))
        id.method <- id$method
        id.n <- if ("identify" %in% id.method) Inf else id$n
        id.cex <- id$cex
        id.col <- id$col
        id.location <- id$location
    }
	term.name <- if (is.character(term.name) & 1==length(term.name)) term.name
			else deparse(substitute(term.name))
	b <- coefficients(model)
	e <- na.omit(residuals(model))
	p <- length(b)
	I.p <- diag(p)
	term.names <- term.names(model)
	term <- which(term.name==term.names)
	if (0==length(term)) stop(paste(term.name,"is not a term in the model."))
	responseName <- responseName(model)
	intercept <- has.intercept(model)
	assign <- model$assign
	X <- model.matrix(model)
	V <- vcov(model, complete=FALSE)
	wt <- if (is.null(weights(model))) rep(1, length(X[,1]))
			else weights(model)
	subs <- which(assign==term-intercept)
	hypothesis.matrix <- I.p[subs,]
	L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
			else hypothesis.matrix
	u <- solve(L %*% V %*% t(L)) %*% L %*% b
	v.x <- X %*% V %*% t(L) %*% u
	v.y <- v.x + e
	if (missing(xlab)) xlab <- paste(term.names[term],"| others")
	if (missing(ylab)) ylab <- paste(responseName," | others")
	plot(v.x, v.y, xlab=xlab, ylab=ylab, type="n", ...)
	if(grid){
		grid(lty=1, equilogs=FALSE)
		box()}
	points(v.x, v.y, col=col, ...)	
	abline(lsfit(v.x, v.y, wt=wt), col=col.lines, lwd=lwd)
	showLabels(v.x, v.y, labels=labels, 
			method=id.method, n=id.n, cex=id.cex, 
			col=id.col, location=id.location)
}

leveragePlot.glm <- function(model, ...){
	stop("Leverage plot requires an lm object")
}

