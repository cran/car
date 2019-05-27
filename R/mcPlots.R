# October 1, 2014 mcPlots, by S. Weisberg and J. Fox
# 'mc' stands for Marginal and Conditional:  for the specified regressor X in model
# The 'marginal' plot is of Y vs X with Y and X both centered
# The 'conditional plot is the added-variable plot e(Y|rest) vs e(X|rest)
# If 'overlaid=TRUE', the default, both plots are overlayed
# If 'overlaid=FALSE', then the plots are side-by-side
# The 'overlaid' plot is similar to the initial and final frame of an ARES plot
# Cook and Weisberg (1989), "Regression diagnostics with dynamic graphics", Technometrics, 31, 277.
# This plot would benefit from animation.
# 2017-02-13: consolidated id and ellipse arguments. J. Fox
# 2017-11-30: substitute carPalette() for palette(). J. Fox
# 2018-12-17: added title argument, if title=FALSE, suppress unchangable main= arguments

mcPlots <- function(model, terms=~., layout=NULL, ask, overlaid=TRUE, ...){
    terms <- if(is.character(terms)) paste("~",terms) else terms
    vform <- update(formula(model),terms)
    if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
        stop("Only predictors in the formula can be plotted.")
    terms.model <- attr(attr(model.frame(model), "terms"), "term.labels")
    terms.vform <- attr(terms(vform), "term.labels")
    terms.used <- match(terms.vform, terms.model)
    mm <- model.matrix(model)
    model.names <- attributes(mm)$dimnames[[2]]
    model.assign <- attributes(mm)$assign
    good <- model.names[!is.na(match(model.assign, terms.used))]
    #	if (intercept) good <- c("(Intercept)", good)
    if(attr(attr(model.frame(model), "terms"), "intercept") == 0)
        stop("Error---the 'lm' object must have an intercept")
    nt <- length(good)
    if (nt == 0) stop("No plots specified")
    if(overlaid){
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
    } else{
        if (nt >= 1 & (is.null(layout) || is.numeric(layout))) {
            if(is.null(layout)){
                layout <- switch(min(nt, 4), c(1, 2), c(2, 2), c(3, 2), c(4, 2))
            }
            ask <- if(missing(ask) || is.null(ask)) layout[1] < nt else ask
            op <- par(mfrow=layout, ask=ask, no.readonly=TRUE,
                      oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
            on.exit(par(op))
        }
    }
    for (term in good) mcPlot(model, term, new=FALSE, overlaid=overlaid, ...)
    #	mtext(side=3,outer=TRUE,main, cex=1.2)
}



mcPlot <-  function(model, ...) UseMethod("mcPlot")

mcPlot.lm <- function(model, variable, id=FALSE,
                      col.marginal=carPalette()[2], col.conditional=carPalette()[3],
                      col.arrows="gray",
                      pch = c(16,1), lwd = 2, grid=TRUE,   
                      ellipse=FALSE,
                      overlaid=TRUE, new=TRUE, title=TRUE, ...){
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
    ellipse.args <- applyDefaults(ellipse, defaults=list(levels=0.5))
    if (!isFALSE(ellipse)) ellipse <- TRUE
    variable <- if (is.character(variable) & 1 == length(variable))
        variable  else deparse(substitute(variable))
    if(new && !overlaid) {
        op <- par(mfrow=c(1,2))
        on.exit(par(op))
    }
    # if(missing(labels))
    #     labels <- names(residuals(model)[!is.na(residuals(model))])
    else deparse(substitute(variable))
    if(attr(attr(model.frame(model), "terms"), "intercept") == 0)
        stop("Error---the 'lm' object must have an intercept")
    mod.mat <- model.matrix(model)
    var.names <- colnames(mod.mat)
    var <- which(variable == var.names)
    if (0 == length(var))
        stop(paste(variable, "is not a column of the model matrix."))
    response <- response(model)
    responseName <- responseName(model)
    if (is.null(weights(model)))
        wt <- rep(1, length(response))
    else wt <- weights(model)
    res0 <- lm(cbind(mod.mat[, var], response) ~ 1, weights=wt)$residual
    res  <- lsfit(mod.mat[, -var], cbind(mod.mat[, var], response),
                  wt = wt, intercept = FALSE)$residuals
    xlab <- paste(var.names[var], "| others")
    ylab <- paste(responseName, " | others")
    xlm <- c( min(res0[, 1], res[, 1]), max(res0[, 1], res[, 1]))
    ylm <- c( min(res0[, 2], res[, 2]), max(res0[, 2], res[, 2]))
    if(overlaid){
        mn <- if(title) paste("Marginal/Conditional plot of", var.names[var]) else NULL
        plot(res[, 1], res[, 2], xlab = xlab, ylab = ylab, type="n",
             main=mn,
             xlim=xlm, ylim=ylm,  ...)
        if(grid){
            grid(lty=1, equilogs=FALSE)
            box()}
        points(res0[, 1], res0[, 2], pch=pch[1], col=col.marginal)
        points(res[, 1], res[, 2], col=col.conditional, pch=pch[2], ...)
        arrows(res0[, 1], res0[, 2], res[, 1], res[, 2], length=0.125, col=col.arrows)
        abline(lsfit(res0[, 1], res0[, 2], wt = wt), col = col.marginal, lwd = lwd)
        abline(lsfit(res[, 1], res[, 2], wt = wt), col = col.conditional, lwd = lwd)
        if (ellipse) {
            ellipse.args1 <- c(list(res0, add=TRUE, plot.points=FALSE, col=col.marginal), ellipse.args)
            do.call(dataEllipse, ellipse.args1)
            ellipse.args1 <- c(list(res, add=TRUE, plot.points=FALSE, col=col.conditional), ellipse.args)
            do.call(dataEllipse, ellipse.args1)
        }
        showLabels(res0[, 1], res0[, 2], labels=labels,
                   method=id.method, n=id.n, cex=id.cex,
                   col=id.col, location=id.location)
        colnames(res) <- c(var.names[var], responseName)
        rownames(res) <- rownames(mod.mat)
        invisible(res)}
    else { # side.by.side plots
        mn <- if(title) paste("Marginal plot of", var.names[var]) else NULL
        plot(res0[, 1], res0[, 2], type="n",
             xlab = paste("Centered", var.names[var], sep=" "),
             ylab = paste("Centered", responseName, sep=" "),
             main=mn,
             xlim=xlm, ylim=ylm,  ...)
        if(grid){
            grid(lty=1, equilogs=FALSE)
            box()}
        points(res0[, 1], res0[, 2], pch=pch[1], col=col.marginal)
        abline(lsfit(res0[, 1], res0[, 2], wt = wt), col = col.marginal, lwd = lwd)
        if (ellipse) {
            ellipse.args1 <- c(list(res0, add=TRUE, plot.points=FALSE, col=col.marginal), ellipse.args)
            do.call(dataEllipse, ellipse.args1)
        }
        showLabels(res0[, 1], res0[, 2], labels=labels,
                   method=id.method, n=id.n, cex=id.cex,
                   col=id.col, location=id.location)
        colnames(res) <- c(var.names[var], responseName)
        rownames(res) <- rownames(mod.mat)
        mn <- if(title) paste("Added-Variable plot of", var.names[var]) else NULL
        plot(res[, 1], res[, 2], xlab = xlab, ylab = ylab, type="n",
             main=mn,
             xlim=xlm, ylim=ylm,  ...)
        if(grid){
            grid(lty=1, equilogs=FALSE)
            box()}
        points(res[, 1], res[, 2], col=col.conditional, pch=pch[2], ...)
        abline(lsfit(res[, 1], res[, 2], wt = wt), col = col.conditional, lwd = lwd)
        if (ellipse) {
            ellipse.args1 <- c(list(res, add=TRUE, plot.points=FALSE, col=col.conditional), ellipse.args)
            do.call(dataEllipse, ellipse.args1)
        }
        showLabels(res[, 1], res[, 2], labels=labels,
                   method=id.method, n=id.n, cex=id.cex,
                   col=id.col, location=id.location)
        invisible(res)}
}

