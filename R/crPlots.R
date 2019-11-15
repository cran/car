# Component + Residual Plots (J. Fox)

# modified 9 October 2009 by J. Fox
# modified 25 November 2009 by S. Weisberg to change 
#   variable specification, layout and point marking
# modified 1 January 2009 by J. Fox
#   to set default id.n=0
# changed showLabels args 15 April 2010 S. Weisberg
# added grid, 10 May 2010
# modified 2 Sept 2010 by S. Weisberg, made colors, axes lables, and
# arguments more consistent with other functions; ... passes args to plot
# and boxplot.
# 16 June 2011 allow layout=NA, in which case the layout is not set in this
#  function, so it is the responsibility of the user
# 14 Sept 2012 use the ScatterplotSmoothers in car
# 19 Sept 2012 restore smooth and span args
# 20 Aug 2013 replace residuals.glm() with residuals(). John
# 2017-02-11: consolidated id and smooth arguments. John
# 2017-11-30: substitute carPalette() for palette(). J. Fox

# these functions to be rewritten; simply renamed for now

crp<-function(...) crPlots(...)

crPlots<-function(model, terms= ~ ., layout=NULL, ask, main, ...){
    terms <- if(is.character(terms)) paste("~", terms) else terms
    vform <- update(formula(model), terms)
    if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
        stop("Only predictors in the formula can be plotted.")
    mf <- attr(model.frame(model), "terms")
    terms <- attr(mf, "term.labels") # this is a list
    vterms <- attr(terms(vform), "term.labels")
    if (any(attr(terms(model),"order")>1)) {
        stop("C+R plots not available for models with interactions.")}
    nt <- length(vterms)
    if (nt == 0) stop("No plots specified")
    if (missing(main)) main <- if (nt == 1) "Component + Residual Plot" else "Component + Residual Plots" 
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
    if(!is.null(class(model$na.action)) && 
       inherits(model$na.action, 'exclude')) class(model$na.action) <- 'omit'
    for(term in vterms) 
        crPlot(model, term, ...)
    mtext(side=3, outer=TRUE, main, cex=1.2)
    invisible(0)
}


crPlot<-function (model, ...) {
    UseMethod("crPlot")
}

crPlot.lm<-function(model, variable, id=FALSE,
                    order=1, line=TRUE, smooth=TRUE, col=carPalette()[1], col.lines=carPalette()[-1],
                    xlab, ylab, pch=1, lwd=2, grid=TRUE, ...) { 
    # method also works for glm objects
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
    smoother.args <- applyDefaults(smooth, defaults=list(smoother=loessLine), type="smooth")
    if (!isFALSE(smoother.args)) {
        smoother <- smoother.args$smoother 
        smoother.args$smoother <- NULL
    }
    else smoother <- "none"
    if(!is.null(class(model$na.action)) && 
       inherits(model$na.action, 'exclude')) class(model$na.action) <- 'omit'
    var<-if (is.character(variable) & 1==length(variable)) variable
    else deparse(substitute(variable))
    xlab <- if(!missing(xlab)) xlab else var
    ylab <- if(!missing(ylab)) ylab else
        paste("Component+Residual(", responseName(model),")", sep="")
    terms<-predictor.names(model)
    if (is.na(match(var, terms))) stop(paste(var,"is not in the model."))
    if (any(attr(terms(model),"order")>1)) {
        stop("C+R plots not available for models with interactions.")
    }
    if (!is.null(model$contrasts[[var]])){
        partial.res<-residuals(model,"partial")
        .x<-model.frame(model)[,var]
        boxplot(partial.res[,var]~.x, xlab=xlab,
                ylab=ylab, ...)
        return(invisible())
    }
    .x<-if (df.terms(model, var)>1) predict(model, type="terms", term=var)
    else model.matrix(model)[,var]
    if (order==1){          # handle first-order separately for efficiency
        partial.res<-residuals(model,"partial")
        plot(.x, partial.res[,var], type="n", xlab=xlab,
             ylab=ylab, ...)
        if(grid){
            grid(lty=1, equilogs=FALSE)
            box()}
        points(.x, partial.res[,var], col=col, pch=pch)
        if (line) abline(lm(partial.res[,var]~.x), lty=2, lwd=lwd, col=col.lines[1])
        if (is.function(smoother)) {
            smoother(.x, partial.res[,var], col=col.lines[2], log.x=FALSE,
                     log.y=FALSE, spread=FALSE, smoother.args=smoother.args)
        }
        showLabels(.x, partial.res[,var], labels=labels, 
                   method=id.method, n=id.n, cex=id.cex,
                   col=id.col, location=id.location)
    }
    else {
        if (df.terms(model, var) > 1) 
            stop(paste("Order", order, "C+R plot not available for a term with > 1 df:", var))
        aug.model<-update(model, 
                          as.formula(paste(".~.-",var,"+poly(",var,",",order,")")))
        partial.res<-residuals(aug.model, "partial")
        last<-ncol(partial.res)
        plot(.x, partial.res[,last], xlab=xlab, 
             ylab=ylab, type="n", ...)
        if(grid){
            grid(lty=1, equilogs=FALSE)
            box()}
        points(.x, partial.res[,last], col=col, pch=pch)
        if (line) abline(lm(partial.res[,last]~.x), lty=2, lwd=lwd, col=col.lines[1])
        if (is.function(smoother)) {
            smoother(.x, partial.res[, last], col=col.lines[2], log.x=FALSE,
                     log.y=FALSE, spread=FALSE, smoother.args=smoother.args)
        }
        showLabels(.x, partial.res[,last], labels=labels, 
                   method=id.method, n=id.n, cex=id.cex,
                   col=id.col, location=id.location)
    }          
}
