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
  if(is.null(layout)){
   layout <- switch(min(nt,9), c(1,1), c(1,2), c(2,2), c(2,2),
                               c(3,2), c(3,2), c(3,3), c(3,3), c(3,3))}
  ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
  if(nt > 1){
     op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
            oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
     on.exit(par(op))}
	if(!is.null(class(model$na.action)) && 
		class(model$na.action) == 'exclude') class(model$na.action) <- 'omit'
  for(term in vterms) 
		crPlot(model, term, ...)
	mtext(side=3, outer=TRUE, main, cex=1.2)
	invisible(0)
	}

	
crPlot<-function (model, ...) {
	UseMethod("crPlot")
}

crPlot.lm<-function(model, variable, 
  id.method = list(abs(residuals(model, type="pearson")), "x"),
  labels, 
  id.n = if(id.method[1]=="identify") Inf else 0,
  id.cex=1, id.col=palette()[1],
  order=1, line=TRUE, smooth=TRUE,
	iter, span=.5, 
  col=palette()[1], col.lines=palette()[-1],
  xlab, ylab, pch=1, lwd=2, grid=TRUE, ...) { 
	# method also works for glm objects
	if(missing(labels)) labels <- names(residuals(model))
	if(!is.null(class(model$na.action)) && 
		class(model$na.action) == 'exclude') class(model$na.action) <- 'omit'
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
		partial.res<-residuals.glm(model,"partial")
		.x<-model.frame(model)[,var]
		boxplot(partial.res[,var]~.x, xlab=xlab,
			ylab=ylab, ...)
		return(invisible())
	}
	if (missing(iter)){
		iter<-if(("glm"==class(model)[1]) &&
				("gaussian"!=as.character(family(model))[1]))
				0
			else 3
	}    # use nonrobust smooth for non-gaussian glm
	.x<-if (df.terms(model, var)>1) predict(model, type="terms", term=var)
		else model.matrix(model)[,var]
	if (order==1){          # handle first-order separately for efficiency
		partial.res<-residuals.glm(model,"partial")
		plot(.x, partial.res[,var], type="n", xlab=xlab,
      ylab=ylab, ...)
	  if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
		points(.x, partial.res[,var], col=col, pch=pch)
		if (line) abline(lm(partial.res[,var]~.x), lty=2, lwd=lwd, col=col.lines[1])
		if (smooth) {
			lines(lowess(.x, partial.res[,var], iter=iter, f=span), lwd=lwd, 
            col=col.lines[2])
		}
		showLabels(.x, partial.res[,var], labels=labels, 
            id.method=id.method, id.n=id.n, id.cex=id.cex,
            id.col=id.col)
	}
	else {
		if (df.terms(model, var)>1) 
			stop(paste("Order", order, "C+R plot not available for a term with > 1 df:", var))
		aug.model<-update(model, 
			as.formula(paste(".~.-",var,"+poly(",var,",",order,")")))
		partial.res<-residuals.glm(aug.model, "partial")
		last<-ncol(partial.res)
		plot(.x, partial.res[,last], xlab=xlab, 
			ylab=ylab, type="n", ...)
	  if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
		points(.x, partial.res[,last], col=col, pch=pch)
		if (line) abline(lm(partial.res[,last]~.x), lty=2, lwd=lwd, col=col.lines[1])
		if (smooth) {
			lines(lowess(.x, partial.res[,last], iter=iter, f=span), lwd=lwd, 
         col=col.lines[2])
		}
		showLabels(.x, partial.res[,last], labels=labels, 
            id.method=id.method, id.n=id.n, id.cex=id.cex,
            id.col=id.col)
	}          
}

crPlot.glm<-function(model, ...){
	crPlot.lm(model, ...)
	}
