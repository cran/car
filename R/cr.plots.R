# Component + Residual Plots (J. Fox)

# last modified 2 Nov 04 by J. Fox

crp<-function(...) cr.plots(...)

cr.plots<-function(model, variable, ask=missing(variable), one.page=!ask, span=.5, ...){
    # last modified 2 Aug 2001 by J. Fox
    if(!is.null(class(model$na.action)) && 
        class(model$na.action) == 'exclude') class(model$na.action) <- 'omit'
    if (!missing(variable)){
        var<-if (is.character(variable) & 1==length(variable)) variable
            else deparse(substitute(variable))
        cr.plot(model, var, span=span, ...)
        }
    else {
        vars<-predictor.names(model)
        if (0==length(vars)) stop("No covariates to plot.")
        else if (ask) {
            repeat{
                selection<-menu(c(paste("Change span = ",span),vars))
                if (selection==0) break
                if (selection==1) {
                    span<-eval(parse(text=readline(prompt="span: ")))
                    if ((!is.numeric(span)) || length(span)>1 || span<0
                        || span>1) stop("Span must be between 0 and 1")
                    }
                else {
                    var<-vars[selection-1]
                    cr.plot(model, var, span=span, ...)
                    }
                }
            }
        else {
            if (one.page){
                save.mfrow <- par(mfrow=mfrow(length(vars)))
                on.exit(par(mfrow=save.mfrow))
                    }
            for (var in vars){ 
                 cr.plot(model, var, span=span, ...)
                }
            }
        }
    }
    
    
cr.plot<-function (model, ...) {
    UseMethod("cr.plot")
    }


cr.plot.lm<-function(model, variable, order=1, line=TRUE, smooth=TRUE,
    iter, span=.5, las=par("las"), col=palette()[2], pch=1, lwd=2,
    main="Component+Residual Plot", ...) {
    # last modified 20 Feb 2002 by J. Fox
    # method also works for glm objects
    if(!is.null(class(model$na.action)) && 
        class(model$na.action) == 'exclude') class(model$na.action) <- 'omit'
    var<-if (is.character(variable) & 1==length(variable)) variable
        else deparse(substitute(variable))
    vars<-predictor.names(model)
    if (is.na(match(var, vars))) stop(paste(var,"is not in the model."))
    if (any(attr(terms(model),"order")>1)) {
        stop("C+R plots not available for models with interactions.")
        }
    if (!is.null(model$contrasts[[var]])){
        partial.res<-residuals.glm(model,"partial")
        .x<-model.frame(model)[,var]
        boxplot(partial.res[,var]~.x, xlab=var,
            ylab=paste("Component+Residual(", responseName(model),")", sep=""),
            main=main)
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
        plot(.x, partial.res[,var], xlab=var, 
            ylab=paste("Component+Residual(", responseName(model),")", sep=""),
            las=las, col=col, pch=pch, main=main)
        if (line) abline(lm(partial.res[,var]~.x), lty=2, lwd=lwd, col=col)
        if (smooth) {
            lines(lowess(.x, partial.res[,var], iter=iter, f=span), lwd=lwd, col=col)
            }
        }
        else {
            if (df.terms(model, var)>1) 
               stop(paste("Order", order, "C+R plot not available for a term with > 1 df:", var))
            aug.model<-update(model, 
                as.formula(paste(".~.-",var,"+poly(",var,",",order,")")))
            partial.res<-residuals.glm(aug.model, "partial")
            last<-ncol(partial.res)
            plot(.x, partial.res[,last], xlab=var, 
                ylab=paste("Component+Residual(", responseName(model),")", sep=""),
                las=las, col=col, pch=pch, main=main)
            if (line) abline(lm(partial.res[,last]~.x), lty=2, lwd=lwd, col=col)
            if (smooth) {
                lines(lowess(.x, partial.res[,last], iter=iter, f=span), lwd=lwd, col=col)
                }
            }          
    }

cr.plot.glm<-function(model, ...){
  # last modified 14 Dec 2000 by J. Fox
  cr.plot.lm(model, ...)
 }
