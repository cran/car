# CERES plots (J. Fox)

ceres.plots<-function(model, variable, ask=missing(variable), one.page=!ask, span=.5, ...){
    # last modified 23 Apr 2001 by J. Fox
    if (!missing(variable)){
        var<-if (is.character(variable) & 1==length(variable)) variable
            else deparse(substitute(variable))
        ceres.plot(model, var, line=line, main=main, ...)
        }
    else {
        vars<-predictor.names(model)
        vars<-if (is.null(model$contrasts)) vars
            else vars[sapply(model$contrasts[vars], is.null)]
        if (0==length(vars)) stop("No covariates to plot.")
        if (any(attr(terms(model),"order")>1)) {
            stop("ceres plots not available for models with interactions.")
            }
        if (ask) {
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
                    ceres.plot(model, var, span=span)
                    }
                }
            }
        else {
            if (one.page){
                save.mfrow <- par(mfrow=mfrow(length(vars)))
                on.exit(par(mfrow=save.mfrow))
                }
            for (var in vars){ 
                 ceres.plot(model, var, span=span, ...)
                 }
            }
        }
    }


ceres.plot<-function (model, ...) {
    UseMethod("ceres.plot")
    }

ceres.plot.lm<-function(model, variable, line=T, smooth=T, span=.5, iter, 
    las=1, col=palette()[2], pch=1, lwd=2, main="Ceres Plot", ...){
    # the lm method works with glm's too
    # last modified 1 Feb 2001 by J. Fox
    if (missing(iter)){
        iter<-if(("glm"==class(model)[1]) &&
                 ("gaussian"!=as.character(family(model))[1]))
                0
                else 3
            }    # use nonrobust smooth for non-gaussian glm
    require(modreg)
    var<-if (is.character(variable) & 1==length(variable)) variable
        else deparse(substitute(variable))
    mod.mat<-model.matrix(model)
    obs<-names(residuals(model))
    all.obs<-if (is.null(model$call$data)) obs else row.names(eval(model$call$data))
    xx<-rep(NA, length(all.obs))
    names(xx)<-all.obs
    vars<-predictor.names(model)
    if (is.na(match(var, vars))) stop(paste(var,"is not in the model."))
    if (!is.null(model$contrasts[[var]])) stop(paste(var,"is a factor."))
    vars<-vars[-match(var,vars)]
    if (any(attr(terms(model),"order")>1)) {
        stop("ceres plot not available for models with interactions.")
        }
    .x<-xvars<-NULL
    for (xvar in vars){
        if (is.null(model$contrasts[[xvar]])){
            xvars<-c(xvars,xvar)
            xx[obs]<-fitted.values(loess(as.formula(paste("mod.mat[,'",xvar,"']~mod.mat[,'",var,"']",sep=""))))
            .x<-cbind(.x, xx)
            }
        }
    if (is.null(xvars)) stop("There are no covariates.")
    n.x<-length(xvars)
    mf<-model.frame(model)
    rownames(.x)<-all.obs
    mf$.x<-.x[obs,]
    aug.model<-update(model, .~.+.x, data=mf)
    aug.mod.mat<-model.matrix(aug.model)
    coef<-coefficients(aug.model)
    k<-length(coef)
    posn<-k:(k-n.x+1)
    partial.res<-residuals.glm(aug.model, "partial")[,var] +
        aug.mod.mat[,posn] %*% as.matrix(coef[posn])
    plot(mod.mat[,var], partial.res, xlab=var, col=col, pch=pch,
        ylab=paste("CERES Residual(",response.name(model),")", sep=""),
        main=main, las=las)
    if (line) abline(lm(partial.res~mod.mat[,var]), lty=2, lwd=lwd, col=col)
    if (smooth) {
        lines(lowess(mod.mat[,var], partial.res, iter=iter, f=span), lwd=lwd, col=col)
        }
    }                    

ceres.plot.glm<-function(model, ...){
  # last modified 14 Dec 2000
  ceres.plot.lm(model, ...)
  }
