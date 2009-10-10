# CERES plots (J. Fox)

# last modified 27 Apr 04 by J. Fox

ceres.plots<-function(model, variable, ask=missing(variable), one.page=!ask, span=.5, ...){
    # last modified 2 Aug 2001 by J. Fox
    if(!is.null(class(model$na.action)) && 
        class(model$na.action) == 'exclude') class(model$na.action) <- 'omit'
    if (!missing(variable)){
        var<-if (is.character(variable) & 1==length(variable)) variable
            else deparse(substitute(variable))
        ceres.plot(model, var, ...)
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

ceres.plot.lm<-function(model, variable, line=TRUE, smooth=TRUE, span=.5, iter, 
    las=par("las"), col=palette()[2], pch=1, lwd=2, main="Ceres Plot", ...){
    # the lm method works with glm's too
    # last modified 27 Apr 2004 by J. Fox
    expand.model.frame <- function (model, extras, envir = environment(formula(model)),
        na.expand = FALSE){  # modified version of R base function
        f <- formula(model)
        data <- eval(model$call$data, envir)
        ff <- foo ~ bar + baz
        if (is.call(extras)) 
            gg <- extras
        else gg <- parse(text = paste("~", paste(extras, collapse = "+")))[[1]]
        ff[[2]] <- f[[2]]
        ff[[3]][[2]] <- f[[3]]
        ff[[3]][[3]] <- gg[[2]]
        if (!na.expand) {
            naa <- model$call$na.action
            subset <- model$call$subset
            rval <- if (is.null(data)) eval(call("model.frame", ff, # modified
                subset = subset, na.action = naa), envir)           #  lines
            else eval(call("model.frame", ff, data = data,          #
                subset = subset, na.action = naa), envir)           #
            }
        else {
            subset <- model$call$subset
            rval <- eval(call("model.frame", ff, data = data, subset = subset, 
                na.action = I), envir)
            oldmf <- model.frame(model)
            keep <- match(rownames(oldmf), rownames(rval))
            rval <- rval[keep, ]
            class(rval) <- "data.frame"
            }
        return(rval)
        }
    if(!is.null(class(model$na.action)) && 
        class(model$na.action) == 'exclude') class(model$na.action) <- 'omit'
    if (missing(iter)){
        iter<-if(("glm"==class(model)[1]) &&
                 ("gaussian"!=as.character(family(model))[1]))
                0
                else 3
            }    # use nonrobust smooth for non-gaussian glm
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
    mf<-na.omit(expand.model.frame(model, all.vars(formula(model))))
    rownames(.x)<-all.obs
    mf$.x<-.x[obs,]
    aug.model <- update(model, . ~ . + .x, data=mf, subset=NULL)
    aug.mod.mat<-model.matrix(aug.model)
    coef<-coefficients(aug.model)
    k<-length(coef)
    posn<-k:(k-n.x+1)
    partial.res<-residuals.glm(aug.model, "partial")[,var] +
        aug.mod.mat[,posn] %*% as.matrix(coef[posn])
    plot(mod.mat[,var], partial.res, xlab=var, col=col, pch=pch,
        ylab=paste("CERES Residual(",responseName(model),")", sep=""),
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
