# Added-Variable plots (J. Fox)

avp<-function(...) av.plots(...)

av.plots<-function(model, variable, ask=missing(variable), one.page=!ask, ...){
    # last modified 23 Apr 2001 by J. Fox
    if (!missing(variable)){
        var<-if (is.character(variable) & 1==length(variable)) variable
            else deparse(substitute(variable))
        av.plot(model, var, ...)
        }
    else {
        vars<-colnames(model.matrix(model))
        if (ask) {
            repeat{
                selection<-menu(vars)
                if (selection==0) break
                else var<-vars[selection]
                av.plot(model, var, ...)
                }
            }
        else {
            if (one.page){
                save.mfrow <- par(mfrow=mfrow(length(vars)))
                on.exit(par(mfrow=save.mfrow))
                }
            for (var in vars) av.plot(model, var, ...)
            }
        }
    }


av.plot<-function (model, ...) {
    UseMethod("av.plot")
    }

av.plot.lm<-function(model, variable, labels=names(residuals(model)[!is.na(residuals(model))]), 
    identify.points=T, las=1, col=palette()[2], pch=1, lwd=2, main="Added-Variable Plot"){
    #last modified 29 July 2001 by J. Fox
    variable<-if (is.character(variable) & 1==length(variable)) variable
        else deparse(substitute(variable))
    mod.mat<-model.matrix(model)
    var.names<-colnames(mod.mat)
    var<-which(variable==var.names)
    if (0==length(var)) stop(paste(variable,"is not a column of the model matrix."))
    response<-response(model)
    response.name<-response.name(model)
    if (is.null(weights(model))) wt<-rep(1, length(response))
        else wt<-weights(model)
    res<-lsfit(mod.mat[,-var], cbind(mod.mat[,var], response), wt=wt,    
        intercept=F)$residuals
    plot(res[,1], res[,2], xlab=paste(var.names[var],"| others"), 
        ylab=paste(response.name," | others"), main=main, las=las, col=col, pch=pch)
    abline(lsfit(res[,1], res[,2], wt=wt), col=col, lwd=lwd)
    if (identify.points) identify(res[,1], res[,2], labels)
    }


av.plot.glm<-function(model, variable, labels=names(residuals(model)[!is.na(residuals(model))]), 
    identify.points=T, las=1, col=palette()[2], pch=1, lwd=2, main="Added-Variable Plot",
    type=c("Wang", "Weisberg")){
    #last modified 14 Dec 2001 by J. Fox
    type<-match.arg(type)
    variable<-if (is.character(variable) & 1==length(variable)) variable
        else deparse(substitute(variable))
    mod.mat<-model.matrix(model)
    var.names<-colnames(mod.mat)
    var<-which(variable==var.names)
    if (0==length(var)) stop(paste(variable,"is not a column of the model matrix."))
    response<-response(model)
    response.name<-response.name(model)
    wt<-model$prior.weights
    mod<-glm(response~mod.mat[,-var]-1, weights=wt, family=family(model))
    res.y<-residuals(mod, type="pearson")
    wt<-if (type=="Wang") wt*model$weights else wt
    res.x<-lsfit(mod.mat[,-var], mod.mat[,var], wt=wt,    
        intercept=F)$residuals
    plot(res.x, res.y, xlab=paste(var.names[var],"| others"), 
        ylab=paste(response.name," | others"), main=main, las=las, col=col, pch=pch)
    abline(lsfit(res.x, res.y, wt=wt), col=col, lwd=lwd)
    if (identify.points) identify(res.x, res.y, labels)
    }
