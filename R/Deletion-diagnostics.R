# Deletion diagonstics (J. Fox)

# last modified 2 Mar 03 by J. Fox

influence.glm<-function(model, do.coef = TRUE, ...) 
# slight changes to lm.influence to get right sigmas for glms,
# to return deviance and Pearson residuals,
# to handle na.exclude, and to work with versions >= 1.7.0
{
    version <- sapply(R.Version()[c("major","minor")], as.numeric) %*% c(10,1)
    if (is.empty.model(model$terms)) {
        warning("Can't compute influence on an empty model")
        return(NULL)
    }
    n <- as.integer(nrow(model$qr$qr))
    k <- as.integer(model$qr$rank)
    e <- na.omit(residuals(model, type="deviance")) # deviance res instead of weighted res
    select <- model$prior.weights != 0
    pear.res <- na.omit(residuals(model, type="pearson"))
    res <- if (version >= 17) { # version 1.7.0 or greater
                c(list(names=names(e)[select], dev.res=e[select], pear.res=pear.res[select]),
                    .Fortran("lminfl",
                    model$qr$qr,
                    n,
                    n,
                    k,
                    as.integer(do.coef),
                    model$qr$qraux,
                    e[select],
                    hat = double(n),
                    coefficients= if(do.coef) matrix(0, n, k) else double(1),
                    sigma = double(n),
                    DUP = FALSE, PACKAGE="base")[c("hat", "coefficients", "sigma")])
                }
            else {
                c(list(names=names(e)[select], dev.res=e[select], pear.res=pear.res[select]),
                    .Fortran("lminfl", model$qr$qr, n, n, k, model$qr$qraux, 
                        e[select], hat = double(n), coefficients = matrix(0, nr = n, nc = k), 
                        sigma = double(n), DUP = FALSE, PACKAGE = "base")[c("hat", "coefficients", "sigma")])
                }               
    if (is.null(model$na.action)) return(res)
    else res <- naresid(model$na.action, cbind(res$dev.res, res$pear.res, 
        res$hat, res$sigma, if (do.coef) res$coefficients))
    list(names=rownames(res), dev.res=res[,1], pear.res=res[,2], hat=res[,3], sigma=res[,4], 
        coefficients=if (do.coef) res[,-(1:4)])
}

influence.lm<-function(model, do.coef = TRUE, ...) 
# slight changes to lm.influence to return weighted residuals
# to handle na.exclude, and to work with versions >= 1.7.0
{
    version <- sapply(R.Version()[c("major","minor")], as.numeric) %*% c(10,1)
    if (is.empty.model(model$terms)) {
        warning("Can't compute influence on an empty model")
        return(NULL)
    }
    n <- as.integer(nrow(model$qr$qr))
    k <- as.integer(model$qr$rank)
    e <- na.omit(weighted.residuals(model))
    res <- if (version >= 17) { # version 1.7.0 or greater
                c(list(names=names(e), wt.res = e), 
                    .Fortran("lminfl",
                    model$qr$qr,
                    n,
                    n,
                    k,
                    as.integer(do.coef),
                    model$qr$qraux,
                    e,
                    hat = double(n),
                    coefficients= if(do.coef) matrix(0, n, k) else double(1),
                    sigma = double(n),
                    DUP = FALSE, PACKAGE="base")[c("hat", "coefficients", "sigma")])
                }
            else {
                c(list(names=names(e), wt.res = e), 
                    .Fortran("lminfl", model$qr$qr, n, n, k, model$qr$qraux, 
                    e, hat = double(n), coefficients = matrix(0, nr = n, nc = k), 
                    sigma = double(n), DUP = FALSE, PACKAGE = "base")[c("hat", "sigma", "coefficients")])
                }
    if (is.null(model$na.action)) return(res)
    else res <- naresid(model$na.action, 
        cbind(res$wt.res, res$hat, res$sigma, if (do.coef) res$coefficients))
    list(names=rownames(res), wt.res=res[,1], hat=res[,2], sigma=res[,3], 
        coefficients=if (do.coef) res[,-(1:3)])
}

hatvalues<-function(model, ...){
    UseMethod("hatvalues")
    }
    
hatvalues.lm<-function(model, infl=influence(model, do.coef=FALSE), names=infl$names, ...){
    hat<-infl$hat
    names(hat)<-names
    hat
    }

rstudent<-function(model, ...){
    UseMethod("rstudent")
    }

rstudent.lm<-function(model, infl=influence(model, do.coef=FALSE), names=infl$names, ...){
    rstud<-infl$wt.res/(infl$sigma*(1 - infl$hat)^.5)
    names(rstud)<-names
    rstud
    }
    
rstudent.glm<-function(model, infl=influence(model, do.coef=FALSE), names=infl$names, ...){
    rstud<-sign(infl$dev)*sqrt(infl$dev.res^2 + (infl$hat * infl$pear.res^2)/(1 - infl$hat))
    rstud<-if (any(family(model)$family == c("binomial", "poisson"))) rstud
            else rstud/infl$sigma
    names(rstud)<-names
    rstud
    }
    
influence<-function(model, ...){
    UseMethod("influence")
    }

cookd<-function(model, ...){
    UseMethod("cookd")
    }
    
cookd.lm<-function(model, infl=influence(model, do.coef=FALSE), sumry=summary(model), names=infl$names, ...){
    sigma<-sumry$sigma
    df<-model$rank
    res<-infl$wt.res/(sigma*(1 - infl$hat)^.5)
    cookd<-(res^2 * infl$hat)/(df*(1 - infl$hat))
    names(cookd)<-names
    cookd
    }

cookd.glm<-function(model, infl=influence(model, do.coef=FALSE), sumry=summary(model), names=infl$names, ...){   
    phi<-sumry$dispersion
    df<-model$rank
    cookd<-(infl$pear.res^2 * infl$hat)/(phi*df*(1 - infl$hat)^2)
    names(cookd)<-names
    cookd
    }
    
dfbeta<-function(model, ...){
    UseMethod("dfbeta")
    }
    
dfbeta.lm<-function(model, infl=influence(model), names=infl$names, ...){
    b<-infl$coefficients
    dimnames(b)<-list(names, variable.names(model))
    b
    }
    
dfbetas<-function(model, ...){
    UseMethod("dfbetas")
    }
        
dfbetas.lm<-function(model, infl=influence(model), sumry=summary(model), names=infl$names, ...){ 
    sb<-(diag(sumry$cov.unscaled))^.5
    b<-dfbeta(model)/(infl$sigma %o% sb)
    dimnames(b)<-list(names, variable.names(model))
    b    
    }
