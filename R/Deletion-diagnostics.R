# Deletion diagonstics (J. Fox)

influence.glm<-function(model) 
# last modified 27/1/2001 by J. Fox
# slight change to lm.influence to get right sigmas for glms
# and to return deviance and Pearson residuals
{
    if (is.empty.model(model$terms)) {
        warning("Can't compute influence on an empty model")
        return(NULL)
    }
    n <- as.integer(nrow(model$qr$qr))
    k <- as.integer(model$qr$rank)
    e <- residuals(model, type="deviance") # deviance res instead of weighted res
    names(e)<-NULL
    c(.Fortran("lminfl", model$qr$qr, n, n, k, model$qr$qraux, 
        e, hat = double(n), coefficients = matrix(0, nr = n, 
            nc = k), sigma = double(n), DUP = FALSE, PACKAGE = "base")[c("hat", 
        "coefficients", "sigma")], 
        list(dev.res=e, pear.res=residuals(model, type="pearson")))
}

influence.lm<-function(model) 
# slight change to lm.influence to return weighted residuals
{
    if (is.empty.model(model$terms)) {
        warning("Can't compute influence on an empty model")
        return(NULL)
    }
    n <- as.integer(nrow(model$qr$qr))
    k <- as.integer(model$qr$rank)
    e <- weighted.residuals(model)
    names(e)<-NULL
    c(.Fortran("lminfl", model$qr$qr, n, n, k, model$qr$qraux, 
        e, hat = double(n), coefficients = matrix(0, nr = n, 
            nc = k), sigma = double(n), DUP = FALSE, PACKAGE = "base")[c("hat", 
        "coefficients", "sigma")], list(wt.res=e))
}

hatvalues<-function(model, ...){
    UseMethod("hatvalues")
    }
    
hatvalues.lm<-function(model, infl=influence(model)){
    hat<-infl$hat
    names(hat)<-case.names(model)
    hat
    }

rstudent<-function(model, ...){
    UseMethod("rstudent")
    }

rstudent.lm<-function(model, infl=influence(model)){
    rstud<-infl$wt.res/(infl$sigma*(1 - infl$hat)^.5)
    names(rstud)<-case.names(model)
    rstud
    }
    
rstudent.glm<-function(model, infl=influence(model)){
    rstud<-sign(infl$dev)*sqrt(infl$dev.res^2 + (infl$hat * infl$pear.res^2)/(1 - infl$hat))
    rstud<-if (any(family(model)$family == c("binomial", "poisson"))) rstud
            else rstud/infl$sigma
    names(rstud)<-case.names(model)
    rstud
    }
    
influence<-function(model, ...){
    UseMethod("influence")
    }

cookd<-function(model, ...){
    UseMethod("cookd")
    }
    
cookd.lm<-function(model, infl=influence(model), sumry=summary(model)){
    sigma<-sumry$sigma
    df<-model$rank
    res<-infl$wt.res/(sigma*(1 - infl$hat)^.5)
    cookd<-(res^2 * infl$hat)/(df*(1 - infl$hat))
    names(cookd)<-case.names(model)
    cookd
    }

cookd.glm<-function(model, infl=influence(model), sumry=summary(model)){   
    phi<-sumry$dispersion
    df<-model$rank
    cookd<-(infl$pear.res^2 * infl$hat)/(phi*df*(1 - infl$hat)^2)
    names(cookd)<-case.names(model)
    cookd
    }
    
dfbeta<-function(model, ...){
    UseMethod("dfbeta")
    }
    
dfbeta.lm<-function(model, infl=influence(model)){
    b<-infl$coefficients
    dimnames(b)<-list(case.names(model), variable.names(model))
    b
    }
    
dfbetas<-function(model, ...){
    UseMethod("dfbetas")
    }
        
dfbetas.lm<-function(model, infl=influence(model), sumry=summary(model)){ 
    sb<-(diag(sumry$cov.unscaled))^.5
    b<-dfbeta(model)/(infl$sigma %o% sb)
    dimnames(b)<-list(case.names(model), variable.names(model))
    b    
    }
