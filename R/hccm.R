# Heteroscedasticity-corrected standard errors (White adjustment) (J. Fox)

hccm<-function(model, ...){
    #last modified 12 Dec 2000 by J. Fox
    UseMethod("hccm")
    }
 
hccm.lm<-function(model, type=c("hc3", "hc0", "hc1", "hc2")) {
    #last modified 12 Dec 2000 by J. Fox
    if (!is.null(weights(model))) stop("requires unweighted lm")
    type<-match.arg(type)
    sumry<-summary(model, corr = FALSE)
    s2<-sumry$sigma^2
    V<-sumry$cov.unscaled
    if (type == FALSE) return(s2*V)
    e<-residuals(model)
    X<-model.matrix(model)
    df.res<-df.residual(model)
    factor<-switch(type,
        hc0=1,
        hc1=length(e)/df.res,
        hc2=1/(1-hat(X)),
        hc3=1/(1-hat(X))^2)
    V %*% t(X) %*% apply(X, 2, "*", (e^2)/factor) %*% V
    }
    
hccm.default<-function(model, ...){
    #last modified 12 Dec 2000 by J. Fox
    stop("requires an lm object")
    }
