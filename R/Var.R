# Variance-Covariance matrices (J. Fox)

Var<-function(object, ...){
    UseMethod("Var")
    }
    
Var.lm<-function(object, diagonal=F){
    summary<-summary(object, corr=F)
    V<-(summary$sigma^2)*summary$cov.unscaled
    if (diagonal) diag(V) else V
    }
    
Var.glm<-function(object, diagonal=F){
    summary<-summary(object, corr = F)
    V<-summary$dispersion*summary$cov.unscaled
    if (diagonal) diag(V) else V
    }

Var.default<-function(object, diagonal=F, ...){
    # last modified 12 Dec 2000 by J. Fox
    V<-var(object, ...)
    if (diagonal) diag(V) else V
    }

    
