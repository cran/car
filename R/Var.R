# Variance-Covariance matrices (J. Fox)

# last modified 9 Nov 02 by J. Fox

Var<-function(object, ...){
    UseMethod("Var")
    }
    
Var.lm<-function(object, diagonal=FALSE, ...){
    summary<-summary(object, corr=FALSE)
    V<-(summary$sigma^2)*summary$cov.unscaled
    if (diagonal) diag(V) else V
    }
    
Var.glm<-function(object, diagonal=FALSE, ...){
    summary<-summary(object, corr = FALSE)
    V<-summary$dispersion*summary$cov.unscaled
    if (diagonal) diag(V) else V
    }

Var.default<-function(object, diagonal=FALSE, ...){
    # last modified 12 Dec 2000 by J. Fox
    V<-var(object, ...)
    if (diagonal) diag(V) else V
    }

    
