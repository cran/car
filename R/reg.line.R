# draw regression line from model to extremes of fit (J. Fox)
 
reg.line<-function(mod, col=palette()[2], lwd=2, lty=1, ...){
    # last modified 1 Feb 2001 by J. Fox
    coef<-coefficients(mod)
    if (length(coef) != 2) error(" Requires simple linear regression.")
    x<-model.matrix(mod)[,2]
    y<-fitted.values(mod)
    min<-which.min(x)
    max<-which.max(x)
    lines(c(x[min],x[max]),c(y[min],y[max]), col=col, lty=lty, lwd=lwd, ...)
    }
