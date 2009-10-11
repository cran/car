# draw regression line from model to extremes of fit (J. Fox)
 
reg.line<-function(mod, col=palette()[2], lwd=2, lty=1, ...){
    # last modified 11 Oct 2007 by J. Fox
    if(!is.null(class(mod$na.action)) && 
        class(mod$na.action) == 'exclude') class(mod$na.action) <- 'omit'
    coef<-coefficients(mod)
    if (length(coef) != 2) stop(" Requires simple linear regression.")
    x<-model.matrix(mod)[,2]
    y<-fitted.values(mod)
    min<-which.min(x)
    max<-which.max(x)
    lines(c(x[min],x[max]),c(y[min],y[max]), col=col, lty=lty, lwd=lwd, ...)
    }
