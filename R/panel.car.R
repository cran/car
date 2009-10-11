# panel function for use with coplot (J. Fox)

# last modified 2 April 02

panel.car<-function(x, y, col, pch, cex=1, span=.5, lwd=2,
    regression.line=lm, lowess.line=TRUE,...){
    # last modified 10 Dec 2001 by J. Fox
    points(x, y, col=col, pch=pch, cex=cex)
    if (is.function(regression.line)) reg.line(regression.line(y~x), 
        lty=2, lwd=lwd, col=col, ...)
    if (lowess.line) lines(lowess(na.omit(as.data.frame(cbind(x,y))), f=span), 
        col=col, lwd=lwd, ...)
    }
