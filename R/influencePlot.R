# moved from Rcmdr 5 December 2006

# the following function adapted from Fox, An R and S-PLUS Companion to Applied Regression
influencePlot <- function(model, ...){
    UseMethod("influencePlot")
    }

influencePlot.lm <- function(model, scale=10, col=c(1,2),
    labels=names(rstud), identify.cex=par("cex"), identify.col=par("col"), ...){
    hatval <- hatvalues(model)
    rstud <- rstudent(model)
    cook <- sqrt(cookd(model))
    scale <- scale/max(cook, na.rm=TRUE)
    p <- length(coef(model))
    n <- length(rstud)
    cutoff <- sqrt(4/(n - p))
    plot(hatval, rstud, xlab='Hat-Values',
        ylab='Studentized Residuals', type='n', ...)
    abline(v=c(2, 3)*p/n, lty=2)
    abline(h=c(-2, 0, 2), lty=2)
    points(hatval, rstud, cex=scale*cook, 
            col=ifelse(cook > cutoff, col[2], col[1]))
    if (labels[1] != FALSE) identify(hatval, rstud, labels, col=identify.col,
        cex=identify.cex)
    }
