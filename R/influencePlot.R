# changed point marking, 25 November 2009 by S. Weisberg
#  deleted the cutoff for Cook's D, and the coloring of the circles
#  inserted default labeling of the id.n largest Cook D.
# 13 January 2009: changed to label points by all of hatvalues,
#  studentized residuals, and Cook's Ds. J. Fox
# 14 April 2010: set id.n = 0. J. Fox
# 23 April 2010: rewrote point marking, S. Weisberg
# 10 May 2010: fixed computation of n
# 2014-04-19: use labels for returned table rownames. J. Fox
# 2015-11-06: now returns Cook's distance, not its square root.  S. Weisberg
# 2017-02-12: consolidated id argument. J. Fox
# 2017-11-30: substitute carPalette() for palette(). J. Fox
# 2019-01-02: added lmerMod method. J. Fox

# moved from Rcmdr 5 December 2006

influencePlot <- function(model, ...){
    UseMethod("influencePlot")
}

influencePlot.lm <- function(model, scale=10,  
                             xlab="Hat-Values", ylab="Studentized Residuals",
                             id=TRUE, ...){
    id <- applyDefaults(id, defaults=list(method="noteworthy", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
    if (isFALSE(id)){
        id.n <- 0
        id.method <- "none"
        labels <- id.cex <- id.col <- id.location <- NULL
    }
    else{
        labels <- id$labels
        if (is.null(labels)) labels <- names(na.omit(residuals(model)))
        id.method <- id$method
        id.n <- if ("identify" %in% id.method) Inf else id$n
        id.cex <- id$cex
        id.col <- id$col
        id.location <- id$location
    }
    
    hatval <- hatvalues(model)
    rstud <- rstudent(model)
    if (missing(labels)) labels <- names(rstud)
    cook <- sqrt(cooks.distance(model))
    scale <- scale/max(cook, na.rm=TRUE)
    p <- length(coef(model))
    n <- sum(!is.na(rstud))
    plot(hatval, rstud, xlab=xlab, ylab=ylab, type='n', ...)
    abline(v=c(2, 3)*p/n, lty=2)
    abline(h=c(-2, 0, 2), lty=2)
    points(hatval, rstud, cex=scale*cook, ...)
    if(id.method == "noteworthy"){
        which.rstud <- order(abs(rstud), decreasing=TRUE)[1:id.n]
        which.cook <- order(cook, decreasing=TRUE)[1:id.n]
        which.hatval <- order(hatval, decreasing=TRUE)[1:id.n]
        which.all <- union(which.rstud, union(which.cook, which.hatval))
        id.method <- which.all
    }
    noteworthy <- if (!isFALSE(id)) showLabels(hatval, rstud, labels=labels, method=id.method, 
                             n=id.n, cex=id.cex, col=id.col, location = id.location)
    else NULL
    if (length(noteworthy > 0)){
        result <- data.frame(StudRes=rstud[noteworthy], Hat=hatval[noteworthy],
                             CookD=cook[noteworthy]^2)
        if (is.numeric(noteworthy)) rownames(result) <- labels[noteworthy]
        return(result)
    }
    else return(invisible(NULL))
}

influencePlot.lmerMod <- function(model, ...){
  influencePlot.lm(model, ...)
}
