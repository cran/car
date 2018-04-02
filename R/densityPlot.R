# checked in 2013-06-05 by J. Fox
# 2014-09-04: J. Fox: empty groups produce warning rather than error
# 2016-10-16: J. Fox: add option for adaptive kernel.
# 2016-11-26: J. Fox: rejig for pure-R adaptive kernel
# 2017-02-12: J. Fox: make adaptive kernel the default, consolidate legend args.
# 2017-11-30: substitute carPalette() for palette(). J. Fox

densityPlot <- function(x, ...){
    UseMethod("densityPlot")
}

densityPlot.default <- function (x, g, method=c("adaptive", "kernel"), 
    bw=if (method == "adaptive") bw.nrd0 else "SJ", adjust=1,
    kernel,
    xlim,
    ylim,
    normalize=FALSE,
    xlab=deparse(substitute(x)), ylab="Density", main="",
    col=carPalette(), lty=seq_along(col), lwd=2, grid=TRUE,
    legend=TRUE, show.bw=FALSE,
    rug=TRUE, ...) {
    norm <- function(x, y){
        n <- length(x)
        x0 <- diff(range(x))/(n - 1)
        y0 <- (y[1:(n-1)] + y[2:n])/2
        exp(log(y) - log(sum(y0*x0)))
    }
    legend <- applyDefaults(legend, defaults=list(location="topright", title=deparse(substitute(g))), type="legend")
    if (isFALSE(legend)){
        legend.title <- ""
        legend.location <- "topleft"
    }
    else{
        legend.title <- legend$title
        legend.location <- legend$location
    }
    method <- match.arg(method)
    if (method == "kernel"){
        kernel <- if (missing(kernel)) "gaussian"
        else match.arg(kernel, c("gaussian", "epanechnikov", "rectangular", 
                                      "triangular", "biweight", "cosine", "optcosine"))
    }
    else{
        if(missing(kernel)) kernel <- dnorm
        if (!is.function(kernel)) stop("for the adaptive kernel estimator, the kernel argument must be a function")
    }
    force(ylab)
    force(xlab)
    if (!is.numeric(x)) stop("argument x must be numeric")
    if (missing(g)) {
        density <- if (method == "adaptive") adaptiveKernel(x, bw=bw, adjust=adjust, ...)
            else density(x, bw=bw, adjust=adjust, kernel=kernel, ...)
        if (normalize) density$y <- norm(density$x, density$y)
        if (missing(xlim)) xlim <- range(density$x)
        if (missing(ylim)) ylim <- c(0, max(density$y))
        if (show.bw) xlab <- paste(xlab, " (bandwidth = ", format(density$bw), ")", sep="")
        plot(xlim, ylim, xlab=xlab, ylab=ylab, main=main, type="n")
        if (rug) rug(x)
        if (grid) grid()
        lines(density, col=col[1], lwd=lwd, lty=lty[1], xlim=xlim, ylim=ylim)
    }
    else {
        if (!is.factor(g)) stop("argument g must be a factor")
        counts <- table(g)
        if (any(counts == 0)){
            levels <- levels(g)
            warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
            g <- factor(g, levels=levels[counts != 0])
        }
        valid <- complete.cases(x, g)
        x <- x[valid]
        g <- g[valid]
        levels <- levels(g)
        if (is.numeric(bw) && length(bw) == 1) bw <- rep(bw, length(levels))
        if (length(adjust) == 1) adjust <- rep(adjust, length(levels))
        if (is.numeric(bw) && length(bw) != length(levels)) stop("number of entries in bw be 1 or must equal number of groups")
        if (length(adjust) != length(levels)) stop("number of entries in adjust must be 1 or must equal number of groups")
        densities <- vector(length(levels), mode="list") 
        names(adjust) <- names(densities) <- levels
        if (is.numeric(bw)) names(bw) <- levels
        for (group in levels){
            densities[[group]] <- if (method == "adaptive") 
                adaptiveKernel(x[g == group], bw=if (is.numeric(bw)) bw[group] else bw, adjust=adjust[group], ...)
                else density(x[g == group], bw=if (is.numeric(bw)) bw[group] else bw, adjust=adjust[group], kernel=kernel, ...)
            if (normalize) densities[[group]]$y <- norm(densities[[group]]$x, densities[[group]]$y)
        }
        if (missing(xlim)){
            xlim <- range(unlist(lapply(densities, function(den) range(den$x))))
        }
        if (missing(ylim)){
            ylim <- c(0, max(sapply(densities, function(den) max(den$y))))
        }
        plot(xlim, ylim, xlab=xlab, ylab=ylab, main=main, type="n")
        if (grid) grid()
        for (i in 1:length(levels)){
            lines(densities[[i]]$x, densities[[i]]$y, lty=lty[i], col=col[i], lwd=lwd)
        }
        if (show.bw){
            bws <- sapply(densities, function(den) den$bw)
            legend.values <- paste(levels, " (bw = ", format(bws), ")", sep="")
        }
        else legend.values <- levels
        if (!isFALSE(legend)) legend(legend.location, legend=legend.values, col=col[1:length(levels)], 
               lty=lty, title=legend.title, inset=0.02)
        abline(h=0, col="gray")
        if (rug){
            for (i in 1:length(levels)) rug(x[g == levels[i]], col=col[i])
        }
    }
    return(invisible(if (missing(g)) density else densities))
}

densityPlot.formula <- function(formula, data=NULL, subset, na.action=NULL, xlab, ylab, main="", legend=TRUE, ...){
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$legend <- m$xlab <- m$ylab <-m$main <- m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    if (missing(ylab)) ylab <- "Density"
    response <- attr(attr(mf, "terms"), "response")
    if (length(formula) == 3){
        legend <- applyDefaults(legend, defaults=list(location="topright", title=names(mf)[-response]), type="legend")
        if (isFALSE(legend)){
            legend.title <- ""
            legend.location <- "topleft"
        }
        else{
            legend.title <- legend$title
            legend.location <- legend$location
        }
        if (missing(xlab)) xlab <- names(mf)[response]
        g <- mf[, -response]
        densityPlot(mf[[response]], g, xlab=xlab, ylab=ylab, main=main,
                    legend=if (isFALSE(legend)) FALSE else list(title=legend.title, location=legend.location), ...)
    }
    else if (length(formula) == 2){
        if (missing(xlab)) xlab <- names(mf)
        densityPlot(mf[[1]], xlab=xlab, ylab=ylab, main=main, ...)
    }
    else stop("improper densityPlot formula")   
}
