# Box-Cox power transformations, with automatic start (J. Fox)

# last modified 1 June 07 by J. Fox (after a suggestion by Henric Nilsson)

bc<-function(x, p, ...) box.cox(x, p, ...)

box.cox<-function(x, p, start=0){
    nx <- length(x)
    np <- length(p)
    min<-min(x, na.rm=TRUE)
    s<-if (missing(start) & (min <= 0)) {
        IQR <- diff(quantile(x, c(.25,.75), na.rm=TRUE))
        if (IQR <= 10*.Machine$double.eps) 
            stop("First and third quartile are equal.\nAutomatic start cannot be computed.")
        nice(-min + 0.05*IQR, "up")
        }
        else start
    if (missing(start) & s != 0) 
        warning(paste("start = ", s, "added to data prior to transformation"))
    x <- x + s
    p.ln.x <- outer(x, p, function(x, p) p*log(x))
    xp <- expm1(p.ln.x)/matrix(p, nx, np, byrow=TRUE)
    sel0 <- abs(p.ln.x) < .Machine$double.eps
    sel0[is.na(sel0)] <- TRUE
    if(any(sel0)) xp[sel0] <- log(matrix(x, nx, np)[sel0])
    if (np > 1) colnames(xp) <- sapply(p, format)
    else xp <- as.vector(xp)
    xp
    }

