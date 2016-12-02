# 2016-11-25: J. Fox: use pure-R code, removed compiled code.

adaptiveKernel <- function(x, kernel=dnorm, bw=bw.nrd0, adjust=1.0, n=500, from, to, cut=3, na.rm=TRUE){
    varname <- deparse(substitute(x))
    if (na.rm) x <- na.omit(x)
    if (!is.numeric(bw)) bw <- bw(x)
    bw <- adjust*bw
    if (missing(from)) from <- min(x) - cut*bw
    if (missing(to)) to <- max(x) + cut*bw
    x0 <- seq(from, to, length=n)
    n.1 <- length(x)
    p <- rep(0, n)
    initialp.x0 <- rep(0, n)
    fac <- 1/(n.1*bw)
    for (i in 1:n) initialp.x0[i] <- fac * sum(kernel((x - x0[i])/bw))
    initialp <- rep(0, n.1)
    for (i in 1:n.1) initialp[i] <- initialp.x0[which.min(abs(x[i] - x0))]
    pbar <- exp((1/n.1)*sum(log(initialp)))
    f <- (initialp/pbar)^-0.5
    for (i in 1:n) p[i] <- fac * sum((1/f)*kernel((x - x0[i])/(f*bw)))
    result <- list(x=x0, y=p, n=n, bw=bw*adjust, call=match.call(), data.name=varname, has.na=FALSE)
    class(result) <- "density"
    result
}
