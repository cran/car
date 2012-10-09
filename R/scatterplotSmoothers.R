# Scatterplot Smoothers (J. Fox and S. Weisberg)

# Sept 17, 2012 moved from scatterplot.R to scatterplotSmoothers.R

default.arg <- function(args.list, arg, default){
    if (is.null(args.list[[arg]])) default else args.list[[arg]]
}

loessLine <- function(x, y, col, log.x, log.y, spread=FALSE, smoother.args,
               draw=TRUE) {
    lty <- default.arg(smoother.args, "lty", 1)
    lwd <- default.arg(smoother.args, "lwd", 2)
    lty.spread <- default.arg(smoother.args, "lty.spread", 2)
    lwd.spread <- default.arg(smoother.args, "lwd.spread", 1)
    span <- default.arg(smoother.args, "span", 0.5)
    family <- default.arg(smoother.args, "family", "symmetric")
    degree <- default.arg(smoother.args, "degree", 2)
    iterations <- default.arg(smoother.args, "iterations", 4)
    if (log.x) x <- log(x)
    if (log.y) y <- log(y)
    valid <- complete.cases(x, y)
    x <- x[valid]
    y <- y[valid]
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    warn <- options(warn=-1)
    on.exit(options(warn))
# mean smooth
    fit <- try(loess(y ~ x, span=span, family=family, degree=degree,
                    control=loess.control(iterations=iterations)), silent=TRUE)
    if (class(fit)[1] != "try-error"){
            if (log.x) x <- exp(x)
            y <- if (log.y) exp(fitted(fit)) else fitted(fit)
            if(draw)lines(x, y, lwd=lwd, col=col, lty=lty) else
               out <- list(x=x, y=y)
            }
    else{ options(warn)
          warning("could not fit smooth")
          return()}
# spread smooth, if requested
    if(spread) {
        res <- residuals(fit)
        pos <- res > 0
        pos.fit <- try(loess(res^2 ~ x, span=span, degree=0, family=family, subset=pos,
                        control=loess.control(iterations=1)),
                        silent=TRUE)
        neg.fit <- try(loess(res^2 ~ x, span=span, degree=0, family=family, subset=!pos,
                        control=loess.control(iterations=1)),
                        silent=TRUE)
        if(class(pos.fit)[1] != "try-error"){
            y.pos <- if (log.y) exp(fitted(fit)[pos] + sqrt(fitted(pos.fit)))
                     else fitted(fit)[pos] + sqrt(fitted(pos.fit))
            if(draw) lines(x[pos], y.pos, lwd=lwd.spread, lty=lty.spread, col=col)
                else {out$x.pos <- x[pos]
                      out$y.pos <- y.pos}
            }
        else{ options(warn)
            warning("coud not fit positive part of the spread")
            }
        if(class(neg.fit)[1] != "try-error"){
            y.neg <- if (log.y) exp(fitted(fit)[!pos] - sqrt(fitted(neg.fit)))
                     else fitted(fit)[!pos] - sqrt(fitted(neg.fit))
            if(draw) lines(x[!pos], y.neg, lwd=lwd.spread, lty=lty.spread, col=col)
                 else {out$x.neg <- x[!pos]
                      out$y.neg <- y.neg}
            }
        else {options(warn)
            warning("cound not fit negative part of the spread") }
        }
    if(!draw) return(out)
    }


gamLine <- function(x, y, col, log.x, log.y, spread=FALSE, smoother.args,
              draw=TRUE) {
    if (!require(mgcv)) stop("mgcv package missing")
    lty <- default.arg(smoother.args, "lty", 1)
    lwd <- default.arg(smoother.args, "lwd", 2)
    lty.spread <- default.arg(smoother.args, "lty.spread", 2)
    lwd.spread <- default.arg(smoother.args, "lwd.spread", 1)
    family <- default.arg(smoother.args, "family", gaussian)
    k <- default.arg(smoother.args, "k", -1)
    bs <- default.arg(smoother.args, "bs", "tp")
    if (is.character(family)) family <- eval(parse(text=family))
    link <- default.arg(smoother.args, "link", NULL)
    weights <- default.arg(smoother.args, "weights", NULL)
    spread <- spread && identical(family, gaussian) && is.null(link)
    if (log.x) x <- log(x)
    if (log.y) y <- log(y)
    valid <- complete.cases(x, y)
    x <- x[valid]
    y <- y[valid]
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    w <-if (is.null(weights)) rep(1, length(y))
    else weights[valid][ord]
    warn <- options(warn=-1)
    on.exit(options(warn))
    fit <- try(if (is.null(link)) gam(y ~ s(x, k=k, bs=bs), weights=w, family=family)
                 else gam(y ~ s(x, k=k, bs=bs), weights=w, family=family(link=link)),
                 silent=TRUE)
    if (class(fit)[1] != "try-error"){
            if (log.x) x <- exp(x)
            y <- if (log.y) exp(fitted(fit)) else fitted(fit)
            if (draw) lines(x, y, lwd=lwd, col=col, lty=lty) else
               out <- list(x=x, y=y)
            }
    else{ options(warn)
          warning("could not fit smooth")
          return()}
    if(spread) {
        res <- residuals(fit)
        pos <- res > 0
        pos.fit <- try(gam(res^2 ~ s(x, k=k, bs=bs), subset=pos), silent=TRUE)
        neg.fit <- try(gam(res^2 ~ s(x, k=k, bs=bs), subset=!pos), silent=TRUE)
        if(class(pos.fit)[1] != "try-error"){
            y.pos <- if (log.y) exp(fitted(fit)[pos] + sqrt(fitted(pos.fit)))
            else fitted(fit)[pos] + sqrt(fitted(pos.fit))
            if(draw) lines(x[pos], y.pos, lwd=lwd.spread, lty=lty.spread, col=col)
               else {out$x.pos <- x[pos]
                     out$y.pos <- y.pos}
            }
        else{ options(warn)
            warning("coud not fit positive part of the spread")
            }
        if(class(neg.fit)[1] != "try-error"){
            y.neg <- if (log.y) exp(fitted(fit)[!pos] - sqrt(fitted(neg.fit)))
            else fitted(fit)[!pos] - sqrt(fitted(neg.fit))
            if(draw) lines(x[!pos], y.neg, lwd=lwd.spread, lty=lty.spread, col=col)
               else {out$x.neg <- x[!pos]
                     out$y.neg <- y.neg}
            }
        else {options(warn)
            warning("cound not fit negative part of the spread") }
        }
    if(!draw) return(out)
    }

quantregLine <- function(x, y, col, log.x, log.y, spread=FALSE, smoother.args,
                   draw=TRUE) {
    if (!require(quantreg)) stop("quantreg package missing")
    lty <- default.arg(smoother.args, "lty", 1)
    lwd <- default.arg(smoother.args, "lwd", 2)
    lty.spread <- default.arg(smoother.args, "lty.spread", 2)
    lwd.spread <- default.arg(smoother.args, "lwd.spread", 1)
    if (log.x) x <- log(x)
    if (log.y) y <- log(y)
    lambda <- default.arg(smoother.args, "lambda", IQR(x))
    valid <- complete.cases(x, y)
    x <- x[valid]
    y <- y[valid]
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    if (!spread){
        fit <- rqss(y ~ qss(x, lambda=lambda))
        if (log.x)  x <- exp(x)
        y <-if (log.y) exp(fitted(fit)) else fitted(fit)
        if(draw) lines(x, y, lwd=lwd, col=col, lty=lty)  else
           out <- list(x=x, y=x)
    }
    else{
        fit <- rqss(y ~ qss(x, lambda=lambda))
        q1fit <- rqss(y ~ qss(x, lambda=lambda), tau=0.25)
        q3fit <- rqss(y ~ qss(x, lambda=lambda), tau=0.75)
        if (log.x) x <- exp(x)
        y <- if (log.y) exp(fitted(fit)) else fitted(fit)
        if(draw) lines(x, y, lwd=lwd, col=col, lty=lty) else
           out <- list(x=x, y=y)
        y.q1 <- if (log.y) exp(fitted(q1fit)) else fitted(q1fit)
        if(draw) lines(x, y.q1, lwd=lwd.spread, lty=lty.spread, col=col) else
           {out$x.neg <- x
            out$y.neg <- y.q1}
        y.q3 <- if (log.y) exp(fitted(q3fit)) else fitted(q3fit)
        if(draw) lines(x, y.q3, lwd=lwd.spread, lty=lty.spread, col=col) else
           {out$x.neg <- x
            out$y.neg <- y.q1}
    }
    if(!draw) return(out)
}
