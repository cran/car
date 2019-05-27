# fancy scatterplots  (J. Fox)

# 2010-09-05: J. Fox: changed color choice
# 2010-09-16: fixed point color when col is length 1
# 2010-12-19: J. Fox: added argument legend.coords to place legend.
# 2011-01-15: J. Fox: If x is a factor, calls Boxplot()
# 2011-03-08: J. Fox: changed col argument
# 2012-04-18: J. Fox: fixed labels argument in scatterplot.formula().
# 2012-04-24: J. Fox: further fix to labels
# 2012-09-12: J. Fox: modified treatment of smoother; added loessLine(), gamLine(), quantregLine().
# 2012-09-17: S. Weisberg:  smoothers moved to scatterplotSmoothers.R, defaults changed
# 2012-09-19: J. Fox: restored smooth and span arguments for backwards compatibility
# 2013-02-07: S. Weisberg:  modifed call to showLabels to work correctly with groups
# 2014-09-04: J. Fox: empty groups produce warning rather than error
# 2015-07-17: J. Fox: improved above-plot legends.
# 2015-08-05: J. Fox: fixed sp()
# 2017-01-09: J. Fox: consolidated many arguments into id, smooth, and legend;
# 2017-02-09: J. Fox: consolidated ellipse arguments; small fixes.
# 2017-02-17: S. Weisberg:  removed many arguments that can be passed via other argument or ...
# 2017-02-22: J. Fox: improvement to col argument
# 2017-02-28: S. Weisberg:  showLabels bug-fix.
# 2016-02-28: S. Weisberg:  added cex arg to the legend
# 2017-04-14: S. Weisberg:  changed default colors so points and corresponding lines always have same color
# 2017-05-08: S. Weisberg changed col=carPalette()
# 2017-11-30: substitute carPalette() for palette(). J. Fox
# 2017-12-07: J. Fox: added fill, fill.alpha subargs to ellipse arg, suggestion of Michael Friendly.
# 2018-03-23: J. Fox: fix ellipses when log-axes used by groups; fix interactive point identification by groups.
# 2018-04-02: J. Fox: warning rather than error for too few colors.
# 2018-04-12: J. Fox: fixed error produced when groups not a factor, reported by Alexandre Courtiol.
# 2018-05-19: J. Fox: fixed bug when legend=FALSE, reported by Castor Guisande.
# 2018-06-25: S. Weisberg  made the argument 'var' an alias of 'spread'
# 2019-01-15: J. Fox: make scatterplot.formula() more robust

reg <- function(reg.line, x, y, col, lwd, lty, log.x, log.y){
  if(log.x) x <- log(x)
  if(log.y) y <- log(y)
  mod <- reg.line(y ~ x)
  y.hat <- fitted.values(mod)
  x <- model.matrix(mod)[, 2]
  min <- which.min(x)
  max <- which.max(x)
  if (!log.x){
    x1 <- x[min]
    x2 <- x[max]
  }
  else {
    x1 <- exp(x[min])
    x2 <- exp(x[max])
  }
  if (!log.y){
    y1 <- y.hat[min]
    y2 <- y.hat[max]
  }
  else {
    y1 <- exp(y.hat[min])
    y2 <- exp(y.hat[max])
  }
  lines(c(x1, x2), c(y1, y2), lwd=lwd, col=col, lty=lty)
}

find.legend.columns <- function(n, target=min(4, n)){
  rem <- n %% target
  if (rem != 0 && rem < target/2) target <- target - 1
  target
}

scatterplot <- function(x, ...){
  UseMethod("scatterplot", x)
}

scatterplot.formula <- function (formula, data, subset, xlab, ylab,
                                 id=FALSE, legend=TRUE, ...) {
  na.save <- options(na.action=na.omit)
  on.exit(options(na.save))
  na.pass <- function(dframe) dframe
  id <- if (is.logical(id)){
    if (isTRUE(id)) list()
    else FALSE
  }
  else as.list(id)
  legend <- applyDefaults(legend, defaults=list(), type="legend")
  m <- match.call(expand.dots=FALSE)
  if (is.matrix(eval(m$data, sys.frame(sys.parent()))))
    m$data <- as.data.frame(data)
  m$na.action <- na.pass
  m$legend <- m$id <- m$xlab <- m$ylab <- m$... <- NULL
  m[[1]] <- as.name("model.frame")
  if (!inherits(formula, "formula") | length(formula) != 3)
    stop("invalid formula")
  formula <- as.character(c(formula))
  formula <- as.formula(sub("\\|", "+", formula))
  m$formula <- formula
  if (missing(data)){
    X <- na.omit(eval(m, parent.frame()))
    if (!isFALSE(id) && is.null(id$labels)) id$labels <- gsub("X", "", row.names(X))
    if (is.factor(X[, 2]) && !is.list(id))  id <- list(labels=gsub("X", "", row.names(X)))
  }
  else{
    X <- eval(m, parent.frame())
    if (!isFALSE(id) && is.null(id$labels)) id$labels <- row.names(X)
    if (is.factor(X[, 2]) && !is.list(id))  id <- list(labels=row.names(X))
  }
  names <- names(X)
  if (missing(xlab)) xlab <- names[2]
  if (missing(ylab)) ylab <- names[1]
  X[, 1] <- as.vector(X[, 1])
  if (!is.factor(X[, 2])) X[, 2] <- as.vector(X[, 2])
  if (ncol(X) == 2) scatterplot(X[,2], X[,1], xlab=xlab, ylab=ylab,
                                id=id, ...)
  else {
    if (!isFALSE(legend)){
      if (is.null(legend$title)) legend$title <- names[3]
    }
    scatterplot(X[,2], X[,1], groups=X[,3], xlab=xlab, ylab=ylab,
                legend=legend, id=id, ...)
  }
}

scatterplot.default <- function(x, y, boxplots=if (by.groups) "" else "xy",
                                regLine=TRUE, legend=TRUE, id=FALSE, ellipse=FALSE, grid=TRUE,
                                smooth=TRUE,
                                groups, by.groups=!missing(groups),
                                xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
                                log="", jitter=list(), cex=par("cex"),
                                col=carPalette()[-1], pch=1:n.groups,
                                reset.par=TRUE, ...){
  force(col)
  id <- applyDefaults(id, defaults=list(method="mahal", n=2, cex=1, col=carPalette()[-1], location="lr"), type="id")
  legend <- applyDefaults(legend, defaults=list(title=deparse(substitute(groups)), inset=0.02, cex=1))
  legend.plot <- !(isFALSE(legend) || missing(groups))
  if (legend.plot){
    legend.title <- legend$title
    legend.cex <- legend$cex
  }
  if (isFALSE(id)){
    id.n <- 0
    id.method <- "mahal"
    labels <- id.cex <- id.col <- id.location <- NULL
  }
  else{
    labels <- id$labels
    id.method <- id$method
    id.n <- if ("identify" %in% id.method) Inf else id$n
    id.cex <- id$cex
    id.col <- if (by.groups) id$col else id$col[1]
    id.location <- id$location
  }
  smoother.args <- applyDefaults(smooth, defaults=list(smoother=loessLine, spread=!by.groups, lty.smooth=2, lty.spread=4), type="smooth")
  if (!isFALSE(smoother.args)) {
# check for an argument 'var' in smoother.args.
    if(!is.null(smoother.args$var)) smoother.args$spread <- smoother.args$var
# end change
    smoother <- smoother.args$smoother
    spread <- if(is.null(smoother.args$spread)) TRUE else smoother.args$spread
    smoother.args$smoother <- NULL
  }
  else smoother <- "none"
  ellipse.args <- applyDefaults(ellipse, defaults=list(levels=c(.5, .95), robust=TRUE, fill=TRUE, fill.alpha=0.2, type="ellipse"))
  if (!is.logical(ellipse)) ellipse <- TRUE
  if (!isFALSE(ellipse.args)){
    levels <- ellipse.args$levels
    robust <- ellipse.args$robust
    fill <- ellipse.args$fill
    fill.alpha <- ellipse.args$fill.alpha
  }
  n.groups <- if (by.groups) {
    if (!is.factor(groups)) groups <- as.factor(groups)
    length(levels(groups))
    }
    else 1
  regLine.args <- applyDefaults(regLine, defaults=list(method=lm, lty=1, lwd=2,
                                                       col=rep(col, n.groups), type="regLine"))
  if(!isFALSE(regLine.args)) {
    if(length(regLine.args$col) < n.groups){
      regLine.args$col <- rep(regLine.args$col, n.groups)
    }
  }
  logged <- function(axis=c("x", "y")){
    axis <- match.arg(axis)
    0 != length(grep(axis, log))
  }
  hbox <- function(x){
    if (logged("x")){
      log.x <- "x"
      .x <- log(x)
    }
    else {
      log.x <- ""
      .x <- x
    }
    plot(x, seq(0, 1, length=length(x)), type="n", axes=FALSE, xlab="", ylab="", log=log.x)
    res <- boxplot.stats(.x, coef = 1.5, do.conf=FALSE)
    if (logged("x")){
      res$stats <- exp(res$stats)
      if (!is.null(res$out)) res$out <- exp(res$out)
    }
    LW <- res$stats[1]
    Q1 <- res$stats[2]
    M <- res$stats[3]
    Q3 <- res$stats[4]
    UW <- res$stats[5]
    lines(c(Q1, Q1, Q3, Q3, Q1), c(0, 1, 1, 0, 0))
    lines(c(M, M), c(0, 1))
    lines(c(LW, Q1), c(.5, .5))
    lines(c(Q3, UW), c(.5, .5))
    if (!is.null(res$out)) points(res$out, rep(.5, length(res$out)), cex=cex)
  }
  vbox <- function(y){
    if (logged("y")){
      log.y <- "y"
      .y <- log(y)
    }
    else {
      log.y <- ""
      .y <- y
    }
    plot(seq(0, 1, length=length(y)), y, type="n", axes=FALSE, xlab="", ylab="", log=log.y)
    res <- boxplot.stats(.y, coef = 1.5, do.conf=FALSE)
    if (logged("y")){
      res$stats <- exp(res$stats)
      if (!is.null(res$out)) res$out <- exp(res$out)
    }
    LW <- res$stats[1]
    Q1 <- res$stats[2]
    M <- res$stats[3]
    Q3 <- res$stats[4]
    UW <- res$stats[5]
    lines(c(0, 1, 1, 0, 0), c(Q1, Q1, Q3, Q3, Q1))
    lines(c(0, 1), c(M, M))
    lines(c(.5, .5), c(LW, Q1))
    lines(c(.5, .5), c(Q3, UW))
    if (!is.null(res$out)) points(rep(.5, length(res$out)), res$out, cex=cex)
  }
#  force(by.groups)
  id <- as.list(id)
  if (is.null(labels)){
    labels <- if (is.null(names(y)))
      seq(along=y)
    else names(y)
  }
  if (length(labels) != length(y)) stop("labels argument is the wrong length")
  if (is.factor(x)) {
    if (!(id.method %in% c("y", "identify", "none"))) id.method <- "y"
    return(Boxplot(y, x, id.method="y", labels=labels, xlab=xlab, ylab=ylab))
  }
  mar <- par("mar")
  mfcol <- par("mfcol")
  if (reset.par) on.exit(par(mar=mar, mfcol=mfcol))
  if( FALSE == boxplots) boxplots <- ""
  if (!missing(groups)){
    data <- na.omit(data.frame(groups, x, y, labels, stringsAsFactors=FALSE))
    groups <- data[ , 1]
#    if (!is.factor(groups)) groups <- as.factor(groups)
    .x <- data[,2]
    .y <- data[,3]
    labels <- data[,4]
    legend.columns <- if (legend.plot) legend$columns else 0
    top <- if (legend.plot && is.null(legend$coords)){
      if (is.null(legend.columns)) legend.columns <- find.legend.columns(nlevels(groups))
      4 + ceiling(nlevels(groups))/legend.columns
    }
    else mar[3]
    if (legend.plot && !is.null(legend$coords) && is.null(legend.columns)){
      legend.columns <- 1
    }
  }
  else {
    .x <- x
    .y <- y
    top <- mar[3]
    groups <- factor(rep(1, length(.x)))
  }
  xbox <- length(grep("x", boxplots)) > 0
  ybox <- length(grep("y", boxplots)) > 0
  if (xbox && ybox)
    layout(matrix(c(1, 0, 3, 2), 2, 2),
           widths = c(5, 95),
           heights= c(95, 5))
  else if (ybox)
    layout(matrix(c(1, 2),1, 2),
           widths = c(5, 95),
           heights= 100)
  else if (xbox)
    layout(matrix(c(2, 1), 2, 1),
           widths = 100,
           heights= c(95, 5))
  else layout (matrix(1, 1, 1),
               widths=100, heights=100)
  par(mar=c(mar[1], 0, top, 0))
  if (ybox > 0) vbox(.y)
  par(mar=c(0, mar[2], 0, mar[4]))
  if (xbox > 0) hbox(.x)
  par(mar=c(mar[1:2], top, mar[4]))
  plot(.x, .y, xlab=xlab, ylab=ylab, log=log, cex=cex,
       type="n", ...)
  if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
  if (n.groups > length(col)) {
    warning("number of groups exceeds number of available colors\n  colors are recycled")
    col <- rep(col, n.groups)
  }
  if (length(col) == 1) col <- rep(col, 3)
  indices <- NULL
  range.x <- if (logged("x")) range(log(.x), na.rm=TRUE) else range(.x, na.rm=TRUE)
  counts <- table(groups)
  if (any(counts == 0)){
    levels <- levels(groups)
    warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
  }
  for (i in 1:n.groups){
    if (counts[i] == 0) next
    subs <- groups == levels(groups)[i]
    points(if (is.null(jitter$x) || jitter$x == 0) .x[subs] else jitter(.x[subs], factor=jitter$x),
           if (is.null(jitter$y) || jitter$y == 0) .y[subs] else jitter(.y[subs], factor=jitter$y),
           pch=pch[i], col=col[if (n.groups == 1) 1 else i], cex=cex)
    if (by.groups){
      if (is.function(smoother))
        smoother(.x[subs], .y[subs], col=col[i],
                 log.x=logged("x"), log.y=logged("y"), spread=spread, smoother.args=smoother.args)
      if (!isFALSE(regLine.args)) reg(regLine.args$method, .x[subs], .y[subs], lty=regLine.args$lty,
                                      lwd=regLine.args$lwd, log.x=logged("x"), log.y=logged("y"), col=
                                        regLine.args$col[i])
      if (ellipse) {
        X <- na.omit(data.frame(x=.x[subs], y=.y[subs]))
        if (logged("x")) X$x <- log(X$x)
        if (logged("y")) X$y <- log(X$y)
        with(X, dataEllipse(x, y, plot.points=FALSE, lwd=1, log=log,
                            levels=levels, col=col[i], robust=robust,
                            fill=fill, fill.alpha=fill.alpha))
      }
      if (id.method[1] != "identify")
        indices <- c(indices,
                     showLabels(.x[subs], .y[subs], labels=labels[subs], method=id.method,
                                n=id.n, cex=id.cex, col=col[i], location=id.location,
                                all=list(labels=labels, subs=subs)))
    }}
  if (!by.groups){
    if (is.function(smoother)) smoother(.x, .y, col=col[1],
                                        log.x=logged("x"), log.y=logged("y"), spread, smoother.args=smoother.args)
    if (!isFALSE(regLine.args)) reg(regLine.args$method, .x, .y, lty=regLine.args$lty, lwd=regLine.args$lwd,
                                    log.x=logged("x"), log.y=logged("y"), col=regLine.args$col[1])
    if (ellipse) {
      X <- na.omit(data.frame(x=.x, y=.y))
      if (logged("x")) X$x <- log(X$x)
      if (logged("y")) X$y <- log(X$y)
      with(X, dataEllipse(x, y, plot.points=FALSE, lwd=1, log=log, levels=levels, col=col[1],
                          robust=robust, fill=fill, fill.alpha=fill.alpha))
    }
    if (id.method[1] != "identify") indices <- showLabels(
      .x, .y, labels=labels,
      method=id.method, n=id.n, cex=id.cex, col=id.col, location=id.location)
  }
  if (legend.plot) {
    xpd <- par(xpd=TRUE)
    on.exit(par(xpd=xpd), add=TRUE)
    usr <- par("usr")
    legend.coords <- if (is.null(legend$coords)){
      legend.x <- if (logged("x")) 10^(usr[1]) else usr[1]
      legend.y <- if (logged("y")) 10^(usr[4] + 1.2*top*strheight("x")) else usr[4] + 1.2*top*strheight("x")
      list(x=legend.x, y=legend.y)
    }
    else legend$coords
    legend(legend.coords, legend=levels(groups)[counts > 0],
           pch=pch[counts > 0], col=col[1:n.groups][counts > 0], cex=legend.cex, #pt.cex=cex, cex=cex, #cex=cex.lab,
           title=legend.title, bg="white", ncol=legend.columns, inset=legend$inset)
  }
  if (id.method[1] == "identify") indices <- showLabels(.x, .y, labels,
                    method=id.method, n=length(.x), cex=id.cex, col="black", id.location=id.location)
  if (is.null(indices)) invisible(indices) else if (is.numeric(indices)) sort(indices) else indices
}

sp <- function(x, ...)  UseMethod("scatterplot", x)

