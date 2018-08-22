# fancy scatterplot matrices (J. Fox)

# 2010-09-04: J. Fox: changed color choice
# 2010-09-16: fixed point color when col is length 1
# 2011-03-08: J. Fox: changed col argument
# 2012-04-18: J. Fox: fixed labels argument in scatterplotMatrix.formula()
# 2012-09-12: J. Fox: smoother now given as function
# 2012-09-19: J. Fox: restored smooth and span args for backwards compatibility
# 2013-02-08: S. Weisberg: bug-fix for showLabels with groups
# 2013-08-26: J. Fox: added use argument
# 2014-08-07: J. Fox: plot univariate distributions by group (except for histogram)
# 2014-08-17: J. Fox: report warning rather than error if not enough points in a group
#                     to compute density
# 2014-09-04: J. Fox: empty groups produce warning rather than error
# 2017-02-14: J. Fox: consolidated smooth, id, legend, and ellipse arguments
# 2017-02-17: S. Weisberg, more changes to arguments
# 2017-02-19: J. Fox: bug fixes and improvement to col argument
# 2017-04-18; S. Weisberg fixed bug in handling id=FALSE with matrix/data frame input.
# 2017-04-18; S. Weisberg changed the default for by.groups to TRUE
# 2017-04-20: S. Weisberg fixed bug with color handling
# 2017-04-20: S. Weisberg the default diagonal is now adaptiveDensity using adaptiveKernel fn
#                  diagonal argument is now a list similar to regLine and smooth
#                  changed arguments and updated man page
# 2017-05-08: S. Weisberg changed col=carPalette()
# 2017-06-22: J. Fox: eliminated extraneous code for defunct labels argument; small cleanup
# 2017-12-07: J. Fox: added fill, fill.alpha subargs to ellipse arg, suggestion of Michael Friendly.
# 2018-02-09: S. Weisberg removed the transform and family arguments from the default method
# 2018-04-02: J. Fox: warning rather than error for too few colors.
# 2018-04-12: J. Fox: clean up handling of groups arg.

scatterplotMatrix <- function(x, ...){
  UseMethod("scatterplotMatrix")
}

scatterplotMatrix.formula <- function (formula, data=NULL, subset, ...) {
  na.save <- options(na.action=na.omit)
  on.exit(options(na.save))
  na.pass <- function(dframe) dframe
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, sys.frame(sys.parent()))))
    m$data <- as.data.frame(data)
  m$id <- m$formula <- m$... <- NULL
  m$na.action <- na.pass
  m[[1]] <- as.name("model.frame")
  if (!inherits(formula, "formula") | length(formula) != 2)
    stop("invalid formula")
  rhs <- formula[[2]]
  if ("|" != deparse(rhs[[1]])){
    groups <- FALSE
  }
  else{
    groups <- TRUE
    formula <- as.character(formula)
    formula <- as.formula(sub("\\|", "+", formula))
  }
  m$formula <-formula
  if (missing(data)){
    X <- na.omit(eval(m, parent.frame()))
 #   if (is.null(labels)) labels <- gsub("X", "", row.names(X))
  }
  else{
    X <- eval(m, parent.frame())
 #   if (is.null(labels)) labels <- rownames(X)
  }
  if (!groups) scatterplotMatrix(X, ...)
  else{
    ncol<-ncol(X)
    scatterplotMatrix.default(X[, -ncol], groups=X[, ncol], ...)
  }
}


scatterplotMatrix.default <-
  function(x, smooth=TRUE, id=FALSE, legend=TRUE,
           regLine=TRUE, ellipse=FALSE,
           var.labels=colnames(x),
           diagonal=TRUE,
           plot.points=TRUE,
           groups=NULL, by.groups=TRUE,
           use=c("complete.obs", "pairwise.complete.obs"),
           col=carPalette()[-1],
           pch=1:n.groups,
           cex=par("cex"), cex.axis=par("cex.axis"),
           cex.labels=NULL, cex.main=par("cex.main"), row1attop=TRUE, ...){
  transform <- FALSE
#  family <- "bcPower"
  force(col)
#  n.groups <- if(by.groups) length(levels(groups)) else 1
  if(isFALSE(diagonal)) diagonal <- "none" else {
    diagonal.args <- applyDefaults(diagonal, defaults=list(method="adaptiveDensity"), type="diag")
    diagonal <- if(!isFALSE(diagonal.args)) diagonal.args$method
    diagonal.args$method <- NULL
  }
# regLine; use old arguments reg.line, lty and lwd
  regLine.args <- applyDefaults(regLine, defaults=list(method=lm, lty=1, lwd=2,
                                                       col=col), type="regLine")
  if(!isFALSE(regLine.args)) {
    reg.line <- regLine.args$method
    lty <- regLine.args$lty
    lwd <- regLine.args$lwd
  } else reg.line <- "none"
  # setup smoother, now including spread
  n.groups <- if(is.null(groups)) 1
    else {
      if (!is.factor(groups)) groups <- as.factor(groups)
      length(levels(groups))
    }
  smoother.args <- applyDefaults(smooth, defaults=list(smoother=loessLine,
                              spread=(n.groups)==1, col=col, lty.smooth=2, lty.spread=4), type="smooth")
  if (!isFALSE(smoother.args)) {
    # check for an argument 'var' in smoother.args.
    if(!is.null(smoother.args$var)) smoother.args$spread <- smoother.args$var
    # end change
    smoother <- smoother.args$smoother
    spread <- if(is.null(smoother.args$spread)) TRUE else smoother.args$spread
    smoother.args$spread <- smoother.args$smoother <- NULL
    if(n.groups==1) smoother.args$col <- col[1]
  }
  else smoother <- "none"
  # setup id
  id <- applyDefaults(id, defaults=list(method="mahal", n=2, cex=1, col=col, location="lr"), type="id")
  if (is.list(id) && "identify" %in% id$method) stop("interactive point identification not permitted")
  if (isFALSE(id)){
    id.n <- 0
    id.method <- "mahal"
    labels <- id.cex <- id.col <- id.location <- NULL
  }
  else{
    labels <- if(!is.null(id$labels)) id$labels else row.names(x)
    id.method <- id$method
    id.n <- id$n
    id.cex <- id$cex
    id.col <- id$col
    id.location <- id$location
  }
  if (is.null(labels)) labels <- as.character(seq(length.out=nrow(x)))
  legend <- applyDefaults(legend, defaults=list(coords=NULL), type="legend")
  if (!(isFALSE(legend) || missing(groups))){
    legend.plot <- TRUE
    legend.pos <- legend$coords
  }
  else {
    legend.plot <- FALSE
    legend.pos <- NULL
  }
  # ellipse
  ellipse <- applyDefaults(ellipse, defaults=list(levels=c(0.5, 0.95), robust=TRUE, fill=TRUE, fill.alpha=0.2), type="ellipse")
  if (isFALSE(ellipse)){
    levels <- NULL
    robust <- NULL
  }
  else{
    levels <- ellipse$levels
    robust <- ellipse$robust
    fill <- ellipse$fill
    fill.alpha <- ellipse$fill.alpha
    ellipse <- TRUE
  }
  # pre 2017 code follows
#  family <- match.arg(family)
  use <- match.arg(use)
  na.action <- if (use == "complete.obs") na.omit else na.pass
  if (!(missing(groups))){
    x <- na.action(data.frame(groups, labels, x, stringsAsFactors=FALSE))
    #      groups <- as.factor(as.character(x[, 1]))
    groups <- x$groups
#    if (!is.factor(groups)) groups <- as.factor(as.character(x[,1]))
    labels <- x[, 2]
    x <- x[, -(1:2)]
  }
  else {
    x <- na.action(data.frame(labels, x, stringsAsFactors=FALSE))
    labels <- x[, 1]
    x <- x[, -1]
    id.col <- id.col[1]
  }
  legendPlot <- function(position="topright"){
    usr <- par("usr")
    legend(position, bg="white",
           legend=levels(groups), pch=pch, col=col[1:n.groups],
           cex=cex)
  }
  do.legend <- legend.plot
####### diagonal panel functions
  # The following panel function adapted from Richard Heiberger
  panel.adaptiveDensity <- function(x, ...){
    args <- applyDefaults(diagonal.args,
        defaults=list(bw=bw.nrd0, adjust=1, kernel=dnorm, na.rm=TRUE))
    if (n.groups > 1){
      levs <- levels(groups)
      for (i in 1:n.groups){
        xx <- x[levs[i] == groups]
        dens.x <- try(adaptiveKernel(xx, adjust = args$adjust, na.rm=args$na.rm,
                              bw=args$bw, kernel=args$kernel), silent=TRUE)
        if (!inherits(dens.x, "try-error")){
          lines(dens.x$x, min(x, na.rm=TRUE) + dens.x$y *
                  diff(range(x, na.rm=TRUE))/diff(range(dens.x$y, na.rm=TRUE)), col=col[i])
        }
        else warning("cannot estimate density for group ", levs[i], "\n",
                     dens.x, "\n")
        rug(xx, col=col[i])
      }
    }
    else {
      dens.x <- adaptiveKernel(x, adjust = args$adjust, na.rm=args$na.rm,
                        bw=args$bw, kernel=args$kernel)
      lines(dens.x$x, min(x, na.rm=TRUE) + dens.x$y * diff(range(x, na.rm=TRUE))/diff(range(dens.x$y, na.rm=TRUE)), col=col[1])
      rug(x)
    }
    if (do.legend) legendPlot(position=if (is.null(legend.pos)) "topright" else legend.pos)
    do.legend <<- FALSE
  }
#
  panel.density <- function(x, ...){
    args <- applyDefaults(diagonal.args,
                          defaults=list(bw="nrd0", adjust=1, kernel="gaussian", na.rm=TRUE))
    if (n.groups > 1){
      levs <- levels(groups)
      for (i in 1:n.groups){
        xx <- x[levs[i] == groups]
        dens.x <- try(density(xx, adjust = args$adjust, na.rm=args$na.rm,
                              bw=args$bw, kernel=args$kernel), silent=TRUE)
        if (!inherits(dens.x, "try-error")){
          lines(dens.x$x, min(x, na.rm=TRUE) + dens.x$y *
                  diff(range(x, na.rm=TRUE))/diff(range(dens.x$y, na.rm=TRUE)), col=col[i])
        }
        else warning("cannot estimate density for group ", levs[i], "\n",
                     dens.x, "\n")
        rug(xx, col=col[i])
      }
    }
    else {
      dens.x <- density(x, adjust = args$adjust, na.rm=args$na.rm,
                        bw=args$bw, kernel=args$kernel)
      lines(dens.x$x, min(x, na.rm=TRUE) + dens.x$y * diff(range(x, na.rm=TRUE))/diff(range(dens.x$y, na.rm=TRUE)), col=col[1])
      rug(x)
    }
    if (do.legend) legendPlot(position=if (is.null(legend.pos)) "topright" else legend.pos)
    do.legend <<- FALSE
  }
  panel.histogram <- function(x, ...){
    par(new=TRUE)
    args <- applyDefaults(diagonal.args, defaults=list(breaks="FD"))
    h.col <- col[1]
    if (h.col == "black") h.col <- "gray"
    hist(x, main="", axes=FALSE, breaks=args$breaks, col=h.col)
    if (do.legend) legendPlot(position=if (is.null(legend.pos)) "topright" else legend.pos)
    do.legend <<- FALSE
  }
  panel.boxplot <- function(x, ...){
    b.col <- col[1:n.groups]
    b.col[b.col == "black"] <- "gray"
    par(new=TRUE)
    if (n.groups == 1) boxplot(x, axes=FALSE, main="", col=col[1])
    else boxplot(x ~ groups, axes=FALSE, main="", col=b.col)
    if (do.legend) legendPlot(position=if (is.null(legend.pos)) "topright" else legend.pos)
    do.legend <<- FALSE
  }
  # The following panel function adapted from Richard Heiberger
  panel.oned <- function(x, ...) {
    range <- range(x, na.rm=TRUE)
    delta <- diff(range)/50
    y <- mean(range)
    if (n.groups == 1) segments(x - delta, x, x + delta, x, col = col[1])
    else {
      segments(x - delta, x, x + delta, x, col = col[as.numeric(groups)])
    }
    if (do.legend) legendPlot(position=if (is.null(legend.pos)) "bottomright" else legend.pos)
    do.legend <<- FALSE
  }
  panel.qqplot <- function(x, ...){
    par(new=TRUE)
    if (n.groups == 1) qqnorm(x, axes=FALSE, xlab="", ylab="", main="", col=col[1])
    else qqnorm(x, axes=FALSE, xlab="", ylab="", main="", col=col[as.numeric(groups)])
    qqline(x, col=col[1])
    if (do.legend) legendPlot(position=if (is.null(legend.pos)) "bottomright" else legend.pos)
    do.legend <<- FALSE
  }
  panel.blank <- function(x, ...){
    if (do.legend) legendPlot(if (is.null(legend.pos)) "topright" else legend.pos)
    do.legend <<- FALSE
  }
  which.fn <- match(diagonal,
                    c("adaptiveDensity", "density", "boxplot", "histogram", "oned", "qqplot", "none"))
  if(is.na(which.fn)) stop("incorrect name for the diagonal argument, see ?scatterplotMatrix")
  diag <- list(panel.adaptiveDensity, panel.density, panel.boxplot, panel.histogram, panel.oned,
               panel.qqplot, panel.blank)[[which.fn]]
  groups <- as.factor(if(missing(groups)) rep(1, length(x[, 1])) else groups)
  counts <- table(groups)
  if (any(counts == 0)){
    levels <- levels(groups)
    warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
    groups <- factor(groups, levels=levels[counts > 0])
  }
#  n.groups <- length(levels(groups))
  if (n.groups > length(col)) {
    warning("number of groups exceeds number of available colors\n  colors are recycled")
    col <- rep(col, n.groups)
  }
  if (length(col) == 1) col <- rep(col, 3)
  labs <- labels
  pairs(x, labels=var.labels,
        cex.axis=cex.axis, cex.main=cex.main, cex.labels=cex.labels, cex=cex,
        diag.panel=diag, row1attop = row1attop,
        panel=function(x, y, ...){
          for (i in 1:n.groups){
            subs <- groups == levels(groups)[i]
            if (plot.points) points(x[subs], y[subs], pch=pch[i], col=col[if (n.groups == 1) 1 else i], cex=cex)
            if (by.groups){
              if (is.function(smoother)) smoother(x[subs], y[subs], col=smoother.args$col[i],
                                                  log.x=FALSE, log.y=FALSE, spread=spread, smoother.args=smoother.args)
              if (is.function(reg.line)) regLine(reg.line(y[subs] ~ x[subs]), lty=lty, lwd=lwd, col=regLine.args$col[i])
              if (ellipse) dataEllipse(x[subs], y[subs], plot.points=FALSE,
                                       levels=levels, col=col[i], robust=robust, lwd=1,
                                       fill=fill, fill.alpha=fill.alpha)
              showLabels(x[subs], y[subs], labs[subs], method=id.method,
                         n=id.n, col=col[i], cex=id.cex, location=id.location,
                         all=list(labels=labs, subs=subs))
            }
          }
          if (!by.groups){
            if (is.function(reg.line)) abline(reg.line(y ~ x), lty=lty, lwd=lwd, col=regLine.args$col[1])
            if (is.function(smoother)) smoother(x, y, col=col[1],
                                                log.x=FALSE, log.y=FALSE, spread=spread, smoother.args=smoother.args)
            if (ellipse) dataEllipse(x, y, plot.points=FALSE, levels=levels, col=smoother.args$col,
                                     robust=robust, lwd=1, fill=fill, fill.alpha=fill.alpha)
            showLabels(x, y, labs, method=id.method,
                       n=id.n, col=id.col, location=id.location, cex=id.cex)
          }
        }, ...
  )
}

spm <- function(x, ...){
  scatterplotMatrix(x, ...)
}
