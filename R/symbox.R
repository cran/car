# 2010-09-05: J. Fox: allow xlab argument, pass through ...
# 2013-08-19: J. Fox: remove loading of stats package
# 2018--7-27: J. Fox: automatically generate start

symbox <- function(x, ...){
	UseMethod("symbox")
}

symbox.formula <- function(formula, data=NULL, subset, na.action=NULL, ylab, ...){
	variable <- all.vars(formula)
	if (length(variable) != 1) stop("the formula must specify one variable")
	m <- match.call(expand.dots = FALSE)
	if (is.matrix(eval(m$data, parent.frame()))) 
		m$data <- as.data.frame(data)
	m$ylab <- m$... <- NULL
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")
	mf <- eval(m, parent.frame())
	if (missing(ylab)) ylab <- paste("Transformations of", variable)
	symbox(as.vector(mf[[1]]), ylab=ylab, ...)
}

symbox.default <- function(x, powers = c(-1, -0.5, 0, 0.5, 1), start, trans=bcPower, 
		xlab="Powers", ylab, ...) {
    if (!(is.vector(x) && is.numeric(x))) stop("x should be a numeric vector.")
	if (missing(ylab)) ylab <- deparse(substitute(x))
	trans.name <- deparse(substitute(trans))
	if (missing(start)){
	  if (trans.name %in% c("bcPower", "bcnPower")){
	    if ((min.x <- min(x, na.rm=TRUE)) <= 0){
	      max.x <- max(x, na.rm=TRUE)
	      start <- abs(min.x) + 0.01*(max.x - min.x)
	      warning("start set to ", format(start))
	    } else {
	      start <- 0
	    }
	  } else {
	    start <- 0
	  }
	}
  x <- x + start
  if (trans.name == "bcnPower") trans <- function(x, lambda) bcnPower(x, lambda, gamma=start)
	result <- lapply(powers, function(lambda) trans(x, lambda))
	names <- as.character(powers)
	names[powers == 0] <- "log"
	names(result) <- names
	result <- as.data.frame(scale(do.call(cbind, result)))
    boxplot(result, xlab=xlab, ylab=ylab, yaxt="n", ...)
}
