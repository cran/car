
# Utility functions (J. Fox)

# 16 March 2010 changed 'vars' argument to 'terms'
# 28 June 2010 added df.terms.surveg and model.matrix.survreg
# 15 November 2010 added squeezeBlanks
# 21 January 2011 added functions to support mixed models
# 2012-04-08 added exists.method
# 2012-06-23: added call to globalVariables(). John
# 2012-12-10: added .carEnv to avoid warnings in R > 2.16.0
# 2013-06020: added .merMod methods to df.residual() and has_intercept(). John
# 2014-05-16: added .multinom method for has_intercept(). John
# 2014-08-19: added package.installed() function, unexported. John
# 2014-11-02: termsToMf fixed, Sandy
# 2015-01-13: fixed model.matrix.lme() to work with model with formula as object. John
# 2015-01-27: .carEnv now lives in the global environment. John
# 2015-09-04: added model.matrix.coxme() and alias.coxme(). John
# 2015-09-11: added some support for VGAM::vglm objects. John
# 2017-02-10: added isFALSE() and applyDefaults() to support plotting functions. John
# 2017-04-14: added carPalette to the end of this file.  May require further work later.
# 2017-11-30: new version of carPalette(). John
# 2017-11-30: substitute carPalette() for palette(). J. Fox
# 2017-12-28: rewrote termsToMf used by residualPlots.  It didn't work right.  SW
# 2018-01-15: df.terms.multinom() now works with response matrix. JF
# 2018-05-23: make model.matrix.lme() more bullet proof, following report by Peter Grossmann. JF
# 2018-11-07: added combineLists(). JF
# 2019-01-02: added na.action.merMod(), removed df.residual.merMod(). JF
# 2019-10-24: include colorblind palette in carPalette(). JF
# 2019-11-14: change class(x) == "y" to inherits(x, "y")
# 2020-02-17: added matchFun() as a replacement for match.fun
# 2020-10-19: added envelope() for plotting confidence/variance envelopes. JF
# 2020-12-03: added getVcov to interpret vcov. argument as matrix or function and return an error otherwise
# 2020-12-18: getVcov() also able to return objects coercible to a matrix such as Matrix objects. JF
# 2021-04-08: added getModelData(), not explorted. JF
# 2024-05-14: format.perc() -> format_perc(), has.intercept -> has_intercept(). J. Fox
# 2024-09-20: added model.matrix.lme()

#if (getRversion() >= "2.15.1") globalVariables(c(".boot.sample", ".boot.indices"))

.carEnv <- new.env(parent=globalenv())

# function to find "nice" numbers

nice <- function(x, direction=c("round", "down", "up"), lead.digits=1){
	direction <- match.arg(direction)
	if (length(x) > 1) return(sapply(x, nice, direction=direction, lead.digits=lead.digits))
	if (x == 0) return(0)
	power.10 <- floor(log(abs(x),10))
	if (lead.digits > 1) power.10 <- power.10 - lead.digits + 1
	lead.digit <- switch(direction,
		round=round(abs(x)/10^power.10),
		down=floor(abs(x)/10^power.10),
		up=ceiling(abs(x)/10^power.10))
	sign(x)*lead.digit*10^power.10
}

has_intercept <- function (model, ...) {
	UseMethod("has_intercept")
}

has_intercept.default <- function(model, ...) any(names(coefficients(model))=="(Intercept)")

has_intercept.multinom <- function(model, ...) {
  nms <- names(coef(model))
  any(grepl("\\(Intercept\\)", nms))
}

term.names <- function (model, ...) {
	UseMethod("term.names")
}

term.names.default <- function (model, ...) {
	term.names <- labels(terms(model))
	if (has_intercept(model)) c("(Intercept)", term.names)
	else term.names
}

predictor.names <- function(model, ...) {
	UseMethod("predictor.names")
}

predictor.names.default <- function(model, ...){
	predictors <- attr(terms(model), "variables")
	as.character(predictors[3:length(predictors)])
}

responseName <- function (model, ...) {
	UseMethod("responseName")
}

responseName.default <- function (model, ...) deparse(attr(terms(model), "variables")[[2]])

response <- function(model, ...) {
	UseMethod("response")
}

response.default <- function (model, ...) model.response(model.frame(model))

is.aliased <- function(model){
	!is.null(alias(model)$Complete)
}

df.terms <- function(model, term, ...){
	UseMethod("df.terms")
}

df.terms.default <- function(model, term, ...){
	if (is.aliased(model)) stop("Model has aliased term(s); df ambiguous.")
	if (!missing(term) && 1 == length(term)){
		assign <- attr(model.matrix(model), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term)) stop(paste(term, "is not in the model."))
		sum(assign == which.term)
	}
	else {
		terms <- if (missing(term)) labels(terms(model)) else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model, term))
		names(result) <- terms
		result
	}
}

df.terms.multinom <- function (model, term, ...){
	nlev <- if (is.null(model$lev)) ncol(model.response(model.frame(model))) else length(model$lev)
	if (!missing(term) && 1 == length(term)) {
		assign <- attr(model.matrix(model), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term))
			stop(paste(term, "is not in the model."))
		sum(assign == which.term) * (nlev - 1)
	}
	else {
		terms <- if (missing(term))
				labels(terms(model))
			else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model,
					term))
		names(result) <- terms
		result
	}
}

df.terms.polr <- function (model, term, ...){
	if (!missing(term) && 1 == length(term)) {
		assign <- attr(model.matrix(model), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term))
			stop(paste(term, "is not in the model."))
		sum(assign == which.term)
	}
	else {
		terms <- if (missing(term))
				labels(terms(model))
			else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model,
					term))
		names(result) <- terms
		result
	}
}

df.terms.survreg <- function(model, term, ...){
	if (is.aliased(model)) stop("Model has aliased term(s); df ambiguous.")
	if (!missing(term) && 1 == length(term)){
		assign <- attr(model.matrix(model, data=model.frame(model)), "assign")
		which.term <- which(term == labels(terms(model)))
		if (0 == length(which.term)) stop(paste(term, "is not in the model."))
		sum(assign == which.term)
	}
	else {
		terms <- if (missing(term)) labels(terms(model)) else term
		result <- numeric(0)
		for (term in terms) result <- c(result, Recall(model, term))
		names(result) <- terms
		result
	}
}

model.matrix.survreg <- function(object, ...) model.matrix.default(object, model.frame(object))

mfrow <- function(n, max.plots=0){
	# number of rows and columns for array of n plots
	if (max.plots != 0 && n > max.plots)
		stop(paste("number of plots =",n," exceeds maximum =", max.plots))
	rows <- round(sqrt(n))
	cols <- ceiling(n/rows)
	c(rows, cols)
}

inv <- function(x) solve(x)

coefnames2bs <- function(g, para.names, parameterPrefix="b"){
	metas <- c("(", ")", "[", "]", "{", "}", ".", "*", "+", "^", "$", ":", "|")
	metas2 <- paste("\\", metas, sep="")
	metas3 <- paste("\\\\", metas, sep="")
	for (i in seq(along=metas))
		para.names <- gsub(metas2[i], metas3[i], para.names) # fix up metacharacters
	para.order <- order(nchar(para.names), decreasing=TRUE)
	para.names <- para.names[para.order] # avoid partial-name substitution
	std.names <- if ("(Intercept)" %in% para.names)
				paste(parameterPrefix, 0:(length(para.names) - 1), sep = "")
			else paste(parameterPrefix, 1:length(para.names), sep = "")
	std.names.ordered <- std.names[para.order]
	for (i in seq(along=para.names)){
		g <- gsub(para.names[i], std.names.ordered[i], g)
	}
	list(g=g, std.names=std.names)
}


showLabelsScatter <- function(x, y, labels, id.var = NULL,
	id.method = c("mahal", "identify", "none"), log="", id.cex=.75, id.n=3, id.col=carPalette()[1],
	range.x=range(.x), show=TRUE) {
	id.method <- match.arg(id.method)
	if (id.method == "none" || id.n == 0 || !show) return(invisible(NULL))
	if(id.n > 0L) {
		if (missing(labels))
			labels <- if (!is.null(id.var)) names(id.var)
				else as.character(seq(along=x))
		getPoints <- function(z) {
			names(z) <- labels
			iid <- seq(length=id.n)
			zs <- z[order(-z)[iid]]
			match(names(zs), labels)
		}
		logged <- function(axis=c("x", "y")){
			axis <- match.arg(axis)
			0 != length(grep(axis, log))
		}
		valid <- complete.cases(x, y)
		x <- x[valid]
		y <- y[valid]
		labels <- labels[valid]
		if (length(id.var) == length(valid))
			id.var <- id.var[valid]
		.x <- if (logged("x")) log(x) else x
		.y <- if (logged("y")) log(y) else y
		ind <- if (!is.null(id.var)) {
				if (length(id.var) == length(x)) order(-abs(id.var))[1L:id.n]
				else if(is.character(id.var)) match(id.var, labels) else id.var
			}
			else switch(id.method,
					x = getPoints(abs(.x - mean(.x))),
					y = getPoints(abs(.y - mean(.y))),
					xy = union(getPoints(abs(.x - mean(.x))),
						getPoints(abs(.y - mean(.y)))),
					mahal= getPoints(rowSums(qr.Q(qr(cbind(1, .x, .y))) ^ 2)))
		ind <- na.omit(ind)
		if (length(ind) == 0) return(invisible(NULL))
		labpos <- c(4, 2)[1 + as.numeric(.x > mean(range.x))]
		text(x[ind], y[ind], labels[ind], cex = id.cex, xpd = TRUE,
			pos = labpos[ind], offset = 0.25, col=id.col)
		return(labels[ind])
	}
}


#  outerLegend, written by S. Weisberg Feb 2010
#  outerLegend function
#  puts a legend in the margin, either at the upper left (margin = 3)
#  the default or upper right side otherwise
#  all the args from legend are used except for x, y, and xpd which are
#  set in the function.
#  offset is a fraction of the plot width or height to locate the legend
outerLegend <- function(..., margin=3, offset=0, adjust=FALSE){
   lims <- par("usr")
   if (margin == 3) {
      x0 <- lims[1] + offset*(lims[2]-lims[1])
      y0 <- lims[4] }
   else {
      x0 <- lims[2] + offset*(lims[2]-lims[1])
      y0 <- lims[4]
   }
   leg <- legend(x0, y0, ... , xpd=TRUE, plot=FALSE)
   if (margin == 3) {
      y0 <- y0 + leg$rect$h
      if(adjust == TRUE) x0 <- x0 - leg$text$x[1]
   }
   legend(x0, y0, ... , xpd=TRUE)
   }

# added by J. Fox 18 Nov 2010

squeezeBlanks <- function(text){
	gsub(" *", "",  text)
}

# added by J. Fox 21 Jan 2011 to support mixed models

df.residual.mer <- function(object, ...) NULL

# df.residual.merMod <- function(object, ...) NULL # no longer needed, now supplied by lme4

df.residual.lme <- function(object, ...) Inf

has_intercept.mer <- function(model, ...){
	any(names(fixef(model))=="(Intercept)")
}

has_intercept.merMod <- function(model, ...){
    any(names(fixef(model))=="(Intercept)")
}

model.matrix.lme <- function(object, ...){
  data <- object$data
  if (is.null(data)){
    NextMethod(formula(object),  data=eval(object$call$data),
               contrasts.arg=object$contrasts)
  } else {
    NextMethod(formula(object), data=data, 
               contrasts.arg=object$contrasts)
  }
} 


# added by J. Fox 2019-01-02:

na.action.merMod <- function(object, ...){
  nms <- names(attributes(model.frame(object)))
  if ("na.action" %in% nms) attributes(model.frame(object))$na.action
  else {
    na.action <- integer(0)
    class(na.action) <- options("na.action")
    na.action
  }
}

# added by J. Fox 2012-04-08 to use in deltaMethod.default()

exists.method <- function(generic, object, default=TRUE, strict=FALSE){
	classes <- class(object)
	if (default) classes <- c(classes, "default")
	if (strict) classes <- classes[1]
	any(paste(generic, ".", classes, sep="") %in%
					as.character(methods(generic)))
}


# Used by marginalModelPlots, residualPlots added 2012-09-24
plotArrayLegend <- function(
      location=c("top", "none", "separate"),
      items, col.items, lty.items, lwd.items, title="legend",
      pch=1:length(items)) {
   if(location== "none") return()
   n <- length(items)
   if(location == "top" ) { # add legend
      usr <- par("usr")
      coords <-list(x=usr[1], y=usr[3])
      leg <- legend( coords, items,
                col=col.items, pch=pch,
                bty="n", cex=1, xpd=NA, plot=FALSE)
      coords <- list(x = usr[1], y=usr[4] + leg$rect$h)
      legend( coords, items,
         col=col.items, pch=pch, bty="n", cex=1, xpd=NA)
  }
  if(location == "separate") {
    plot(0:1, 0:1, xaxt="n", yaxt="n", xlab="", ylab="", type="n")
    bg <- par()$bg
    legend("center", items,
          lty=lty.items, lwd=lwd.items, fill=col.items, border=col.items,,
          col=col.items, box.col=par()$bg,
          title=title)
  }
}

termsToMf <- function(model, terms){
  gform <- function(formula) {
    if (is.null(formula))
      return(list(vars=formula, groups=NULL))
    # is formula one-sided?
    if(length(formula) == 3) stop("terms must be a one-sided formula")
    rhs <- formula[[2]]
    # is '|' present in the formula?
    if("|" %in% all.names(rhs)){
      if(length(rhs[[3]]) > 1) stop("only one conditional variable permitted")
      groups <- as.formula(paste("~ ", deparse(rhs[[3]])))
      vars <- as.formula(paste("~", deparse(rhs[[2]])))} else{
        groups <- NULL
        vars <- formula
      }
    list(vars=vars, groups=groups)
  }
  terms <- gform(as.formula(terms))
  mf.vars <- try(update(model, terms$vars, method="model.frame"),
                 silent=TRUE)
  # This second test is used for models like m1 <- lm(longley) which
  # fail the first test because update doesn't work
  if(inherits(mf.vars, "try-error"))
    mf.vars <- try(update(model, terms$vars,
                          method="model.frame", data=model.frame(model)), silent=TRUE)
  if(inherits(mf.vars, "try-error")) stop("argument 'terms' not interpretable.")
  if(!is.null(terms$groups)){
    mf.groups <- try(update(model, terms$groups, method="model.frame"), silent=TRUE)
    if(inherits(mf.groups, "try-error"))
      mf.groups <- try(update(model, terms$groups,
                              method="model.frame", data=model.frame(model)), silent=TRUE)
    if(inherits(mf.groups, "try-error")) stop("argument 'terms' not interpretable.")
  } else {mf.groups <- NULL}
  list(mf.vars=mf.vars, mf.groups=mf.groups)
}

# the following function isn't exported, tests for existance of a package:

package.installed <- function(package){
  package <- as.character(substitute(package))
  result <- try(find.package(package), silent=TRUE)
  !inherits(result,  "try-error")
}

# support for coxme objects

model.matrix.coxme <- function(object, ...){
    if (!requireNamespace("survival")) stop("survival package is missing")
    class(object) <- "coxph"
    model.matrix(object)
}

alias.coxme <- function(object, ...){
    if(any(which <- is.na(coef(object)))) return(list(Complete=which))
    else list()
}

# to make linearHypothesis() work again and to make Anova() work with VGAM:"vglm" objects

# df.residual.vglm <- function(object, ...) object@df.residual

# vcov.vglm <- function(object, ...) vcovvlm(object, ...)

# coef.vglm <- function(object, ...) coefvlm(object, ...)

has_intercept.vlm <- function(model, ...) any(grepl("^\\(Intercept\\)", names(coef(model))))

# formula.vglm <- function(x, ...) formulavlm(x = x, ...)

# model.matrix.vglm <- function(object, ...) model.matrixvlm(object, ...)

# for plotting functions, not exported:

isFALSE <- function(x) length(x) == 1 && is.logical(x) && !isTRUE(x)

applyDefaults <- function(args, defaults, type=""){
    if (isFALSE(args)) return(FALSE)
    names <- names(args)
    names <- names[names != ""]
    if (!isTRUE(args) && !is.null(args) && length(names) != length(args)) warning("unnamed ", type, " arguments, will be ignored")
    if (isTRUE(args) || is.null(names)) defaults
    else defaults[names] <- args[names]
    as.list(defaults)
}

# carPal <- function(){
#     car.palette <- default <- c("black", "blue", "magenta", "cyan", "orange", "gray", "green3", "red")
#     function(palette){
#         if (missing(palette)) return(car.palette)
#         else{
#             previous <- car.palette
#             car.palette <<- if( palette[1] == "default") default else palette
#             return(invisible(previous))
#         }
#     }
# }

carPal <- function(){
  car.palette <- default <- c("black", "blue", "magenta", "cyan", "orange", "gray", "green3", "red")
  colorblind <- rgb(red = c(0, 230, 86, 0, 240, 0, 213, 204),
                    green = c(0, 159, 180, 158, 228, 114, 94, 121),
                    blue  = c(0, 0, 233, 115, 66, 178, 0, 167),
                    names = c("black", "orange", "sky.blue", "bluish.green", "yellow", 
                              "blue", "vermillion", "reddish.purple"),
                    maxColorValue = 255)
  # colorblind palette from https://jfly.uni-koeln.de/color/
  function(palette){
    if (missing(palette)) return(car.palette)
    else{
      previous <- car.palette
      car.palette <<- if (palette[1] %in% c("default", "car")) {
        default
      } else if (palette[1] == "colorblind") {
        colorblind
      } else if (palette[1] == "R"){
        palette()
      } else {
        palette
      }
      return(invisible(previous))
    }
  }
}

carPalette <- carPal()

# the following function borrowed from stats:::format.perc(), not exported
format_perc <- function (probs, digits){
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

# the following unexported function is useful for combining results of parallel computations

combineLists <- function(..., fmatrix="list", flist="c", fvector="rbind", 
                         fdf="rbind", recurse=FALSE){
    # combine lists of the same structure elementwise
    
    # ...: a list of lists, or several lists, each of the same structure
    # fmatrix: name of function to apply to matrix elements
    # flist: name of function to apply to list elements
    # fvector: name of function to apply to data frame elements
    # recurse: process list element recursively
    
    frecurse <- function(...){
        combineLists(..., fmatrix=fmatrix, fvector=fvector, fdf=fdf, 
                     recurse=TRUE)
    }
    
    if (recurse) flist="frecurse"
    list.of.lists <- list(...)
    if (length(list.of.lists) == 1){
        list.of.lists <- list.of.lists[[1]]
        list.of.lists[c("fmatrix", "flist", "fvector", "fdf")] <- 
            c(fmatrix, flist, fvector, fdf)
        return(do.call("combineLists", list.of.lists))
    }
    if (any(!sapply(list.of.lists, is.list))) 
        stop("arguments are not all lists")
    len <- sapply(list.of.lists, length)
    if (any(len[1] != len)) stop("lists are not all of the same length")
    nms <- lapply(list.of.lists, names)
    if (any(unlist(lapply(nms, "!=", nms[[1]])))) 
        stop("lists do not all have elements of the same names")
    nms <- nms[[1]]
    result <- vector(len[1], mode="list")
    names(result) <- nms
    for(element in nms){
        element.list <- lapply(list.of.lists, "[[", element)
#        clss <- sapply(element.list, class)
        clss <- lapply(element.list, class)
#        if (any(clss[1] != clss)) stop("list elements named '", element,
        if (!all(vapply(clss, function(e) all(e == clss[[1L]]), NA)))
          stop("list elements named '", element, "' are not all of the same class")
        
        is.df <- is.data.frame(element.list[[1]])
        fn <- if (is.matrix(element.list[[1]])) fmatrix 
        else if (is.list(element.list[[1]]) && !is.df) flist 
        else if (is.vector(element.list[[1]])) fvector
        else if (is.df) fdf
        else stop("list elements named '", element, 
                  "' are not matrices, lists, vectors, or data frames")
        result[[element]] <- do.call(fn, element.list)
    }
    result
}

matchFun <- function(name){
  object <- getFromNamespace(name, ns = "car")
  if (!is.function(object)) stop("'", name, "' is not a function")
  object
}

envelope <- function(x.low, x.up=x.low, lower, upper, col=1, lty=1, lwd=1, 
                     alpha=0.15, border=TRUE){
  color <- as.vector(col2rgb(col))/255
  polygon(c(x.up, rev(x.low)), c(upper, rev(lower)), 
          col=rgb(red=color[1], green=color[2], blue=color[3], alpha=alpha),
          border=if (border) rgb(red=color[1], green=color[2], blue=color[3]) else NA,
          lty=lty, lwd=lwd)
}

getVcov <- function(v, mod, ...){
  if(missing(v)) return(vcov(mod, ...)) 
  if(inherits(v, "matrix")) return(v)
  if(is.function(v)) return(v(mod, ...)) 
  if(is.null(v)) return(vcov(mod, ...))
  v <- try(as.matrix(v), silent=TRUE)
  if (is.matrix(v)) return(v)
  stop("vcov. must be a matrix or a function")
}

getModelData <- function(model) {
  # returns a data frame with the data to which the model was fit
  # model: a statistical model object that responds to model.frame() and formula() 
  data1 <- data <- model.frame(model)
  vars <- all.vars(formula(model))
  if ("pi" %in% vars) {
    vars <- setdiff(vars, "pi")
    message("the symbol 'pi' is treated as a numeric constant in the model formula")
  }
  cols <- colnames(data)
  check <- vars %in% cols
  if (!(all(check))) {
    missing.cols <- !check
    data1 <- expand.model.frame(model, vars[missing.cols])
  }
  missing.cols <- !cols %in% colnames(data1)
  if (any(missing.cols)) {
    data1 <- cbind(data1, data[missing.cols])
  }
  cols <- colnames(data1)
  valid <- make.names(cols) == cols | grepl("^\\(.*\\)$", cols)
  data1[valid]
}
