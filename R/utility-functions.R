
# Utility functions (J. Fox)

# 16 March 2010 changed 'vars' argument to 'terms'
# 28 June 2010 added df.terms.surveg and model.matrix.survreg
# 15 November 2010 added squeezeBlanks
# 21 January 2011 added functions to support mixed models
# 2012-04-08 added exists.method
# 20120-06-23: added call to globalVariables(). John

if (getRversion() >= "2.15.1") globalVariables(c(".boot.sample", ".boot.indices"))

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

has.intercept <- function (model, ...) {
	UseMethod("has.intercept")
}

has.intercept.default <- function(model, ...) any(names(coefficients(model))=="(Intercept)")

term.names <- function (model, ...) {
	UseMethod("term.names")
}

term.names.default <- function (model, ...) {
	term.names <- labels(terms(model))
	if (has.intercept(model)) c("(Intercept)", term.names)
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
	nlev <- length(model$lev)
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
	id.method = c("mahal", "identify", "none"), log="", id.cex=.75, id.n=3, id.col=palette()[1], 
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

df.residual.lme <- function(object, ...) Inf

has.intercept.mer <- function(model){
	any(names(fixef(model))=="(Intercept)")
}
	
model.matrix.lme <- function(object, ...){
	model.matrix(as.formula(object$call$fixed), eval(object$call$data))
}

# added by J. Fox 2012-04-08 to use in deltaMethod.default()

exists.method <- function(generic, object, default=TRUE, strict=FALSE){
	classes <- class(object)
	if (default) classes <- c(classes, "default")
	if (strict) classes <- classes[1]
	any(paste(generic, ".", classes, sep="") %in%
					as.character(methods(generic)))
}
