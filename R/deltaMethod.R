#-------------------------------------------------------------------------------
# Revision history:
# 2009-10-29: renamed var argument to .vcov; tidied code. John
# 2010-07-02; added method for survreg and coxph objects.
# 2010-07-02; rewrote default method to permit prarmeter names to have
#   meta-characters 
#-------------------------------------------------------------------------------

deltaMethod <- function (object, ...) {
	UseMethod("deltaMethod")
}

 deltaMethod.default <- function (object, g, vcov., func=g, ...) {
	if (!is.character(g)) 
		stop("The argument 'g' must be a character string")
	metas <- c("(", ")", "[", "]", "{", "}", ".", "*", "+", "^", "$", ":", "|")
	metas2 <- paste("\\", metas, sep="")
	metas3 <- paste("\\\\", metas, sep="")
	func <- func
	para <- object
	para.names <- names(para)
	for (i in seq(along=metas))
		para.names <- gsub(metas2[i], metas3[i], para.names) # fix up metacharacters
	para.order <- order(nchar(para.names), decreasing=TRUE) 
	para.names <- para.names[para.order] # avoid partial-name substitution
	std.names <- paste("Param", 1:length(para), sep = "")
	std.names.ordered <- std.names[para.order]
	for (i in seq(along=para.names)){
		g <- gsub(para.names[i], std.names.ordered[i], g) 
	}
	names(para) <- std.names
	g <- parse(text = g)
	q <- length(para)
	for (i in 1:q) {
		assign(names(para)[i], para[i])
	}   
	est <- eval(g)
	names(est) <- NULL
	gd <- NULL
	for (i in 1:q) {
		gd <- c(gd, eval(D(g, names(para)[i])))
	}
	se.est <- as.vector(sqrt(t(gd) %*% vcov. %*% gd))
	data.frame(Estimate = est, SE = se.est, row.names = c(func))
}

deltaMethod.lm <- function (object, g, vcov. = vcov, parameterPrefix = "b", ...) {
	metas <- c("(", ")", "[", "]", "{", "}", ".", "*", "+", "^", "$", ":", "|")
	metas2 <- paste("\\", metas, sep="")
	metas3 <- paste("\\\\", metas, sep="")
	para <- coef(object)
	para.names <- names(para)
	for (i in seq(along=metas))
		para.names <- gsub(metas2[i], metas3[i], para.names) # fix up metacharacters
	para.order <- order(nchar(para.names), decreasing=TRUE) 
	para.names <- para.names[para.order] # avoid partial-name substitution
	std.names <- if ("(Intercept)" %in% names(para)) 
			paste(parameterPrefix, 0:(length(para) - 1), sep = "")
		else paste(parameterPrefix, 1:length(para), sep = "")
	std.names.ordered <- std.names[para.order]
	func <- g
	for (i in seq(along=para.names)){
		g <- gsub(para.names[i], std.names.ordered[i], g) 
	}
	vcov. <- if (is.function(vcov.)) 
			vcov.(object)
		else vcov.
	names(para) <- std.names
	deltaMethod.default(para, g, vcov., func)
}

# nls has named parameters so parameterPrefix is ignored
deltaMethod.nls <- function(object, g, vcov.=vcov,...){
	vcov. <- if(is.function(vcov.)) vcov.(object)
	deltaMethod.default(coef(object), g, vcov.)
}

deltaMethod.multinom <- function(object, g, vcov.=vcov, parameterPrefix="b", ...){
	out <- NULL
	coefs <- coef(object)
	if (!is.matrix(coefs)) {
		nn <- names(coefs)
		coefs <- matrix(coefs, nrow=1)
		colnames(coefs) <- nn
	}
	nc <- dim(coefs)[2]
	for (i in 1:dim(coefs)[1]){
		para <- coefs[i,]
		names(para) <- if ("(Intercept)" %in% names(para))
				paste(parameterPrefix, 0:(length(para)-1), sep="") else
				paste(parameterPrefix, 1:length(para), sep="")
		ans <- deltaMethod.default(para, g, vcov.(object)[(i - 1) + 1:nc, (i - 1) + 1:nc])
		rownames(ans)[1] <- paste(rownames(coefs)[i], rownames(ans)[1])
		out <- rbind(out,ans)
	}
	out}

deltaMethod.polr <- function(object,g,vcov.=vcov,...){
	sel <- 1:(length(coef(object)))
	vcov. <- if(is.function(vcov.)) vcov.(object)[sel, sel]
	deltaMethod.lm(object, g, vcov., ...)
}

# method for survreg objects.
deltaMethod.survreg <- function (object, g, vcov. = vcov,  ...) {
  para <- c(coef(object), object$icoef[2])
 	vcov. <- if (is.function(vcov.)) 
			vcov.(object)
		else vcov. 
  deltaMethod.default(para, g, vcov.)
  }

 # method for coxph objects.
deltaMethod.coxph<- function(object, g, vcov.=vcov,...){
      sel <- 1:(length(coef(object)))
      vcov. <- if(is.function(vcov.)) vcov.(object)[sel, sel]
      deltaMethod(coef(object), g, vcov., ...)
}