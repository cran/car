#-------------------------------------------------------------------------------
# Revision history:
# 2009-10-29: renamed var argument to .vcov; tidied code. John
# 2010-07-02; added method for survreg and coxph objects.
# 2010-07-02; rewrote default method to permit parameter names to have
#   meta-characters 
# 2011-07028  Removed meta-character checks; removed parameterPrefix because
#   it didn't work and caused problems; added parameterNames to restore the
#   utility of parameterPrefix
# 2011-10-02 Fixed bugs in the .survreg and .coxph methods so parameterNames
#   works correctly
# 2012-03-02: fixed abbreviation of envir argument. J. Fox
# 2012-04-08: modfied deltaMethod.default() to use coef and vcov
# 2012-12-10: removed the 'deltaMethodMessageFlag'
# 2013-06-20: added deltaMethod.merMod(). J. Fox
# 2013-06-20: tweaks for lme4. J. Fox
# 2013-07-01: New 'constants' argument for use when called from within a function.
# 2013-07-18: fixed a bug in passing the 'func' argument
# 2016-03-31: added level argument and report CIs. J. Fox
# 2017-11-09: made compatible with vcov() in R 2.5.0. J. Fox
# 2017-11-29: further fixes for vcov() and vcov.(). J. Fox
# 2017-12-01: fix bug in handling vcov. arg in some methods. J. Fox
# 2019-01-16: changed g arg to g. to allow variable named "g". J. Fox
# 2019-06-03: introduction of environment to hold coefficients and constants. Pavel Krivitsky
# 2019-06-05: option for hypothesis test. J. Fox
# 2019-06-07: move handling intercepts to default method, suggestion of Pavel Krivitsky. J. Fox 
#-------------------------------------------------------------------------------

deltaMethod <- function (object, ...) {
	UseMethod("deltaMethod")
}

deltaMethod.default <- function (object, g., vcov., func = g., constants, level=0.95, rhs=NULL, ..., envir=parent.frame()) {
  if (!is.character(g.)) 
    stop("The argument 'g.' must be a character string")
  if ((exists.method("coef", object, default=FALSE) ||
       (!is.atomic(object) && !is.null(object$coefficients))) 
      && exists.method("vcov", object, default=FALSE)){
    if (missing(vcov.)) vcov. <- vcov(object, complete=FALSE)
    object <- coef(object)
  }
	para <- object         
	para.names <- names(para)
	para.names[1] <- gsub("\\(Intercept\\)", "Intercept", para.names[1])
	g. <- parse(text = g.)
	q <- length(para)

	envir <- new.env(parent=envir)
	for (i in 1:q) {
	    assign(para.names[i], para[i], envir)
	}
	if(!missing(constants)){
     for (i in seq_along(constants)) assign(names(constants[i]), constants[[i]], envir)}
	est <- eval(g., envir)
	names(est) <- NULL
	gd <- rep(0, q)
	for (i in 1:q) {
		gd[i] <- eval(D(g., names(para)[i]), envir)
	}
	se.est <- as.vector(sqrt(t(gd) %*% vcov. %*% gd))
	result <- data.frame(Estimate = est, SE = se.est, row.names = c(func))
	p <- (1 - level)/2
	z <- - qnorm(p)
	lower <- est - z*se.est
	upper <- est + z*se.est
	pct <- paste(format(100*c(p, 1 - p), trim=TRUE, scientific=FALSE, digits=3), "%")
	result <- cbind(result, lower, upper)
	names(result)[3:4] <- pct
	if (!is.null(rhs)){
	    z <- (est - rhs)/se.est
	    p <- 2*(pnorm(abs(z), lower.tail=FALSE))
	    result <- cbind(result, "Hypothesis"=rhs, "z value"=z, "Pr(>|z|)"=p)
	}
	class(result) <- c("deltaMethod", class(result))
	result
}

print.deltaMethod <- function(x, ...){
    if (ncol(x) == 3) print(x, ...) else printCoefmat(x, ...)
    invisible(x)
}

deltaMethod.lm <- function (object, g., vcov. = vcov(object, complete=FALSE), 
           parameterNames = names(coef(object)), ..., envir=parent.frame()) {
	para <- coef(object)
	para.names <- parameterNames
	names(para) <- para.names
	vcov. <- if (is.function(vcov.)) 
			vcov.(object)
		else vcov.
	deltaMethod.default(para, g., vcov.,  ..., envir=envir)
}


# nls has named parameters so parameterNames is ignored
deltaMethod.nls <- function(object, g., vcov.=vcov(object, complete=FALSE), ..., envir=parent.frame()){
	vcov. <- if(is.function(vcov.)) vcov.(object) else vcov.
	deltaMethod.default(coef(object), g., vcov., ..., envir=envir)   
}

deltaMethod.polr <- function(object,g.,vcov.=vcov(object, complete=FALSE), ..., envir=parent.frame()){
	sel <- 1:(length(coef(object)))
	vcov. <- if(is.function(vcov.)) vcov.(object)[sel, sel] else vcov.[sel, sel]
	deltaMethod.lm(object, g., vcov., ..., envir=envir)
}

deltaMethod.multinom <- function(object, g., vcov.=vcov(object, complete=FALSE), 
   parameterNames = if(is.matrix(coef(object)))
     colnames(coef(object)) else names(coef(object)), ..., envir=parent.frame()){
  vcov. <- if(is.function(vcov.)) vcov.(object) else vcov.
	out <- NULL
	coefs <- coef(object)
	if (!is.matrix(coefs)) { coefs <- t(as.matrix(coefs)) }
	colnames(coefs) <- parameterNames
	nc <- dim(coefs)[2]
	for (i in 1:dim(coefs)[1]){
		para <- coefs[i, ]
		ans <- deltaMethod(para, g., vcov.[(i - 1) + 1:nc, (i - 1) + 1:nc], ..., envir=envir)
		rownames(ans)[1] <- paste(rownames(coefs)[i], rownames(ans)[1])
		out <- rbind(out,ans)
	}
	out}

# method for survreg objects. 
deltaMethod.survreg <- function(object, g., vcov. = vcov(object, complete=FALSE), 
           parameterNames = names(coef(object)), ..., envir=parent.frame()) {
 deltaMethod.lm(object, g., vcov., parameterNames , ..., envir=envir) }


 # method for coxph objects.
deltaMethod.coxph <- function(object, g., vcov. = vcov(object, complete=FALSE), 
           parameterNames = names(coef(object)), ..., envir=parent.frame()) {
 deltaMethod.lm(object, g., vcov.,  parameterNames, ..., envir=envir) }
           
# lmer

deltaMethod.merMod <- function(object, g., vcov. = vcov(object, complete=FALSE),
                            parameterNames = names(fixef(object)), ..., envir=parent.frame()) {
    deltaMethod.mer(object=object, g.=g., vcov.=vcov, 
                       parameterNames=parameterNames, ..., envir=envir)
}
    
deltaMethod.mer <- function(object, g., vcov. = vcov(object, complete=FALSE),
           parameterNames = names(fixef(object)), ..., envir=parent.frame()) {
  para <- fixef(object)
  names(para) = parameterNames
 	vcov. <- if (is.function(vcov.)) 
			vcov.(object)
		else vcov.
  deltaMethod(para, g., vcov., ..., envir=envir)
  }


#lme
deltaMethod.lme <- function(object, g., vcov. = vcov(object, complete=FALSE),
           parameterNames = names(fixef(object)), ..., envir=parent.frame()) {
  para <- fixef(object)
  names(para) = parameterNames
 	vcov. <- if (is.function(vcov.)) 
			vcov.(object)
		else vcov.
  deltaMethod(para, g., vcov., ..., envir=envir)
  }
  
# nlsList  lsList
deltaMethod.lmList <- function(object, g., ..., envir=parent.frame()) {
  out <- t(sapply(object, function(x) deltaMethod(x, g., ..., envir=envir)))
  rownames(out) <- paste(rownames(out), g.)
  out
  }
  

