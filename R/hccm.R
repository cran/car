#-------------------------------------------------------------------------------
# Revision history:
# 2009-09-16: optionally allow models with aliased coefficients J. Fox
# 2012-04-04: modified to allow weighted linear models. J. Fox
# 2020-06-25: Fix bug in hccm.lm() when model matrix includes just one column 
#             (reported by Justin Yap). J. Fox
# 2021-07-29: Report error when any hatvalue = 1 for all but hc0 and hc1
#             (following report of problem reported by Peng Ding) J. Fox
#-------------------------------------------------------------------------------

# Heteroscedasticity-corrected standard errors (Huber/White adjustment) (J. Fox)

hccm <- function(model, ...){
	UseMethod("hccm")
}

hccm.lm <-function (model, type = c("hc3", "hc0", "hc1", "hc2", "hc4"), 
		singular.ok = TRUE, ...) {
	e <- na.omit(residuals(model))
	removed <- attr(e, "na.action")
	wts <- if (is.null(weights(model))) 1 
			else weights(model)
	type <- match.arg(type)
	if (any(aliased <- is.na(coef(model))) && !singular.ok) 
		stop("there are aliased coefficients in the model")
	sumry <- summary(model, corr = FALSE)
	s2 <- sumry$sigma^2
	V <- sumry$cov.unscaled
	if (type == FALSE) 
		return(s2 * V)
	h <- hatvalues(model)
	if (!is.null(removed)){
		wts <- wts[-removed]
		h <- h[-removed]
	}
	X <- model.matrix(model)[, !aliased, drop=FALSE]
	df.res <- df.residual(model)
	n <- length(e)
	e <- wts*e
	p <- ncol(X)
	factor <- switch(type, hc0 = 1, hc1 = df.res/n,
			hc2 = 1 - h, hc3 = (1 - h)^2, hc4 = (1 - h)^pmin(4, n * h/p))
	V <- V %*% t(X) %*% apply(X, 2, "*", (e^2)/factor) %*% V
	bad <- h > 1 - sqrt(.Machine$double.eps)
	if ((n.bad <- sum(bad)) > 0 && !(type %in% c("hc0", "hc1"))) {
	  nms <- names(e)
	  bads <- if (n.bad <= 10) {
	    paste(nms[bad], collapse=", ")
	  } else {
	    paste0(paste(nms[bad[1:10]], collapse=", "), ", ...")
	  }
	  if (any(is.nan(V))){
	    stop("cannot proceed because of ", n.bad, if (n.bad == 1) " case " else " cases ",
	       "with hatvalue = 1:\n   ", bads)
	  } else {
	    warning("adjusted coefficient covariances may be unstable because of ", n.bad, 
	            if (n.bad == 1) " case " else " cases ",
	            "with hatvalue near 1:\n   ", bads)
	  }
	}
	V
}

hccm.default<-function(model, ...){
	stop("requires an lm object")
}
