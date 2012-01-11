#-------------------------------------------------------------------------------
# Revision history:
# 2009-01-05: bug fix in Anova.II.lm(). J. Fox
# 2009-01-16: Cox models with clusters now handled. J. Fox
# 2009-09-16: reworked glm and lm methods to handle aliased parameters. J. Fox
# 2009-09-30: renamed "Anova" to "Analysis of Deviance" in output for some methods. J. Fox
# 2009-12-22: modified Anova.mlm() to handle a user-supplied within-subject model matrix. J. Fox
# 2009-12-28: named the components of P in Anova.III.mlm(). John
# 2010-01-01: Anova.II.mlm() now hands off (again) to Anova.III.mlm() when there
#             is only an intercept in the between-subjects model
# 2010-02-17: Fixed bug that caused some models with aliased coefficients to fail. J. Fox
# 2010-06-14: added wcrossprod and allow use of observation weights in Anova.mlm()
# 2010-06-28: Fixed Anova() tables for coxph and survreg models 
#             (failed because of changes in survival package.
# 2011-01-21: Added functions for mixed models. J. Fox
# 2011-01-25: Fixed Anova.polr() and Anova.multinom() to work with models with only one term. J. Fox
# 2011-05-19: local fixef() to avoid nlme/lme4 issues. J. Fox
# 2011-05-11: changed order of columns in ANOVA tables for mixed models. J. Fox
# 2011-11-27: added Anova.svyglm(). J. Fox
# 2011-12-31: fixed bug in Anova.II(and III).F.glm() when na.exclude used. J. Fox
#-------------------------------------------------------------------------------

# Type II and III tests for linear, generalized linear, and other models (J. Fox)

ConjComp <- function(X, Z = diag( nrow(X)), ip = diag(nrow(X))) {
	# This function by Georges Monette
	# finds the conjugate complement of the proj of X in span(Z) wrt
	#    inner product ip
	# - assumes Z is of full column rank
	# - projects X conjugately wrt ip into span Z
	xq <- qr(t(Z) %*% ip %*% X)
	if (xq$rank == 0) return(Z)
	Z %*% qr.Q(xq, complete = TRUE) [ ,-(1:xq$rank)] 
}

relatives <- function(term, names, factors){
	is.relative <- function(term1, term2) {
		all(!(factors[,term1]&(!factors[,term2])))
	}
	if(length(names) == 1) return(NULL)
	which.term <- which(term==names)
	(1:length(names))[-which.term][sapply(names[-which.term], 
					function(term2) is.relative(term, term2))]
}


Anova <- function(mod, ...){
	UseMethod("Anova", mod)
}

# linear models

Anova.lm <- function(mod, error, type=c("II","III", 2, 3), 
		white.adjust=c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"), 
		singular.ok, ...){
	type <- as.character(type)
	white.adjust <- as.character(white.adjust)
	type <- match.arg(type)
	white.adjust <- match.arg(white.adjust)
	if (missing(singular.ok)){
		singular.ok <- type == "2" || type == "II"
	}
	if (has.intercept(mod) && length(coef(mod)) == 1 
			&& (type == "2" || type == "II")) {
		type <- "III"
		warning("the model contains only an intercept: Type III test substituted")
	}
	if (white.adjust != "FALSE"){
		if (white.adjust == "TRUE") white.adjust <- "hc3" 
		return(Anova.default(mod, type=type, vcov.=hccm(mod, type=white.adjust), test="F", 
						singular.ok=singular.ok))
	}
	switch(type,
			II=Anova.II.lm(mod, error, singular.ok=singular.ok, ...),
			III=Anova.III.lm(mod, error, singular.ok=singular.ok, ...),
			"2"=Anova.II.lm(mod, error, singular.ok=singular.ok, ...),
			"3"=Anova.III.lm(mod, error, singular.ok=singular.ok,...))
}

Anova.aov <- function(mod, ...){
	class(mod) <- "lm"
	Anova.lm(mod, ...)
}

Anova.II.lm <- function(mod, error, singular.ok=TRUE, ...){
	if (!missing(error)){
		sumry <- summary(error, corr=FALSE)
		s2 <- sumry$sigma^2
		error.df <- error$df.residual
		error.SS <- s2*error.df
	}
	SS.term <- function(term){
		which.term <- which(term == names)
		subs.term <- which(assign == which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives) 
			subs.relatives <- c(subs.relatives, which(assign == relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov(mod)))
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(SS=NA, df=0))
		lh <- linearHypothesis(mod, hyp.matrix.term, 
				singular.ok=singular.ok, ...)
		abs(c(SS=lh$"Sum of Sq"[2], df=lh$Df[2]))
	}
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(mod$terms, "factors")
	intercept <- has.intercept(mod)
	I.p <- diag(length(coefficients(mod)))
	assign <- mod$assign
	assign[!not.aliased] <- NA
	names <- term.names(mod)
	if (intercept) names <-names[-1]
	n.terms <- length(names)
	p <- df <- f <- SS <- rep(0, n.terms + 1)
	sumry <- summary(mod, corr = FALSE)
	SS[n.terms + 1] <- if (missing(error)) sumry$sigma^2*mod$df.residual 
			else error.SS   
	df[n.terms + 1] <- if (missing(error)) mod$df.residual else error.df
	p[n.terms + 1] <- f[n.terms + 1] <- NA
	for (i in 1:n.terms){
		ss <- SS.term(names[i])
		SS[i] <- ss["SS"]
		df[i] <- ss["df"]
		f[i] <- df[n.terms+1]*SS[i]/(df[i]*SS[n.terms + 1])
		p[i] <- pf(f[i], df[i], df[n.terms + 1], lower.tail = FALSE)
	}    
	result <- data.frame(SS, df, f, p)
	row.names(result) <- c(names,"Residuals")
	names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Anova Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

# type III

Anova.III.lm <- function(mod, error, singular.ok=FALSE, ...){
	if (!missing(error)){
		error.df <- df.residual(error)
		error.SS <- deviance(error)
	}
	else {
		error.df <- df.residual(mod)
		error.SS <- deviance(mod)
	}
	intercept <- has.intercept(mod)
	I.p <- diag(length(coefficients(mod)))
	Source <- term.names(mod)
	n.terms <- length(Source)
	p <- df <- f <- SS <- rep(0, n.terms + 1)
	assign <- mod$assign
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	for (term in 1:n.terms){
		subs <- which(assign == term - intercept)
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix) == 0){
			SS[term] <- NA
			df[term] <- 0
			f[term] <- NA
			p[term] <- NA
		}
		else {
			test <- if (missing(error)) linearHypothesis(mod, hyp.matrix, 
								singular.ok=singular.ok, ...)
					else linearHypothesis(mod, hyp.matrix, error.SS=error.SS, error.df=error.df, 
								singular.ok=singular.ok, ...)
			SS[term] <- test$"Sum of Sq"[2]
			df[term] <- test$"Df"[2]
			f[term] <- test$"F"[2]
			p[term] <- test$"Pr(>F)"[2]
		}
	}
	Source[n.terms + 1] <- "Residuals"
	SS[n.terms + 1] <- error.SS
	df[n.terms + 1] <- error.df
	p[n.terms + 1] <- f[n.terms + 1] <- NA
	result <- data.frame(SS, df, f, p)
	row.names(result) <- Source
	names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
	result
}

# generalized linear models

Anova.glm <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald", "F"), 
		error, error.estimate=c("pearson", "dispersion", "deviance"), singular.ok, ...){
	type <- as.character(type)
	type <- match.arg(type)
	if (has.intercept(mod) && length(coef(mod)) == 1 
			&& (type == "2" || type == "II")) {
		type <- "III"
		warning("the model contains only an intercept: Type III test substituted")
	}
	if (missing(singular.ok)){
		singular.ok <- type == "2" || type == "II"
	}
	test.statistic <- match.arg(test.statistic)
	error.estimate <- match.arg(error.estimate)
	switch(type,
			II=switch(test.statistic,
					LR=Anova.II.LR.glm(mod, singular.ok=singular.ok),
					Wald=Anova.default(mod, type="II", singular.ok=singular.ok),
					F=Anova.II.F.glm(mod, error, error.estimate, singular.ok=singular.ok)),
			III=switch(test.statistic,
					LR=Anova.III.LR.glm(mod, singular.ok=singular.ok),
					Wald=Anova.default(mod, type="III", singular.ok=singular.ok),
					F=Anova.III.F.glm(mod, error, error.estimate, singular.ok=singular.ok)),
			"2"=switch(test.statistic,
					LR=Anova.II.LR.glm(mod, singular.ok=singular.ok),
					Wald=Anova.default(mod, type="II", singular.ok=singular.ok),
					F=Anova.II.F.glm(mod, error, error.estimate, singular.ok=singular.ok)),
			"3"=switch(test.statistic,
					LR=Anova.III.LR.glm(mod, singular.ok=singular.ok),
					Wald=Anova.default(mod, type="III", singular.ok=singular.ok),
					F=Anova.III.F.glm(mod, error, error.estimate, singular.ok=singular.ok)))
}


# type III

# LR test

Anova.III.LR.glm <- function(mod, singular.ok=FALSE, ...){
	if (!singular.ok && any(is.na(coef(mod))))
		stop("there are aliased coefficients in the model")
	Source <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(Source)
	p <- df <- LR <- rep(0, n.terms)
	dispersion <- summary(mod, corr = FALSE)$dispersion
	deviance <- deviance(mod)/dispersion
	for (term in 1:n.terms){
		mod.1 <- drop1(mod, scope=eval(parse(text=paste("~",Source[term]))))
		LR[term] <- (mod.1$Deviance[2]/dispersion)-deviance
		df[term] <- mod.1$Df[2]
		p[term] <- pchisq(LR[term], df[term], lower.tail = FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- Source
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova","data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", paste("Response:", responseName(mod)))
	result
}

# F test

Anova.III.F.glm <- function(mod, error, error.estimate, singular.ok=FALSE, ...){
	if (!singular.ok && any(is.na(coef(mod))))
		stop("there are aliased coefficients in the model")
	fam <- family(mod)$family
	if (fam == "binomial" || fam == "poisson") 
		warning("dispersion parameter estimated from the Pearson residuals, not taken as 1")
	if (missing(error)) error <- mod
	df.res <- df.residual(error)
	error.SS <- switch(error.estimate,
			pearson=sum(residuals(error, "pearson")^2, na.rm=TRUE),
			dispersion=df.res*summary(error, corr = FALSE)$dispersion,
			deviance=deviance(error))
	Source <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(Source)
	p <- df <- f <- SS <-rep(0, n.terms+1)
	f[n.terms+1] <- p[n.terms+1] <- NA
	df[n.terms+1] <- df.res
	SS[n.terms+1] <- error.SS
	dispersion <- error.SS/df.res
	deviance <- deviance(mod)
	for (term in 1:n.terms){
		mod.1 <- drop1(mod, scope=eval(parse(text=paste("~",Source[term]))))
		df[term] <- mod.1$Df[2]
		SS[term] <- mod.1$Deviance[2] - deviance
		f[term] <- (SS[term]/df[term])/dispersion
		p[term] <- pf(f[term], df[term], df.res, lower.tail = FALSE)
	}
	result <- data.frame(SS, df, f, p)
	row.names(result) <- c(Source, "Residuals")
	names(result) <- c("SS", "Df", "F", "Pr(>F)")
	class(result) <- c("anova","data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", paste("Response:", responseName(mod)))
	result
}

# type II

# LR test

Anova.II.LR.glm <- function(mod, singular.ok=TRUE, ...){
	if (!singular.ok && any(is.na(coef(mod))))
		stop("there are aliased coefficients in the model")
	# (some code adapted from drop1.glm)
	which.nms <- function(name) which(asgn == which(names == name))
	fac <- attr(mod$terms, "factors")
	names <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- mod$y
	if (is.null(y)) y <- model.response(model.frame(mod), "numeric")
	wt <- mod$prior.weights
	if (is.null(wt)) wt <- rep(1, length(y))
	asgn <- attr(X, 'assign')
	df <- p <- LR <- rep(0, n.terms)
	dispersion <- summary(mod, corr = FALSE)$dispersion
	for (term in 1:n.terms){
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
		mod.1 <- glm.fit(X[, -exclude.1, drop = FALSE], y, wt, offset = mod$offset, 
				family = mod$family, control = mod$control)
		dev.1 <- deviance(mod.1)
		mod.2 <- if (length(rels) == 0) mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					glm.fit(X[, -exclude.2, drop = FALSE], y, wt, offset = mod$offset, 
							family = mod$family, control = mod$control)
				}
		dev.2 <- deviance(mod.2)
		df[term] <- df.residual(mod.1) - df.residual(mod.2)
		if (df[term] == 0) LR[term] <- p[term] <- NA
		else {
			LR[term] <- (dev.1 - dev.2)/dispersion
			p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
		}
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- 
			c("Analysis of Deviance Table (Type II tests)\n", paste("Response:", responseName(mod)))
	result
}


# F test

Anova.II.F.glm <- function(mod, error, error.estimate, singular.ok=TRUE, ...){
	# (some code adapted from drop1.glm)
	if (!singular.ok && any(is.na(coef(mod))))
		stop("there are aliased coefficients in the model")
	fam <- family(mod)$family
	if (fam == "binomial" || fam == "poisson") 
		warning("dispersion parameter estimated from the Pearson residuals, not taken as 1")
	which.nms <- function(name) which(asgn == which(names == name))
	if (missing(error)) error <- mod
	df.res <- df.residual(error)
	error.SS <- switch(error.estimate,
			pearson = sum(residuals(error, "pearson")^2, na.rm=TRUE),
			dispersion = df.res*summary(error, corr = FALSE)$dispersion,
			deviance = deviance(error))
	fac <- attr(mod$terms, "factors")
	names <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- mod$y
	if (is.null(y)) y <- model.response(model.frame(mod), "numeric")
	wt <- mod$prior.weights
	if (is.null(wt)) wt <- rep(1, length(y))
	asgn <- attr(X, 'assign')
	p <- df <- f <- SS <- rep(0, n.terms+1)
	f[n.terms+1] <- p[n.terms+1] <- NA
	df[n.terms+1] <- df.res
	SS[n.terms+1] <- error.SS
	dispersion <- error.SS/df.res
	for (term in 1:n.terms){
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
		mod.1 <- glm.fit(X[, -exclude.1, drop = FALSE], y, wt, offset = mod$offset, 
				family = mod$family, control = mod$control)
		dev.1 <- deviance(mod.1)
		mod.2 <- if (length(rels) == 0) mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					glm.fit(X[, -exclude.2, drop = FALSE], y, wt, offset = mod$offset, 
							family = mod$family, control = mod$control)
				}
		dev.2 <- deviance(mod.2)
		df[term] <- df.residual(mod.1) - df.residual(mod.2)
		if (df[term] == 0) SS[term] <- f[term] <- p[term] <- NA
		else {
			SS[term] <- dev.1 - dev.2
			f[term] <- SS[term]/(dispersion*df[term])
			p[term] <- pf(f[term], df[term], df.res, lower.tail=FALSE)
		}
	}
	result <- data.frame(SS, df, f, p)
	row.names(result) <- c(names, "Residuals")
	names(result) <- c("SS", "Df", "F", "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

# multinomial logit models (via multinom in the nnet package)

Anova.multinom <-
		function (mod, type = c("II", "III", 2, 3), ...)
{
	type <- as.character(type)
	type <- match.arg(type)
	if (has.intercept(mod) && length(coef(mod)) == 1 
			&& (type == "2" || type == "II")) {
		type <- "III"
		warning("the model contains only an intercept: Type III test substituted")
	}
	switch(type,
			II = Anova.II.multinom(mod, ...),
			III = Anova.III.multinom(mod, ...),
			"2" = Anova.II.multinom(mod, ...),
			"3" = Anova.III.multinom(mod, ...))
}

Anova.II.multinom <- function (mod, ...)
{
	which.nms <- function(name) which(asgn == which(names ==
								name))
	fac <- attr(mod$terms, "factors")
	names <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- model.response(model.frame(mod))
	wt <- mod$weights
	asgn <- attr(X, "assign")
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	for (term in 1:n.terms) {
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
								which.nms)))
		mod.1 <-if (n.terms > 1) multinom(y ~ X[, -c(1, exclude.1)], weights=wt, trace=FALSE)
				else multinom(y ~ 1, weights=wt, race=FALSE)
		dev.1 <- deviance(mod.1)
		mod.2 <- if (length(rels) == 0)
					mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					multinom(y ~ X[, -c(1, exclude.2)], weights=wt, trace=FALSE)
				}
		dev.2 <- deviance(mod.2)
		LR[term] <- dev.1 - dev.2
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n",
			paste("Response:", responseName(mod)))
	result
}

Anova.III.multinom <- function (mod, ...)
{
	names <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- model.response(model.frame(mod))
	wt <- mod$weights
	asgn <- attr(X, "assign")
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	deviance <- deviance(mod)
	for (term in 1:n.terms) {
		mod.1 <- if (n.terms > 1) multinom(y ~ X[, term != asgn][, -1], weights=wt, trace=FALSE)
				else multinom(y ~ 1, weights=wt, trace=FALSE)
		LR[term] <- deviance(mod.1) - deviance
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n",
			paste("Response:", responseName(mod)))
	result
}


# proportional-odds logit models (via polr in the MASS package)

Anova.polr <- function (mod, type = c("II", "III", 2, 3), ...)
{
	type <- as.character(type)
	type <- match.arg(type)
	if (has.intercept(mod) && length(coef(mod)) == 1 
			&& (type == "2" || type == "II")) {
		type <- "III"
		warning("the model contains only an intercept: Type III test substituted")
	}
	switch(type,
			II = Anova.II.polr(mod, ...),
			III = Anova.III.polr(mod, ...),
			"2" = Anova.II.polr(mod, ...),
			"3" = Anova.III.polr(mod, ...))
}

Anova.II.polr <- function (mod, ...)
{
	which.nms <- function(name) which(asgn == which(names ==
								name))
	fac <- attr(mod$terms, "factors")
	names <- term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- model.response(model.frame(mod))
	wt <- model.weights(model.frame(mod))
	asgn <- attr(X, "assign")
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	for (term in 1:n.terms) {
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
								which.nms)))
		mod.1 <- if (n.terms > 1) polr(y ~ X[, -c(1, exclude.1)], weights=wt)
				else polr(y ~ 1, weights=wt)
		dev.1 <- deviance(mod.1)
		mod.2 <- if (length(rels) == 0)
					mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					polr(y ~ X[, -c(1, exclude.2)], weights=wt)
				}
		dev.2 <- deviance(mod.2)
		LR[term] <- dev.1 - dev.2
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n",
			paste("Response:", responseName(mod)))
	result
}

Anova.III.polr <- function (mod, ...)
{
	names <- term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- model.response(model.frame(mod))
	wt <- model.weights(model.frame(mod))
	asgn <- attr(X, "assign")
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	deviance <- deviance(mod)
	for (term in 1:n.terms) {
		mod.1 <- if (n.terms > 1) polr(y ~ X[, term != asgn][, -1], weights=wt)
				else polr(y ~ 1, weights=wt)
		LR[term] <- deviance(mod.1) - deviance
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n",
			paste("Response:", responseName(mod)))
	result
}

# multivariate linear models

has.intercept.mlm <- function (model, ...) 
	any(row.names(coefficients(model)) == "(Intercept)")


Anova.mlm <- function(mod, type=c("II","III", 2, 3), SSPE, error.df, idata, 
		idesign, icontrasts=c("contr.sum", "contr.poly"), imatrix,
		test.statistic=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),...){
	wts <- if (!is.null(mod$weights)) mod$weights else rep(1, nrow(model.matrix(mod)))
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (missing(SSPE)) SSPE <- wcrossprod(residuals(mod), w=wts)
	if (missing(idata)) {
		idata <- NULL
		idesign <- NULL
	}
	if (missing(imatrix)) imatrix <- NULL
	error.df <- if (missing(error.df)) df.residual(mod)
			else error.df
	switch(type,
			II=Anova.II.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
			III=Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
			"2"=Anova.II.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
			"3"=Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...))
}

#Anova.mlm <- function(mod, type=c("II","III", 2, 3), SSPE, error.df, idata, 
#   idesign, icontrasts=c("contr.sum", "contr.poly"), imatrix,
#   test.statistic=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),...){
#   type <- as.character(type)
#   type <- match.arg(type)
#   test.statistic <- match.arg(test.statistic)
#   if (missing(SSPE)) SSPE <- crossprod(residuals(mod))
#   if (missing(idata)) {
#       idata <- NULL
#       idesign <- NULL
#   }
#   if (missing(imatrix)) imatrix <- NULL
#   error.df <- if (missing(error.df)) df.residual(mod)
#       else error.df
#   switch(type,
#       II=Anova.II.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
#       III=Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
#       "2"=Anova.II.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
#       "3"=Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...))
#}

Anova.III.mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...){
	intercept <- has.intercept(mod)
	V <- solve(crossprod(model.matrix(mod)))
	p <- nrow(coefficients(mod))
	I.p <- diag(p)
	terms <- term.names(mod)
	n.terms <- length(terms)
	assign <- mod$assign
	if (is.null(idata) && is.null(imatrix)){
		if ((n.terms == 0) && intercept) {
			Test <- linearHypothesis(mod, 1, SSPE=SSPE, ...)
			result <- list(SSP=Test$SSPH, SSPE=SSPE, df=1, error.df=error.df,
					terms="(Intercept)", repeated=FALSE, type="III", test=test)
			class(result) <- "Anova.mlm"
			return(result)
		}
		SSP <- as.list(rep(0, n.terms))
		df <- rep(0, n.terms)
		names(df) <- names(SSP) <- terms
		for (term in 1:n.terms){
			subs <- which(assign == term - intercept)
			hyp.matrix <- I.p[subs,,drop=FALSE]
			Test <- linearHypothesis(mod, hyp.matrix, SSPE=SSPE, ...)
			SSP[[term]] <- Test$SSPH
			df[term]<- length(subs)
		}
		result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
				repeated=FALSE, type="III", test=test)
	}
	else {
		if (!is.null(imatrix)){
			X.design <- do.call(cbind, imatrix)
			ncols <- sapply(imatrix, ncol)
			end <- cumsum(ncols)
			start <- c(1, (end + 1))[-(length(end) + 1)]
			cols <- mapply(seq, from=start, to=end)
			iterms <- names(end)
			names(cols) <- iterms
			check.imatrix(X.design, iterms)
		}
		else {
			if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
			for (i in 1:length(idata)){
				if (is.null(attr(idata[,i], "contrasts"))){
					contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
							else icontrasts[1]
				}
			}
			X.design <- model.matrix(idesign, data=idata)
			i.intercept <- has.intercept(X.design)
			iterms <- term.names(idesign)
			if (i.intercept) iterms <- c("(Intercept)", iterms)
			check.imatrix(X.design)
		}
		df <- rep(0, n.terms*length(iterms))
		hnames <- rep("", length(df))
		P <- SSPEH <- SSP <- as.list(df)
		i <- 0
		for (iterm in iterms){
			for (term in 1:n.terms){
				subs <- which(assign == term - intercept)
				hyp.matrix <- I.p[subs,,drop=FALSE]
				i <- i + 1
				Test <- linearHypothesis(mod, hyp.matrix, SSPE=SSPE, 
						idata=idata, idesign=idesign, icontrasts=icontrasts, iterms=iterm, 
						check.imatrix=FALSE, P=imatrix[[iterm]], ...)
				SSP[[i]] <- Test$SSPH
				SSPEH[[i]] <- Test$SSPE
				P[[i]] <- Test$P
				df[i] <- length(subs)
				hnames[i] <- if (iterm == "(Intercept)") terms[term]
						else if (terms[term] == "(Intercept)") iterm
						else paste(terms[term], ":", iterm, sep="")
			}
		}
		names(df) <- names(SSP) <- names(SSPEH) <- names(P) <- hnames
		result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
				terms=hnames, repeated=TRUE, type="III", test=test, 
				idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix)       
	}
	class(result) <- "Anova.mlm"
	result
}

Anova.II.mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...){
	wts <- if (!is.null(mod$weights)) mod$weights else rep(1, nrow(model.matrix(mod)))
	V <- solve(wcrossprod(model.matrix(mod), w=wts))
	SSP.term <- function(term, iterm){
		which.term <- which(term == terms)
		subs.term <- which(assign == which.term)
		relatives <- relatives(term, terms, fac)
		subs.relatives <- NULL
		for (relative in relatives) subs.relatives <- c(subs.relatives, which(assign==relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives, subs.term),,drop=FALSE]
		if (missing(iterm)){
			SSP1 <- if (length(subs.relatives) == 0) 0 
					else linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, ...)$SSPH
			SSP2 <- linearHypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, ...)$SSPH
			return(SSP2 - SSP1)
		}
		else {
			SSP1 <- if (length(subs.relatives) == 0) 0 
					else linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
								idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, P=imatrix[[iterm]], ...)$SSPH
			lh2 <- linearHypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, 
					idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, P=imatrix[[iterm]], ...)
			return(list(SSP = lh2$SSPH - SSP1, SSPE=lh2$SSPE, P=lh2$P))
		}
	}
	fac <- attr(mod$terms, "factors")
	intercept <- has.intercept(mod)
	p <- nrow(coefficients(mod))
	I.p <- diag(p)
	assign <- mod$assign
	terms <- term.names(mod)
	if (intercept) terms <- terms[-1]
	n.terms <- length(terms)
	if (n.terms == 0){
		message("Note: model has only an intercept; equivalent type-III tests substituted.")
		return(Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...))
	}
	if (is.null(idata) && is.null(imatrix)){
		SSP <- as.list(rep(0, n.terms))
		df <- rep(0, n.terms)
		names(df) <- names(SSP) <- terms
		for (i in 1:n.terms){
			SSP[[i]] <- SSP.term(terms[i])
			df[i]<- df.terms(mod, terms[i])
		}    
		result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
				repeated=FALSE, type="II", test=test)
	}
	else {
		if (!is.null(imatrix)){
			X.design <- do.call(cbind, imatrix)
			ncols <- sapply(imatrix, ncol)
			end <- cumsum(ncols)
			start <- c(1, (end + 1))[-(length(end) + 1)]
			cols <- mapply(seq, from=start, to=end)
			iterms <- names(end)
			names(cols) <- iterms
			check.imatrix(X.design, iterms)     
		}
		else {
			if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
			for (i in 1:length(idata)){
				if (is.null(attr(idata[,i], "contrasts"))){
					contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
							else icontrasts[1]
				}
			}
			X.design <- model.matrix(idesign, data=idata)
			iintercept <- has.intercept(X.design)
			iterms <- term.names(idesign)
			if (iintercept) iterms <- c("(Intercept)", iterms)
			check.imatrix(X.design)
		}
		df <- rep(0, (n.terms + intercept)*length(iterms))
		hnames <- rep("", length(df))
		P <- SSPEH <- SSP <- as.list(df)
		i <- 0
		for (iterm in iterms){
			if (intercept){
				i <- i + 1
				hyp.matrix.1 <- I.p[-1,,drop=FALSE]
				SSP1 <- linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
						idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, 
						check.imatrix=FALSE, P=imatrix[[iterm]],...)$SSPH
				lh2 <- linearHypothesis(mod, I.p, SSPE=SSPE, V=V, 
						idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, 
						check.imatrix=FALSE, P=imatrix[[iterm]], ...)
				SSP[[i]] <- lh2$SSPH - SSP1
				SSPEH[[i]] <- lh2$SSPE
				P[[i]] <- lh2$P
				df[i] <- 1
				hnames[i] <- iterm
			}
			for (term in 1:n.terms){
				subs <- which(assign == term)
				i <- i + 1
				Test <- SSP.term(terms[term], iterm)
				SSP[[i]] <- Test$SSP
				SSPEH[[i]] <- Test$SSPE
				P[[i]] <- Test$P
				df[i]<- length(subs)
				hnames[i] <- if (iterm == "(Intercept)") terms[term]
						else paste(terms[term], ":", iterm, sep="")
			}
		}
		names(df) <- names(P) <- names(SSP) <- names(SSPEH) <- hnames
		result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
				terms=hnames, repeated=TRUE, type="II", test=test,
				idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix)       
	}
	class(result) <- "Anova.mlm"
	result
}

#Anova.II.mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...){
#   V <- solve(crossprod(model.matrix(mod)))
#   SSP.term <- function(term, iterm){
#       which.term <- which(term == terms)
#       subs.term <- which(assign == which.term)
#       relatives <- relatives(term, terms, fac)
#       subs.relatives <- NULL
#       for (relative in relatives) subs.relatives <- c(subs.relatives, which(assign==relative))
#       hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
#       hyp.matrix.2 <- I.p[c(subs.relatives, subs.term),,drop=FALSE]
#       if (missing(iterm)){
#           SSP1 <- if (length(subs.relatives) == 0) 0 
#               else linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, ...)$SSPH
#           SSP2 <- linearHypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, ...)$SSPH
#           return(SSP2 - SSP1)
#       }
#       else {
#           SSP1 <- if (length(subs.relatives) == 0) 0 
#               else linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
#                       idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, P=imatrix[[iterm]], ...)$SSPH
#           lh2 <- linearHypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, 
#               idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, P=imatrix[[iterm]], ...)
#           return(list(SSP = lh2$SSPH - SSP1, SSPE=lh2$SSPE, P=lh2$P))
#       }
#   }
#   fac <- attr(mod$terms, "factors")
#   intercept <- has.intercept(mod)
#   p <- nrow(coefficients(mod))
#   I.p <- diag(p)
#   assign <- mod$assign
#   terms <- term.names(mod)
#   if (intercept) terms <- terms[-1]
#   n.terms <- length(terms)
#   if (n.terms == 0){
#       message("Note: model has only an intercept; equivalent type-III tests substituted.")
#       return(Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...))
#   }
#   if (is.null(idata) && is.null(imatrix)){
##      if ((n.terms == 0) && intercept) {
##          Test <- linearHypothesis(mod, 1, SSPE=SSPE, ...)
##          result <- list(SSP=list(Test$SSPH), SSPE=SSPE, df=1, error.df=error.df,
##              terms="(Intercept)", repeated=FALSE, type="II", test=test)
##          class(result) <- "Anova.mlm"
##          return(result)
##      }
#       SSP <- as.list(rep(0, n.terms))
#       df <- rep(0, n.terms)
#       names(df) <- names(SSP) <- terms
#       for (i in 1:n.terms){
#           SSP[[i]] <- SSP.term(terms[i])
#           df[i]<- df.terms(mod, terms[i])
#       }    
#       result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
#           repeated=FALSE, type="II", test=test)
#   }
#   else {
#       if (!is.null(imatrix)){
#           X.design <- do.call(cbind, imatrix)
#           ncols <- sapply(imatrix, ncol)
#           end <- cumsum(ncols)
#           start <- c(1, (end + 1))[-(length(end) + 1)]
#           cols <- mapply(seq, from=start, to=end)
#           iterms <- names(end)
#           names(cols) <- iterms
#           check.imatrix(X.design, iterms)     
#       }
#       else {
#           if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
#           for (i in 1:length(idata)){
#               if (is.null(attr(idata[,i], "contrasts"))){
#                   contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
#                       else icontrasts[1]
#               }
#           }
#           X.design <- model.matrix(idesign, data=idata)
#           iintercept <- has.intercept(X.design)
#           iterms <- term.names(idesign)
#           if (iintercept) iterms <- c("(Intercept)", iterms)
#           check.imatrix(X.design)
#       }
#       df <- rep(0, (n.terms + intercept)*length(iterms))
#       hnames <- rep("", length(df))
#       P <- SSPEH <- SSP <- as.list(df)
#       i <- 0
#       for (iterm in iterms){
#           if (intercept){
#               i <- i + 1
#               hyp.matrix.1 <- I.p[-1,,drop=FALSE]
#               SSP1 <- linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
#                   idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, 
#                   check.imatrix=FALSE, P=imatrix[[iterm]],...)$SSPH
#               lh2 <- linearHypothesis(mod, I.p, SSPE=SSPE, V=V, 
#                   idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, 
#                   check.imatrix=FALSE, P=imatrix[[iterm]], ...)
#               SSP[[i]] <- lh2$SSPH - SSP1
#               SSPEH[[i]] <- lh2$SSPE
#               P[[i]] <- lh2$P
#               df[i] <- 1
#               hnames[i] <- iterm
#           }
#           for (term in 1:n.terms){
#               subs <- which(assign == term)
#               i <- i + 1
#               Test <- SSP.term(terms[term], iterm)
#               SSP[[i]] <- Test$SSP
#               SSPEH[[i]] <- Test$SSPE
#               P[[i]] <- Test$P
#               df[i]<- length(subs)
#               hnames[i] <- if (iterm == "(Intercept)") terms[term]
#                   else paste(terms[term], ":", iterm, sep="")
#           }
#       }
#       names(df) <- names(P) <- names(SSP) <- names(SSPEH) <- hnames
#       result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
#           terms=hnames, repeated=TRUE, type="II", test=test,
#           idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix)       
#   }
#   class(result) <- "Anova.mlm"
#   result
#}

print.Anova.mlm <- function(x, ...){
	test <- x$test
	repeated <- x$repeated
	ntests <- length(x$terms)
	tests <- matrix(NA, ntests, 4)
	if (!repeated) SSPE.qr <- qr(x$SSPE) 
	for (term in 1:ntests){
		# some of the code here adapted from stats:::summary.manova
		eigs <- Re(eigen(qr.coef(if (repeated) qr(x$SSPE[[term]]) else SSPE.qr,
								x$SSP[[term]]), symmetric = FALSE)$values)
		tests[term, 1:4] <- switch(test,
				Pillai = stats:::Pillai(eigs, x$df[term], x$error.df),
				Wilks = stats:::Wilks(eigs, x$df[term], x$error.df),
				"Hotelling-Lawley" = stats:::HL(eigs, x$df[term], x$error.df),
				Roy = stats:::Roy(eigs, x$df[term], x$error.df))
	}
	ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
	ok <- !is.na(ok) & ok
	tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4], 
					lower.tail = FALSE))
	rownames(tests) <- x$terms
	colnames(tests) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")
	tests <- structure(as.data.frame(tests), 
			heading = paste("\nType ", x$type, if (repeated) " Repeated Measures",
					" MANOVA Tests: ", test, " test statistic", sep=""), 
			class = c("anova", "data.frame"))
	print(tests)      
	invisible(x)
}

summary.Anova.mlm <- function(object, test.statistic, multivariate=TRUE, univariate=TRUE, 
		digits=getOption("digits"), ...){
	GG <- function(SSPE, P){ # Greenhouse-Geisser correction
		p <- nrow(SSPE)
		if (p < 2) return(NA) 
		lambda <- eigen(SSPE %*% solve(t(P) %*% P))$values
		lambda <- lambda[lambda > 0]
		((sum(lambda)/p)^2)/(sum(lambda^2)/p)
	}
	HF <- function(gg, error.df, p){ # Huynh-Feldt correction
		((error.df + 1)*p*gg - 2)/(p*(error.df - p*gg))
	}
	mauchly <- function (SSD, P, df) {
		# most of this function borrowed from stats:::mauchly.test.SSD
		if (nrow(SSD) < 2) return(c(NA, NA))
		Tr <- function (X) sum(diag(X))
		p <- nrow(P)
		I <- diag(p)
		Psi <- t(P) %*% I %*% P 
		B <- SSD 
		pp <- nrow(SSD) 
		U <- solve(Psi, B)
		n <- df 
		logW <- log(det(U)) - pp * log(Tr(U/pp))
		rho <- 1 - (2 * pp^2 + pp + 2)/(6 * pp * n)
		w2 <- (pp + 2) * (pp - 1) * (pp - 2) * (2 * pp^3 + 6 * pp^2 + 
					3 * p + 2)/(288 * (n * pp * rho)^2)
		z <- -n * rho * logW
		f <- pp * (pp + 1)/2 - 1
		Pr1 <- pchisq(z, f, lower.tail = FALSE)
		Pr2 <- pchisq(z, f + 4, lower.tail = FALSE)
		pval <- Pr1 + w2 * (Pr2 - Pr1)
		c(statistic = c(W = exp(logW)), p.value = pval)
	}        
	if (missing(test.statistic)) test.statistic <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
	test.statistic <- match.arg(test.statistic, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
			several.ok=TRUE)
	nterms <- length(object$terms)
	if (multivariate || !object$repeated){       
		cat(paste("\nType ", object$type, if (object$repeated) " Repeated Measures",
						" MANOVA Tests:\n", sep=""))
		if (!object$repeated){ 
			cat("\nSum of squares and products for error:\n")
			print(object$SSPE, digits=digits)
		}
		for (term in 1:nterms){
			cat(paste("\n------------------------------------------\n",
							"\nTerm:", object$terms[term], "\n"))
			hyp <- list(SSPH=object$SSP[[term]], 
					SSPE=if (object$repeated) object$SSPE[[term]] else object$SSPE,
					P=if (object$repeated) object$P[[term]] else NULL, 
					test=test.statistic, df=object$df[term], 
					df.residual=object$error.df, title=object$terms[term])
			class(hyp) <- "linearHypothesis.mlm"
			print(hyp, digits=digits, SSPE=object$repeated, ...)
		}
	}
	if (object$repeated && univariate){
		error.df <- object$error.df
		table <- matrix(0, nterms, 6)
		table2 <- matrix(0, nterms, 4)
		table3 <- matrix(0, nterms, 2)
		rownames(table3) <- rownames(table2) <- rownames(table) <- object$terms
		colnames(table) <- c("SS", "num Df", "Error SS", "den Df", "F", "Pr(>F)")
		colnames(table2) <- c("GG eps", "Pr(>F[GG])",  "HF eps", "Pr(>F[HF])")
		colnames(table3) <- c("Test statistic", "p-value")
		for (term in 1:nterms){
			SSP <- object$SSP[[term]]
			SSPE <- object$SSPE[[term]]
			P <- object$P[[term]]
			p <- ncol(P)
			PtPinv <- solve(t(P) %*% P)
			gg <- GG(SSPE, P)
			table[term, "SS"] <- sum(diag(SSP %*% PtPinv))
			table[term, "Error SS"] <- sum(diag(SSPE %*% PtPinv))
			table[term, "num Df"] <- object$df[term] * p
			table[term, "den Df"] <- error.df * p
			table[term, "F"] <-  (table[term, "SS"]/table[term, "num Df"])/
					(table[term, "Error SS"]/table[term, "den Df"])
			table[term, "Pr(>F)"] <- pf(table[term, "F"], table[term, "num Df"],
					table[term, "den Df"], lower.tail=FALSE)
			table2[term, "GG eps"] <- gg
			table2[term, "HF eps"] <- HF(gg, error.df, p)
			table3[term,] <- mauchly(SSPE, P, object$error.df)
		}
		cat("\nUnivariate Type", object$type, 
				"Repeated-Measures ANOVA Assuming Sphericity\n\n")
		print.anova(table)
		table3 <- na.omit(table3)
		if (nrow(table3) > 0){
			cat("\n\nMauchly Tests for Sphericity\n\n")
			print.anova(table3)
			cat("\n\nGreenhouse-Geisser and Huynh-Feldt Corrections\n",
					"for Departure from Sphericity\n\n")
			table2[,"Pr(>F[GG])"] <- pf(table[,"F"], table2[,"GG eps"]*table[,"num Df"],
					table2[,"GG eps"]*table[,"den Df"], lower.tail=FALSE)
			table2[,"Pr(>F[HF])"] <- pf(table[,"F"], 
					pmin(1, table2[,"HF eps"])*table[,"num Df"],
					pmin(1, table2[,"HF eps"])*table[,"den Df"], lower.tail=FALSE)
			table2 <- na.omit(table2)
			print.anova(table2[,1:2, drop=FALSE])
			cat("\n")
			print.anova(table2[,3:4, drop=FALSE])
			if (any(table2[,"HF eps"] > 1)) 
				warning("HF eps > 1 treated as 1")
		}
	}
	invisible(object)
}

Anova.manova <- function(mod, ...){
	class(mod) <- c("mlm", "lm")
	Anova(mod, ...)
}

Manova <- function(mod, ...){
	UseMethod("Manova")
}

Manova.mlm <- function(mod, ...){
	Anova(mod, ...)
}

# Cox regression models

df.residual.coxph <- function(object, ...){    
	object$n - sum(!is.na(coef(object)))
}

alias.coxph <- function(model){
	if(any(which <- is.na(coef(model)))) return(list(Complete=which))
	else list()
}

logLik.coxph <- function(object, ...) object$loglik[2]

Anova.coxph <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald"), ...){
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (length((mod$rscore) > 0) && (test.statistic == "LR")){ 
		warning("LR tests unavailable with robust variances\nWald tests substituted")
		test.statistic <- "Wald"
	}
	switch(type,
			II=switch(test.statistic,
					LR=Anova.II.LR.coxph(mod),
					Wald=Anova.default(mod, type="II", test="Chisq", vcov.=vcov(mod))),
			III=switch(test.statistic,
					LR=Anova.III.LR.coxph(mod),
					Wald=Anova.default(mod, type="III", test="Chisq", vcov.=vcov(mod))),
			"2"=switch(test.statistic,
					LR=Anova.II.LR.coxph(mod),
					Wald=Anova.default(mod, type="II", test="Chisq", vcov.=vcov(mod))),
			"3"=switch(test.statistic,
					LR=Anova.III.LR.coxph(mod),
					Wald=Anova.default(mod, type="III", test="Chisq", vcov.=vcov(mod))))
}

Anova.II.LR.coxph <- function(mod, ...){
	which.nms <- function(name) which(asgn == which(names == name))
	fac <-attr(terms(mod), "factors")
	names <- term.names(mod)
	n.terms <- length(names)
	if (n.terms < 2) return(anova(mod, test="Chisq"))
	X <- model.matrix(mod)
	asgn <- attr(X, 'assign')
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	for (term in 1:n.terms){
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
		mod.1 <- coxph(mod$y ~ X[, -exclude.1, drop = FALSE])
		loglik.1 <- logLik(mod.1)
		mod.2 <- if (length(rels) == 0) mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					coxph(mod$y ~ X[, -exclude.2, drop = FALSE])
				}
		loglik.2 <- logLik(mod.2)
		LR[term] <- -2*(loglik.1 - loglik.2)
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- "Analysis of Deviance Table (Type II tests)"
	result
}

Anova.III.LR.coxph <- function(mod, ...){
	which.nms <- function(name) which(asgn == which(names == name))
	fac <-attr(terms(mod), "factors")
	names <- term.names(mod)
	n.terms <- length(names)
	if (n.terms < 2) return(anova(mod, test="Chisq"))
	X <- model.matrix(mod)
	asgn <- attr(X, 'assign')
	df <- df.terms(mod)
	LR <- p <- rep(0, n.terms)
	loglik1 <- logLik(mod)
	for (term in 1:n.terms){
		mod.0 <- coxph(mod$y ~ X[, -which.nms(names[term])])
		LR[term] <- -2*(logLik(mod.0) - loglik1)
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df","Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result,"heading") <- "Analysis of Deviance Table (Type III tests)"
	result
}

# parametric survival regression models

alias.survreg <- function(model){
	if(any(which <- diag(vcov(model)) < 1e-10)) return(list(Complete=which))
	else list()
}

logLik.survreg <- function(object, ...) object$loglik[2]

Anova.survreg <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald"), ...){
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (length((mod$rscore) > 0) && (test.statistic == "LR")){ 
		warning("LR tests unavailable with robust variances\nWald tests substituted")
		test.statistic <- "Wald"
	}
	switch(type,
			II=switch(test.statistic,
					LR=Anova.II.LR.survreg(mod),
					Wald=Anova.II.Wald.survreg(mod)),
			III=switch(test.statistic,
					LR=Anova.III.LR.survreg(mod),
					Wald=Anova.III.Wald.survreg(mod)),
			"2"=switch(test.statistic,
					LR=Anova.II.LR.survreg(mod),
					Wald=Anova.II.Wald.survreg(mod)),
			"3"=switch(test.statistic,
					LR=Anova.III.LR.survreg(mod),
					Wald=Anova.III.Wald.survreg(mod)))
}

Anova.II.LR.survreg <- function(mod, ...){
	which.nms <- function(name) which(asgn == which(names == name))
	fac <-attr(terms(mod), "factors")
	names <- term.names(mod)
	X <- model.matrix(mod)
	asgn <- attr(X, 'assign')
	asgn <- asgn[asgn != 0]
	if (has.intercept(mod)){
		int <- which(names == "(Intercept)")
		X <- X[, -int]
		names <- names[-int]
	}
	n.terms <- length(names)
	if (n.terms < 2) return(anova(mod))
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	y <- model.frame(mod)[,1]
	for (term in 1:n.terms){
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
		mod.1 <- survreg(y ~ X[, -exclude.1, drop = FALSE])
		loglik.1 <- logLik(mod.1)
		mod.2 <- if (length(rels) == 0) mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					survreg(y ~ X[, -exclude.2, drop = FALSE])
				}
		loglik.2 <- logLik(mod.2)
		LR[term] <- -2*(loglik.1 - loglik.2)
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- "Analysis of Deviance Table (Type II tests)"
	result
}

Anova.III.LR.survreg <- function(mod, ...){
	which.nms <- function(name) which(asgn == which(names == name))
	fac <-attr(terms(mod), "factors")
	names <- term.names(mod)
	X <- model.matrix(mod)
	asgn <- attr(X, 'assign')
	asgn <- asgn[asgn != 0]
	if (has.intercept(mod)){
		int <- which(names == "(Intercept)")
		X <- X[, -int]
		names <- names[-int]
	}
	n.terms <- length(names)
	if (n.terms < 2) return(anova(mod))
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	y <- model.frame(mod)[,1]
	loglik1 <- logLik(mod)
	for (term in 1:n.terms){
		mod.0 <- survreg(y ~ X[, -which.nms(names[term])])
		LR[term] <- -2*(logLik(mod.0) - loglik1)
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df","Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result,"heading") <- "Analysis of Deviance Table (Type III tests)"
	result
}

Anova.II.Wald.survreg <- function(mod){
	V <- vcov(mod)
	p <- nrow(V)
	V <- V[-p, -p]
	Anova.II.default(mod, V, test="Chisq")
}

Anova.III.Wald.survreg <- function(mod){
	V <- vcov(mod)
	p <- nrow(V)
	V <- V[-p, -p]
	Anova.III.default(mod, V, test="Chisq")
}

# Default Anova() method: requires methods for vcov() (if vcov. argument not specified) and coef().

Anova.default <- function(mod, type=c("II","III", 2, 3), test.statistic=c("Chisq", "F"), 
		vcov.=vcov(mod), singular.ok, ...){
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (missing(singular.ok))
		singular.ok <- type == "2" || type == "II"
	switch(type,
			II=Anova.II.default(mod, vcov., test.statistic, singular.ok=singular.ok),
			III=Anova.III.default(mod, vcov., test.statistic, singular.ok=singular.ok),
			"2"=Anova.II.default(mod, vcov., test.statistic, singular.ok=singular.ok),
			"3"=Anova.III.default(mod, vcov., test.statistic, singular.ok=singular.ok))
}

Anova.II.default <- function(mod, vcov., test, singular.ok=TRUE, ...){
	hyp.term <- function(term){
		which.term <- which(term==names)
		subs.term <- which(assign==which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives) 
			subs.relatives <- c(subs.relatives, which(assign==relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(statistic=NA, df=0))            
		hyp <- linearHypothesis.default(mod, hyp.matrix.term, 
				vcov.=vcov., test=test, singular.ok=singular.ok, ...)
		if (test=="Chisq") c(statistic=hyp$Chisq[2], df=hyp$Df[2])
		else c(statistic=hyp$F[2], df=hyp$Df[2])
	}
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(mod$terms, "factors")
	intercept <- has.intercept(mod)
	p <- length(coefficients(mod))
	I.p <- diag(p)
	assign <- attr(model.matrix(mod), "assign")
	assign[!not.aliased] <- NA
	names <- term.names(mod)
	if (intercept) names <- names[-1]
	n.terms <- length(names)
	df <- c(rep(0, n.terms), df.residual(mod))
	if (inherits(mod, "coxph")){
		assign <- assign[assign != 0]
		clusters <- grep("^cluster\\(", names)
		if (length(clusters) > 0) {
			names <- names[-clusters]
			df <- df[-clusters]
		}
	}
	p <- teststat <- rep(0, n.terms + 1)
	teststat[n.terms + 1] <- p[n.terms + 1] <- NA
	for (i in 1:n.terms){
		hyp <- hyp.term(names[i])
		teststat[i] <- abs(hyp["statistic"])
		df[i] <- abs(hyp["df"])
		p[i] <- if (test == "Chisq") 
					pchisq(teststat[i], df[i], lower.tail=FALSE) 
				else pf(teststat[i], df[i], df[n.terms + 1], lower.tail=FALSE)
	}    
	result <- data.frame(df, teststat, p)
	row.names(result) <- c(names,"Residuals")
	names(result) <- c ("Df", test, if (test == "Chisq") "Pr(>Chisq)" 
					else "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

Anova.III.default <- function(mod, vcov., test, singular.ok=FALSE, ...){
	intercept <- has.intercept(mod)
	p <- length(coefficients(mod))
	I.p <- diag(p)
	names <- term.names(mod)
	n.terms <- length(names)
	assign <- attr(model.matrix(mod), "assign")
	df <- c(rep(0, n.terms), df.residual(mod))
	if (inherits(mod, "coxph")){
		if (intercept) names <- names[-1]
		assign <- assign[assign != 0]
		clusters <- grep("^cluster\\(", names)
		if (length(clusters) > 0) {
			names <- names[-clusters]
			df <- df[-clusters]
		}
	}
	if (intercept) df[1] <- 1
	teststat <- rep(0, n.terms + 1)
	p <- rep(0, n.terms + 1)
	teststat[n.terms + 1] <- p[n.terms + 1] <- NA
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	for (term in 1:n.terms){
		subs <- which(assign == term - intercept)        
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
		if (nrow(hyp.matrix) == 0){
			teststat[term] <- NA
			df[term] <- 0
			p[term] <- NA
		}
		else {
			hyp <- linearHypothesis.default(mod, hyp.matrix, 
					vcov.=vcov., test=test, singular.ok=singular.ok, ...)
			teststat[term] <- if (test=="Chisq") hyp$Chisq[2] else hyp$F[2]
			df[term] <- abs(hyp$Df[2])
			p[term] <- if (test == "Chisq") 
						pchisq(teststat[term], df[term], lower.tail=FALSE) 
					else pf(teststat[term], df[term], df[n.terms + 1], lower.tail=FALSE)
		}
	}
	result <- data.frame(df, teststat, p)
	row.names(result) <- c(names, "Residuals")
	names(result) <- c ("Df", test, if (test == "Chisq") "Pr(>Chisq)" 
					else "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

## functions for mixed models

# the following function, not exported, to make nlme play better with lme4

fixef <- function (object){
	if (isS4(object)) object@fixef else object$coefficients$fixed
}

Anova.mer <- function(mod, type=c("II","III", 2, 3),
		vcov.=vcov(mod), singular.ok, ...){
	type <- as.character(type)
	type <- match.arg(type)
	if (missing(singular.ok))
		singular.ok <- type == "2" || type == "II"
	switch(type,
			II=Anova.II.mer(mod, vcov., singular.ok=singular.ok),
			III=Anova.III.mer(mod, vcov., singular.ok=singular.ok),
			"2"=Anova.II.mer(mod, vcov., singular.ok=singular.ok),
			"3"=Anova.III.mer(mod, vcov., singular.ok=singular.ok))
}

Anova.II.mer <- function(mod, vcov., singular.ok=TRUE, ...){
	hyp.term <- function(term){
		which.term <- which(term==names)
		subs.term <- which(assign==which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives) 
			subs.relatives <- c(subs.relatives, which(assign==relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(statistic=NA, df=0))            
		hyp <- linearHypothesis(mod, hyp.matrix.term, 
				vcov.=vcov., singular.ok=singular.ok, ...)
		c(statistic=hyp$Chisq[2], df=hyp$Df[2])
	}
	not.aliased <- !is.na(fixef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(terms(mod), "factors")
	intercept <- has.intercept(mod)
	p <- length(fixef(mod))
	I.p <- diag(p)
	assign <- attr(model.matrix(mod), "assign")
	assign[!not.aliased] <- NA
	names <- term.names(mod)
	if (intercept) names <- names[-1]
	n.terms <- length(names)
	p <- teststat <- df <- rep(0, n.terms)
	for (i in 1:n.terms){
		hyp <- hyp.term(names[i])
		teststat[i] <- abs(hyp["statistic"])
		df[i] <- abs(hyp["df"])
		p[i] <- pchisq(teststat[i], df[i], lower.tail=FALSE) 
	}    
	result <- data.frame(teststat, df, p)
	row.names(result) <- names
	names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

Anova.III.mer <- function(mod, vcov., singular.ok=FALSE, ...){
	intercept <- has.intercept(mod)
	p <- length(fixef(mod))
	I.p <- diag(p)
	names <- term.names(mod)
	n.terms <- length(names)
	assign <- attr(model.matrix(mod), "assign")
	df <- rep(0, n.terms)
	if (intercept) df[1] <- 1
	p <- teststat <- rep(0, n.terms)
	not.aliased <- !is.na(fixef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	for (term in 1:n.terms){
		subs <- which(assign == term - intercept)        
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
		if (nrow(hyp.matrix) == 0){
			teststat[term] <- NA
			df[term] <- 0
			p[term] <- NA
		}
		else {
			hyp <- linearHypothesis(mod, hyp.matrix, 
					vcov.=vcov., singular.ok=singular.ok, ...)
			teststat[term] <-  hyp$Chisq[2] 
			df[term] <- abs(hyp$Df[2])
			p[term] <- pchisq(teststat[term], df[term], lower.tail=FALSE) 
		}
	}
	result <- data.frame(teststat, df, p)
	row.names(result) <- names
	names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

Anova.lme <- function(mod, type=c("II","III", 2, 3),
		vcov.=vcov(mod), singular.ok, ...){
	type <- as.character(type)
	type <- match.arg(type)
	if (missing(singular.ok))
		singular.ok <- type == "2" || type == "II"
	switch(type,
			II=Anova.II.lme(mod, vcov., singular.ok=singular.ok),
			III=Anova.III.lme(mod, vcov., singular.ok=singular.ok),
			"2"=Anova.II.lme(mod, vcov., singular.ok=singular.ok),
			"3"=Anova.III.lme(mod, vcov., singular.ok=singular.ok))
}

Anova.II.lme <- function(mod, vcov., singular.ok=TRUE, ...){
	hyp.term <- function(term){
		which.term <- which(term==names)
		subs.term <- which(assign==which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives) 
			subs.relatives <- c(subs.relatives, which(assign==relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(statistic=NA, df=0))            
		hyp <- linearHypothesis(mod, hyp.matrix.term, 
				vcov.=vcov., singular.ok=singular.ok, ...)
		c(statistic=hyp$Chisq[2], df=hyp$Df[2])
	}
	not.aliased <- !is.na(fixef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(terms(mod), "factors")
	intercept <- has.intercept(mod)
	p <- length(fixef(mod))
	I.p <- diag(p)
	assign <- attr(model.matrix(mod), "assign")
	assign[!not.aliased] <- NA
	names <- term.names(mod)
	if (intercept) names <- names[-1]
	n.terms <- length(names)
	p <- teststat <- df <- rep(0, n.terms)
	for (i in 1:n.terms){
		hyp <- hyp.term(names[i])
		teststat[i] <- abs(hyp["statistic"])
		df[i] <- abs(hyp["df"])
		p[i] <- pchisq(teststat[i], df[i], lower.tail=FALSE) 
	}    
	result <- data.frame(teststat, df, p)
	row.names(result) <- names
	names(result) <- c("Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

Anova.III.lme <- function(mod, vcov., singular.ok=FALSE, ...){
	intercept <- has.intercept(mod)
	p <- length(coefficients(mod))
	I.p <- diag(p)
	names <- term.names(mod)
	n.terms <- length(names)
	assign <- attr(model.matrix(mod), "assign")
	df <- rep(0, n.terms)
	if (intercept) df[1] <- 1
	p <- teststat <-rep(0, n.terms)
	not.aliased <- !is.na(fixef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	for (term in 1:n.terms){
		subs <- which(assign == term - intercept)        
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
		if (nrow(hyp.matrix) == 0){
			teststat[term] <- NA
			df[term] <- 0
			p[term] <- NA
		}
		else {
			hyp <- linearHypothesis(mod, hyp.matrix, 
					vcov.=vcov., singular.ok=singular.ok, ...)
			teststat[term] <-  hyp$Chisq[2] 
			df[term] <- abs(hyp$Df[2])
			p[term] <- pchisq(teststat[term], df[term], lower.tail=FALSE) 
		}
	}
	result <- data.frame(teststat, df, p)
	row.names(result) <- names
	names(result) <- c ("Chisq",  "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

Anova.svyglm <- function(mod, ...) Anova.default(mod, ...)
