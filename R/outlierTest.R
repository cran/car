#-------------------------------------------------------------------------------
# Revision history:
# 2009-09-28 by J. Fox (renamed)
# 2010-04-14 by J. Fox fixed error in reporting largest abs rstudent
# 2012-12-12 by J. Fox fixed handling of labels argument
# 2019-01-02 by J. Fox added lmerMod method
# 2019-05-12 by J. Fox fixed spelling of "Bonferroni"
#-------------------------------------------------------------------------------

# Bonferroni test for an outlier (J. Fox)

outlierTest <- function(model, ...){
	UseMethod("outlierTest")
}

outlierTest.lm <- function(model, cutoff=0.05, n.max=10, order=TRUE, labels=names(rstudent), ...){
	rstudent <- rstudent(model)
	if (length(rstudent) != length(labels))
		stop("Number of labels does not correspond to number of residuals.")
    else names(rstudent) <- labels
	df <- df.residual(model) - 1
	rstudent <- rstudent[!is.na(rstudent)]
	n <- length(rstudent)
	p <- if (class(model)[1] == "glm")
			2*(pnorm(abs(rstudent), lower.tail=FALSE))
		else 2*(pt(abs(rstudent), df, lower.tail=FALSE))
	bp <- n*p
	ord <- if (order) order(bp) else 1:n
	ord <- ord[bp[ord] <= cutoff]
	result <- if (length(ord) == 0){
			which <- which.max(abs(rstudent))
			list(rstudent=rstudent[which], p=p[which], bonf.p=bp[which], signif=FALSE, cutoff=cutoff)
		}
		else {
			if (length(ord) > n.max) ord <- ord[1:n.max]
			result <- list(rstudent=rstudent[ord], p=p[ord], bonf.p=bp[ord], signif=TRUE, cutoff=cutoff)
		}
	class(result)<-"outlierTest"
	result
}

outlierTest.lmerMod <- function(model, ...){
  outlierTest.lm(model, ...)
}

print.outlierTest<-function(x, digits=5, ...){
	if (!x$signif){
		cat("No Studentized residuals with Bonferroni p <", x$cutoff)
		cat("\nLargest |rstudent|:\n")
	}
	bp <- x$bonf
	bp[bp > 1] <- NA
	table <- data.frame(rstudent=x$rstudent,
		"unadjusted p-value"=signif(x$p, digits), "Bonferroni p"=signif(bp, digits),
		check.names=FALSE)
	rownames(table) <- names(x$rstudent)
	print(table)
	invisible(x)
}
