# last modified 8 March 2008 by J. Fox

# moved from Rcmdr 5 December 2006

# the following function adapted from Fox, An R and S-PLUS Companion to Applied Regression
influencePlot <- function(model, ...){
    UseMethod("influencePlot")
    }

	influencePlot.lm <- function(model, scale=10, col=c(1,2), identify=c(TRUE, FALSE, "auto"),
		labels=names(rstud), cex.identify=par("cex"), col.identify=par("col"), ...){
		identify <- as.character(identify)
		identify <- identify[1]
		if(!identify %in% c(TRUE, FALSE, "auto")) stop('argument identify must be one of TRUE, FALSE, "auto"')
		hatval <- hatvalues(model)
		rstud <- rstudent(model)
		cook <- sqrt(cookd(model))
		scale <- scale/max(cook, na.rm=TRUE)
		p <- length(coef(model))
		n <- length(rstud)
		cutoff <- sqrt(4/(n - p))
		plot(hatval, rstud, xlab='Hat-Values',
			ylab='Studentized Residuals', type='n', ...)
		abline(v=c(2, 3)*p/n, lty=2)
		abline(h=c(-2, 0, 2), lty=2)
		points(hatval, rstud, cex=scale*cook, 
			col=ifelse(noteworthy <- cook > cutoff, col[2], col[1]))
		if (identify == "TRUE") identify(hatval, rstud, labels, col=col.identify,
				cex=cex.identify)
		else if (identify == "auto"){
			pos <- ifelse((hatval - sum(rev(range(hatval)))/2) <= 0, 4, 2)
			text(hatval[noteworthy], rstud[noteworthy], labels[noteworthy],  pos=pos[noteworthy],
				cex=cex.identify, col=col.identify)
		}
	}
