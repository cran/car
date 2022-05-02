# logit transformation of proportion or percent (J. Fox)

# 2022-02-04: Fix range tests (bug reported by Collin Erickson)
#             Print message if percents converted to proportions

logit <- function(p, percents, adjust){
	range.p <- range(p, na.rm=TRUE)
	if (missing(percents) && range.p[2] > 1){
	  percents <- TRUE
	  message("Note: largest value of p > 1 so values of p interpreted as percents")
	} else {
	  percents <- FALSE
	}
	if (percents){
		if (range.p[1] < 0 || range.p[2] > 100) stop("p must be in the range 0 to 100")
		p <- p/100
		range.p <- range.p/100
	}
	else if (range.p[1] < 0 || range.p[2] > 1) stop("p must be in the range 0 to 1")
	a <-if (missing(adjust)) {
				if (isTRUE(all.equal(range.p[1], 0)) || isTRUE(all.equal(range.p[2], 1))) .025 else 0
			}
			else adjust
	if (missing(adjust) && a != 0) warning(paste("proportions remapped to (", a, ", ", 1-a, ")", sep=""))
	a <- 1 - 2*a
	log((0.50 + a*(p - 0.50))/(1 - (0.50 + a*(p - 0.50))))
}
    
