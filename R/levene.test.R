# moved from Rcmdr 27 April 2004

# the following function slightly modified from Brian Ripley via R-help

levene.test <- function(y, group) {
    meds <- tapply(y, group, median, na.rm=TRUE)
    resp <- abs(y - meds[group])
    table <- anova(lm(resp ~ group))
    rownames(table)[2] <- " "
    cat("Levene's Test for Homogeneity of Variance\n\n")
    table[,c(1,4,5)]
    } 
