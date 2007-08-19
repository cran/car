# moved from Rcmdr 13 July 2004
# last modified 26 March 2007 by J. Fox

# the following function slightly modified from Brian Ripley via R-help

levene.test <- function (y, group) {
    if (!is.numeric(y)) 
        stop(deparse(substitute(y)), " is not a numeric variable")
    if (!is.factor(group)) {
        warning(deparse(substitute(group)), " coerced to factor.")
        group <- as.factor(group)
        }
    meds <- tapply(y, group, median, na.rm = TRUE)
    resp <- abs(y - meds[group])
    table <- anova(lm(resp ~ group))[, c(1, 4, 5)]
    rownames(table)[2] <- " "
    attr(table, "heading") <- "Levene's Test for Homogeneity of Variance"
    table 
    }
