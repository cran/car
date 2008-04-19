# moved from Rcmdr 13 July 2004
# last modified 19 April 2008 by J. Fox

# levene.test.default function slightly modified from Brian Ripley via R-help
# the generic version here was contributed by Derek Ogle

levene.test <- function (y, ...) {
  UseMethod("levene.test") 
}

levene.test.default <- function (y, group, ...) { # original levene.test
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


levene.test.formula <- function(y, ...) {
  form <- y
  mf <- model.frame(form, ...)
  if (any(sapply(2:dim(mf)[2], function(j) is.numeric(mf[[j]])))) stop("Levene's test is not appropriate with quantitative explanatory variables.")
  y <- mf[,1]
  if(dim(mf)[2]==2) group <- mf[,2]
    else {
      if (length(grep("\\+ | \\| | \\^ | \\:",form))>0) stop("Model must be completely crossed formula only.")
      group <- interaction(mf[,2:dim(mf)[2]])
  }
  levene.test.default(y=y,group=group)
}


levene.test.lm <- function(y, ...) {
  levene.test.formula(formula(y), data=model.frame(y))
}
