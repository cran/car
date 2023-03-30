# Modified Nov. 24, 2009 by S. Weisberg to use showLabels
#   rather than showExtremes
# 11 & 20 January 2010: changed lty=3 to lty=1 for fitted curve. J. Fox
# 14 April 2010: set id.n = 0. J. Fox
# 15 April 2010; rewrite showLabels
# 25 May 2010 added grid() to plots, S. Weisberg
# 15 August 2010, fixed so col= works correctly with plot, but not Boxplot
# 15 August 2010, deleted pch= argument, as it wasn't used
# 17 January 2011, allow spline terms; plot against
#   predict(model, type="terms")[[term.name]]
# 1 February 2011 default for AsIs changed to TRUE
# 31 March 2011 tukeyNonaddTest updated to check that yhat^2 is not 
#   a linear combination of other predictors (as in 1-way anova).
# 6 April 2011 omit printing lack-of-fit if no lack-of-fit test is possible
# 16 June 2011 allow layout=NA, in which case the layout is not set in this
#  function, so it is the responsibility of the user
# 10 Feb 2013:  adjusted colinearity check in tukeyNonaddTest
# 21 March 2013:  fixed nonconstant variance test with missing values for glms
# 11 July 2013:  wording changes
# 11 July 2013:  'groups' arg for residualPlot and residualPlots.
# 19 July 2014:  type='rstudent' fixed
# 7 October 2014: trapped error resulting from groups= when n<3
# 25 April 2016: checks for na.action=na.exclude and changes it to na.omit for compatibility with Rcmdr. sw
# 2017-02-13: consolidated id and smooth arguments. John
# 2017-11-30: substitute carPalette() for palette(). J. Fox
# 2019-11-14: change class(x) == "y" to inherits(x, "y")
# 2018-08-06: enabled spread and var for smoothers. J. Fox
# 2022-11-03: convert one-column matrix regressor to vector. J. Fox
# 2023-03-10: the 'lwd' argument was ignored, now it is used.  S. Weisberg
# 2023-03-10: the 'lty' argument was removed, as it was previously unused.  S. Weisberg
# 2023-03-14: grouping now works correctly.  the groups variable if used must be a quoted string
# 2023-03-15: AsIs argument removed

residualPlots <- function(model, ...){UseMethod("residualPlots")}


residualPlots.default <- function(model, terms= ~ . , 
                                  layout=NULL, ask, main="", 
                                  fitted=TRUE, AsIs=TRUE, plot=TRUE, tests=TRUE, groups, ...){
  decodeGroups <- function(terms){
    rhs <- terms[[2]]
    # is '|' present in the formula?
    if("|" %in% all.names(rhs)){
      if(length(rhs[[3]]) > 1) stop("only one conditional variable permitted")
      groups <- as.character(rhs[[3]])
      terms <- as.formula(paste("~", deparse(rhs[[2]])))} 
    else {
      groups <- NULL
      terms = terms}
    list(terms=terms, groups=groups)
  }
  gps <- decodeGroups(terms)
  if(!is.null(gps$groups)) {
    terms <- gps$terms
    groups <- gps$groups}
  if(missing(groups)) groups <- NULL
# Added for compatibility with Rcmdr
    if(inherits(model$na.action, "exclude")) model <- update(model, na.action=na.omit)
# End addition
# construct list of names of horizontal terms
  model.terms <- names(update(model, method="model.frame"))[-1] # remove response
  use.terms <- names(update(model, terms, method="model.frame"))[-1] # remove response
  offsets <- which(substr(use.terms, 1, 7) == "offset(") 
  if(any(offsets)) use.terms <- use.terms[-offsets] # remove offsets
  wts <- which(use.terms == "(weights)")
  if(any(wts)) use.terms <- use.terms[-wts] # remove weights
  AsIss <- which(substr(use.terms, 1, 2) == "I(" ) #)
  if(any(AsIss) & !AsIs) use.terms <- use.terms[-AsIss] # remove I() is AsIs=FALSE
  if( !all(use.terms %in% model.terms))
        stop("Only terms in the formula can be plotted.")
  nt <- length(use.terms) + fitted
  nr <- 0  
  if (nt == 0) stop("No plots specified")
  if (nt > 1 & plot == TRUE & (is.null(layout) || is.numeric(layout))) {
        if(is.null(layout)){
            layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2), 
                             c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))
        }
        ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
        op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
                  oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
        on.exit(par(op))
    }
    ans <- NULL   
    for (term in use.terms){
            nr <- nr + 1
            qtest <- if(is.null(groups))
                residualPlot(model, term, plot=plot, ...) else
                    residualPlot(model, term, plot=plot, groups=groups, ...)
            if(!is.null(qtest)){
                ans <- rbind(ans, qtest)
                row.names(ans)[nr] <- term}
        }    
    # Tukey's test
    if (fitted == TRUE){      
        tuk <- if(is.null(groups))
            residualPlot(model, "fitted", plot=plot, ...) else
                residualPlot(model, "fitted", plot=plot, groups=groups, ...)
        if (!is.null(tuk)  & class(model)[1] == "lm"){
            ans <- rbind(ans, tuk)
            row.names(ans)[nr + 1] <- "Tukey test"
            ans[nr + 1, 2] <- 2*pnorm(abs(ans[nr + 1, 1]), lower.tail=FALSE)}} 
    if(plot == TRUE) mtext(side=3, outer=TRUE, main, cex=1.2)
    if(!is.null(ans)) {
        dimnames(ans)[[2]] <- c("Test stat", "Pr(>|Test stat|)")
        return(if(tests == FALSE | !is.null(groups)) invisible(ans) else
            if(all(is.na(ans))) warning("No possible lack-of-fit tests") else 
                printCoefmat(ans, has.Pvalue=TRUE, na.print="")) } else
                    invisible(NULL)
}

residualPlots.lm <- function(model, ...) {
    residualPlots.default(model, ...)
}

residualPlots.glm <- function(model, ...) {
    residualPlots.default(model, ...)
}

residualPlot <- function(model, ...) UseMethod("residualPlot")

residualPlot.default <- function(model, variable = "fitted", type = "pearson",
                                 groups, 
                                 plot = TRUE,
                                 linear = TRUE,     
                                 quadratic = if(missing(groups)) TRUE else FALSE, 
                                 smooth=FALSE, id=FALSE,
                                 col = carPalette()[1], col.quad = carPalette()[2],
                                 pch=1,
                                 xlab, ylab, lwd=1, 
                                 grid=TRUE, key=!missing(groups), ...) {
    id <- applyDefaults(id, defaults=list(method="r", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
    if (isFALSE(id)){
        id.n <- 0
        id.method <- "none"
        labels <- id.cex <- id.col <- id.location <- NULL
    }
    else{
        labels <- id$labels
        if (is.null(labels)) labels <- names(na.omit(residuals(model)))
        id.method <- id$method
        id.n <- if ("identify" %in% id.method) Inf else id$n
        id.cex <- id$cex
        id.col <- id$col
        id.location <- id$location
    }
    
    smoother.args <- applyDefaults(smooth, defaults=list(smoother=loessLine, span=2/3, 
                                                         var=FALSE, col=carPalette()[3]), type="smooth")
    if (!isFALSE(smoother.args)) {
        smoother <- smoother.args$smoother 
        col.smooth <- smoother.args$col
        smoother.args$smoother <- smoother.args$col <- NULL
        if (is.null(smoother.args$spread)) smoother.args$spread <- smoother.args$var
    }
    else smoother <- NULL
    string.capitalize <- function(string) {
        paste(toupper(substring(string, 1, 1)), substring(string, 2), sep="")}
    # if(missing(labels)) 
    #     labels <- names(residuals(model)[!is.na(residuals(model))])
    ylab <- if(!missing(ylab)) ylab else
        paste(string.capitalize(type), "residuals")
    nb <- function(x) gsub(" ", "", x) # deletes all blanks from a column name
    column <- match(nb(variable), nb(names(model$model)))
# the groups argument should be a quoted name from names(model$model), corresponding usually to
# a factor, a character variable, or a numeric variable with few levels
    if(!missing(groups)){
      groups.name <- groups
      groups <- model$model[[match(nb(groups), nb(names(model$model)))]]
      if(is.null(groups)){stop("groups must be  in the model")}
      groups <- if(class(groups)[1] == "factor") groups else factor(groups, ordered=FALSE)
      if(length(levels(groups)) > length(groups)/4) stop("groups has too many distinct values")
      if(key){ 
        mar3 <- 1.1 + length(levels(groups))
        op <- par(mar=c(5.1, 4.1, mar3, 2.1))
        on.exit(par(op))
      }   
      colors <- if(length(col) >=length(levels(groups))) col else carPalette()
      col <- colors[as.numeric(groups)]
      pchs <- if(length(pch) >= length(levels(groups))) pch else 1:length(levels(groups))
      pch <-  pchs[as.numeric(groups)] 
    }
# Added for compatibility with Rcmdr
    if(inherits(model$na.action, "exclude")) model <- update(model, na.action=na.omit)
# End addition
    if(is.na(column) && variable != "fitted")
        stop(paste(variable, "is not a regressor in the mean function"))
    horiz <- if(variable == "fitted") predict(model) else model$model[[column]]
    lab <- if(variable == "fitted") {
        if(inherits(model, "glm")) 
            "Linear Predictor" else "Fitted values"} else variable
    lab <- if(!missing(xlab)) xlab else lab
    if(class(horiz)[1] == "ordered") horiz <- factor(horiz, ordered=FALSE)
    if(is.character((horiz))) horiz <- factor(horiz)
    if (is.matrix(horiz) && ncol(horiz) == 1) horiz <- as.vector(horiz)
    ans <-
        if(inherits(horiz, "poly")) {
            horiz <- horiz[ , 1]
            lab <- paste("Linear part of", lab)
            c(NA, NA)}
    else if (inherits(horiz, "matrix")) {
        horiz <- try(predict(model, type="terms"), silent=TRUE)
        if(inherits(horiz, "try-error"))
            stop("Could not plot spline terms") 
        warning("Splines replaced by a fitted linear combination")
        horiz <- horiz[ , variable]
        c(NA, NA)
    }
    else if (inherits(horiz, "factor")) c(NA, NA)
    else residCurvTest(model, variable)
    theResiduals <- switch(type, "rstudent"=rstudent(model), 
                           "rstandard"=rstandard(model), residuals(model, type=type))
    if(plot==TRUE){
        if(inherits(horiz, "factor")) {
            idm <- if(is.list(id.method)) {
                lapply(id.method, function(x) if(x[1]=="xy") "y" else x)} else {
                    if(id.method[1] == "xy") "y"}    
            Boxplot(theResiduals, horiz, xlab=lab, ylab=ylab, labels=labels, 
                    id.method=idm, id.n=id.n, id.cex=id.cex,  
                    id.col=id.col, id.location=id.location, ...) 
            abline(h=0, lty=2) } else 
            {    
                plot(horiz, theResiduals, xlab=lab, ylab=ylab, type="n", ...)
                if(grid){
                    grid(lty=1, equilogs=FALSE)
                    box()}
                points(horiz, theResiduals, col=col, pch=pch, ...)
                if(linear){
                    if(missing(groups)){abline(h=0, lty=2, lwd=1)} else {
                        for (g in 1:length(levels(groups)))
                            try(abline(lm(theResiduals ~ horiz, 
                                          subset=groups==levels(groups)[g]), lty=2, lwd=1,
                                       col=colors[g]), silent=TRUE)
                    }}
                if(quadratic){
                    new <- seq(min(horiz), max(horiz), length=200)
                    if(missing(groups)){
                        if(length(unique(horiz)) > 2){
                            lm2 <- lm(theResiduals ~ poly(horiz, 2))
                            lines(new, predict(lm2, list(horiz=new)), lty=1, lwd=lwd, col=col.quad)
                        }} else {
                            for (g in 1:length(levels(groups))){
                                if(length(unique(horiz)) > 2){
                                    lm2 <- lm(theResiduals~poly(horiz, 2),
                                              subset=groups==levels(groups)[g])
                                    lines(new, predict(lm2, list(horiz=new)), 
                                          lty=1, lwd=lwd, col=colors[g])
                                }}}}
                if(is.function(smoother))
                    if(missing(groups)){
                        smoother(horiz, theResiduals, col.smooth, log.x=FALSE, log.y=FALSE,
                                 spread=smoother.args$spread, smoother.args=smoother.args)} else
                                     for (g in 1:length(levels(groups))){
                                         sel <- groups == levels(groups)[g]
                                         smoother(horiz[sel], theResiduals[sel], colors[g], 
                                                  log.x=FALSE, log.y=FALSE,
                                                  spread=smoother.args$spread, 
                                                  smoother.args=smoother.args)}
                if(key & !missing(groups)){
                    items <- paste(groups.name, levels(groups), sep= " = ")
                    plotArrayLegend("top", items=items, col.items=colors, pch=pchs)
                }
                showLabels(horiz, theResiduals, labels=labels, 
                           method=id.method, n=id.n, cex=id.cex, 
                           col=id.col, location=id.location)  
            }
    }  
    invisible(ans)}

residCurvTest <- function(model, variable) {UseMethod("residCurvTest")}
residCurvTest.lm <- function(model, variable) {
    if(variable == "fitted") tukeyNonaddTest(model) else {
        if(is.na(match(variable, attr(model$terms, "term.labels"))))
            stop(paste(variable, "is not a term in the mean function")) else {
                xsqres <- qr.resid(model$qr, model.frame(model)[[variable]]^2)
                r <- residuals(model, type="pearson")
                m1 <- lm(r ~ xsqres, weights=weights(model))
                df.correction <- sqrt((df.residual(model)-1) / df.residual(m1))
                test <- summary(m1)$coef[2, 3] * df.correction
                c(Test=test, Pvalue=2 * pt(-abs(test), df.residual(model)-1))
            }}}

residCurvTest.glm <- function(model, variable) {
    if(variable == "fitted") c(NA, NA) else {
        if(is.na(match(variable, attr(model$terms, "term.labels"))))
            stop(paste(variable, "is not a term in the mean function")) else {
                newmod <- paste(" ~ . + I(", variable, "^2)")
                m2 <- update(model, newmod, start=NULL)
                c(Test= test<-deviance(model)-deviance(m2), Pvalue=1-pchisq(test, 1))
            }}}

residCurvTest.negbin <- function(model, variable) {
    if(variable == "fitted") c(NA, NA) else {
        if(is.na(match(variable, attr(model$terms, "term.labels"))))
            stop(paste(variable, "is not a term in the mean function")) else {
                newmod <- paste(" ~ . + I(", variable, "^2)")
                m2 <- update(model, newmod, start=NULL)
                c(Test= test<-m2$twologlik - model$twologlik, Pvalue=1-pchisq(test, 1))
            }}}

tukeyNonaddTest <- function(model){
    tol <- model$qr$tol
    qr <- model$qr
    fitsq <- predict(model, type="response")^2
    fitsq <- qr.resid(qr, fitsq/sqrt(sum(fitsq^2)))
    if(sd(fitsq) < tol) {
        return(c(Test=NA, Pvalue=NA))
    } else {
        r <- residuals(model, type="pearson")
        m1 <- lm(r ~ fitsq, weights=weights(model))
        df.correction <- sqrt((df.residual(model) - 1)/df.residual(m1))
        tukey <- summary(m1)$coef[2, 3] * df.correction
        c(Test=tukey, Pvalue=2*pnorm(-abs(tukey)))
    }
}



residualPlot.lm <- function(model, ...) {
    residualPlot.default(model, ...)
}

residualPlot.glm <- function(model, variable = "fitted", type = "pearson", 
                             plot = TRUE, quadratic = FALSE, 
                             smooth=TRUE, ...){
    residualPlot.default(model, variable=variable, type=type, plot=plot, 
                         quadratic=quadratic, smooth=smooth, ...)
}

