## last modified 26 March 03 by J. Fox

effect <- function (term, mod, xlevels=list(), default.levels=10, se=TRUE, 
    confidence.level=.95, transformation=family(mod)$linkinv, typical=mean){
    expand.model.frame <- function (model, extras, envir = environment(formula(model)),
    na.expand = FALSE){  # modified version of R base function
        f <- formula(model)
        data <- eval(model$call$data, envir)
        ff <- foo ~ bar + baz
        if (is.call(extras)) 
            gg <- extras
        else gg <- parse(text = paste("~", paste(extras, collapse = "+")))[[1]]
        ff[[2]] <- f[[2]]
        ff[[3]][[2]] <- f[[3]]
        ff[[3]][[3]] <- gg[[2]]
        if (!na.expand) {
            naa <- model$call$na.action
            subset <- model$call$subset
            rval <- if (is.null(data)) eval(call("model.frame", ff, # modified
                subset = subset, na.action = naa), envir)           #  lines
            else eval(call("model.frame", ff, data = data,          #
                subset = subset, na.action = naa), envir)           #
            }
        else {
            subset <- model$call$subset
            rval <- eval(call("model.frame", ff, data = data, subset = subset, 
                na.action = I), envir)
            oldmf <- model.frame(model)
            keep <- match(rownames(oldmf), rownames(rval))
            rval <- rval[keep, ]
            class(rval) <- "data.frame"
            }
        return(rval)
        }
    subscripts <- function(index){
        subs <- function(dims, index){
            dim <- length(dims)
            if (dim == 0) return(NULL)
            cum <- c(1,cumprod(dims))[dim]
            i <- index %/% cum
            if (index %% cum != 0) i <- i + 1
            c(i, subs(dims[-dim], index - (i - 1)*cum))
            }
        rev(subs(dims, index))
        }
    matrix.to.df <- function(matrix){
        on.exit(options(warn = opt[[1]]))
        opt <- options(warn = -1)
        ncol <- ncol(matrix)
        colnames <- colnames(matrix)
        result <- list()
        for (j in 1:ncol){
            numbers <- as.numeric(matrix[,j])
            result[[colnames[j]]] <-
                if(all(is.na(numbers))) matrix[,j] else numbers
            }
        as.data.frame(result)
        }
    is.relative <- function(term1, term2, factors) {
        all(!(factors[,term1]&(!factors[,term2])))
        }
    ancestors <- function(term, mod){
        names <- term.names(mod)
        if (has.intercept(mod)) names <- names[-1]
        if(length(names)==1) return(NULL)
        which.term<-which(term==names)
        if (length(which.term) == 0){
            factors <- attr(terms(mod.aug), "factors")
            result<-(1:length(names))[sapply(names,
                function(term2) is.relative(term2, term, factors))]
            if (0 ==  length(result)) which.term else result
            }
        else {
            factors <- attr(mod$terms, "factors")        
            result<-(1:length(names))[-which.term][sapply(names[-which.term],
                function(term2) is.relative(term2, term, factors))]
            if (0 ==  length(result)) which.term else result
            }
        }
    first.order.ancestors <- function(term, mod){
        ancestors <- ancestors(term, mod)
        ancestors[attr(mod$terms, 'order')[ancestors]==1]
        }
    descendants<-function(term, mod){
        names <- term.names(mod)
        if (has.intercept(mod)) names <- names[-1]
        if(length(names)==1) return(NULL)
        which.term<-which(term==names)
        if (length(which.term) == 0){
            factors <- attr(terms(mod.aug), "factors")
            (1:length(names))[sapply(names,
                function(term2) is.relative(term, term2, factors))]
            }
        else {
            factors <- attr(mod$terms, "factors")
            (1:length(names))[-which.term][sapply(names[-which.term],
                function(term2) is.relative(term, term2, factors))]
            }
        }
    is.high.order.term <- function(term, mod){
        0 == length(descendants(term, mod))
        }
    strangers <- function(term, mod){
        names <- term.names(mod)
        if (has.intercept(mod)) names <- names[-1]
        self <- which(names==term)
        ancestors <- ancestors(term, mod)
        descendants <- descendants(term, mod)
        sort(setdiff(1:ncol(attr(mod$terms, "factors")),
            union(union(ancestors, descendants), self)))
        }
    term <- gsub("\\*", ":", term)
    intercept <- has.intercept(mod)
    terms <- term.names(mod)
    if (intercept) terms <- terms[-1]
    which.term <- which(term==terms)
    if (length(which.term) == 0){
        warning(paste(term,"does not appear in the model"))
        mod.aug <- update(formula(mod), eval(parse(text=paste(". ~ . +", term))))
        }
    if (!is.high.order.term(term, mod))
        warning(paste(term, 'is not a high-order term in the model'))
    basic.vars <- first.order.ancestors(term, mod)
    all.vars <- (1:nrow(attr(mod$terms, 'factors')))[
            0 != apply(attr(mod$terms, 'factors'), 1, sum) ]
    if (intercept) all.vars <- all.vars - 1
    excluded.vars <- setdiff(all.vars, basic.vars)
    all.vars <- all.vars(as.formula(paste ("~", paste(terms[all.vars], collapse="+"))))
    basic.vars <- all.vars(as.formula(paste ("~", paste(terms[basic.vars], collapse="+"))))
    excluded.vars <- if (length(excluded.vars) > 0) 
        all.vars(as.formula(paste ("~", paste(terms[excluded.vars], collapse="+"))))
        else NULL
    X.mod <- model.matrix(mod)
    cnames <- colnames(X.mod)
    factor.cols <- rep(FALSE, length(cnames))
    names(factor.cols) <- cnames
    X <- model.frame(mod)
    for (name in all.vars){
        if (is.factor(X[[name]])) factor.cols[grep(paste("^", name, sep=""), cnames)] <- TRUE
        }
    factor.cols[grep(":", cnames)] <- FALSE   
    X <- na.omit(expand.model.frame(mod, all.vars))
    x<-list()
    factor.levels <- list()
    for (name in basic.vars){
        levels <- mod$xlevels[[name]]
        fac <- !is.null(levels)
        if (!fac) {
            levels <- if (is.null(xlevels[[name]]))
                    seq(min(X[, name]), max(X[,name]), length=default.levels)
                else if (length(xlevels[[name]]) == 1) 
                    seq(min(X[, name]), max(X[,name]), length=xlevels[[name]])
                else xlevels[[name]]
                }
            else factor.levels[[name]] <- levels
        x[[name]] <- list(name=name, is.factor=fac, levels=levels)
        }
    x.excluded <- list()
    for (name in excluded.vars){
        levels <- mod$xlevels[[name]]
        fac <- !is.null(levels)
        level <- if (fac) levels[1] else typical(X[, name])
        if (fac) factor.levels[[name]] <- levels
        x.excluded[[name]] <- list(name=name, is.factor=fac,
                                    level=level)
        }
    dims <- sapply(x, function(x) length(x$levels))
    len <- prod(dims)
    n.basic <- length(basic.vars)
    n.excluded <- length(excluded.vars)
    n.vars <- n.basic + n.excluded
    predict.data <-matrix('', len, n.vars)
    excluded <- sapply(x.excluded, function(x) x$level)
    for (i in 1:len){
        subs <- subscripts(i)
        for (j in 1:n.basic){
            predict.data[i,j] <- x[[j]]$levels[subs[j]]
            }
        if (n.excluded > 0)
            predict.data[i, (n.basic+1):n.vars] <- excluded
        }
    colnames(predict.data) <- c(sapply(x, function(x) x$name),
                                sapply(x.excluded, function(x) x$name))
    predict.data <- matrix.to.df(predict.data)
    formula.rhs <- formula(mod)[c(1,3)]   
    nrow.X <- nrow(X)
    mf <- model.frame(formula.rhs, data=rbind(X[,names(predict.data)], predict.data), 
        xlev=factor.levels)
    mod.matrix.all <- model.matrix(formula.rhs, data=mf, contrasts.arg=mod$contrasts)
    mod.matrix <- mod.matrix.all[-(1:nrow.X),]
    fit.1 <- predict(mod)
    wts <- mod$weights
    if (is.null(wts)) wts <- rep(1, length(fit.1))
    mod.2 <- lm.wfit(mod.matrix.all[1:nrow.X,], fit.1, wts)
    discrepancy <- 100*sqrt(mean(mod.2$residuals^2)/mean(mod$residuals^2))
    if (discrepancy > 1e-3) warning(paste("There is a discrepancy of", round(discrepancy, 3),
        "percent \n     in the 'safe' predictions used to generate effect", term))
    attr(mod.matrix, "assign") <- attr(mod.matrix.all, "assign")
    stranger.cols <- factor.cols & 
        apply(outer(strangers(term, mod), attr(mod.matrix,'assign'), '=='), 2, any)
    if (has.intercept(mod)) stranger.cols[1] <- TRUE
    if (any(stranger.cols)) mod.matrix[,stranger.cols] <- 
        matrix(apply(as.matrix(X.mod[,stranger.cols]), 2, mean), 
            nrow=nrow(mod.matrix), ncol=sum(stranger.cols),byrow=TRUE)
    for (name in cnames){
        components <- unlist(strsplit(name, ':'))
        if (length(components) > 1) 
            mod.matrix[,name] <- apply(mod.matrix[,components], 1, prod)
        }
    effect <- mod.matrix %*% mod.2$coefficients
    result <- list(term=term, formula=formula(mod), response=response.name(mod),
        variables=x, effect=effect, fit=transformation(effect), 
        x=predict.data[,1:n.basic, drop=FALSE], model.matrix=mod.matrix, 
        data=X, discrepancy=discrepancy)
    if (se){
        dispersion <- if (any(family(mod)$family == c('binomial', 'poisson'))) 1
            else sum(wts * mod$residuals^2)/mod$df.residual
        mod.2$terms <- mod$terms
        V <- dispersion * summary.lm(mod.2)$cov
        var <- diag(mod.matrix %*% V %*% t(mod.matrix))
        result$se <- sqrt(var)
        z <- qnorm(1 - (1 - confidence.level)/2)
        result$lower <- transformation(effect - z*result$se)
        result$upper <- transformation(effect + z*result$se)
        result$confidence.level <- confidence.level
        }
    class(result)<-'effect'
    result
    }

summary.effect <- function(object, ...){
    cat(paste("\n", gsub(":", "*", object$term), 'effect\n'))
    table <- array(object$fit,     
        dim=sapply(object$variables, function(x) length(x$levels)),
        dimnames=lapply(object$variables, function(x) x$levels))
    print(table)
    if (!is.null(object$se)){
        cat(paste('\n Lower', 100*object$confidence.level, 
            'Percent Confidence Limits\n'))
        table <- array(object$lower,   
            dim=sapply(object$variables, function(x) length(x$levels)),
            dimnames=lapply(object$variables, function(x) x$levels))
        print(table)
        cat(paste('\n Upper', 100*object$confidence.level,
            'Percent Confidence Limits\n'))
        table <- array(object$upper,   
            dim=sapply(object$variables, function(x) length(x$levels)),
            dimnames=lapply(object$variables, function(x) x$levels))
        print(table)
        }
    if (object$discrepancy > 1e-3) cat(paste("\nWarning: There is an average discrepancy of", 
        round(object$discrepancy, 3),
        "percent \n     in the 'safe' predictions for effect", object$term, '\n'))
    invisible(NULL)
    }

print.effect <- function(x, ...){
    cat(paste("\n", gsub(":", "*", x$term), 'effect\n'))
    table <- array(x$fit,     
        dim=sapply(x$variables, function(x) length(x$levels)),
        dimnames=lapply(x$variables, function(x) x$levels))
    print(table)
    if (x$discrepancy > 1e-3) cat(paste("\nWarning: There is an average discrepancy of", round(x$discrepancy, 3),
        "percent \n     in the 'safe' predictions for effect", x$term, '\n'))
    invisible(x)
    }
    
all.effects <- function(mod, ...){
    descendants<-function(term, mod){
        names <- term.names(mod)
        if (has.intercept(mod)) names <- names[-1]
        factors <- attr(mod$terms, "factors")
        if(length(names)==1) return(NULL)
        which.term<-which(term==names)
        (1:length(names))[-which.term][sapply(names[-which.term],
            function(term2) is.relative(term, term2, factors))]
        }
    is.relative <- function(term1, term2, factors) {
        all(!(factors[,term1]&(!factors[,term2])))
        }
    high.order.terms <- function(mod){
        names <- term.names(mod)
        if (has.intercept(mod)) names<-names[-1]
        rel <- lapply(names, descendants, mod=mod)
        (1:length(names))[sapply(rel, function(x) length(x)==0)]
        }
    names <- term.names(mod)
    if (has.intercept(mod)) names <- names[-1]
    terms <- names[high.order.terms(mod)]
    result <- lapply(terms, effect, mod=mod, ...)
    names(result) <- terms
    class(result) <- 'effect.list'
    result
    }
    
print.effect.list <- function(x, ...){
    cat(" model: ")
    print(x[[1]]$formula)
    for (effect in names(x)){
        print(x[[effect]])
        }
    invisible(x) 
    }

summary.effect.list <- function(object, ...){
    cat(" model: ")
    print(object[[1]]$formula)
    for (effect in names(object)){
        summary(object[[effect]])
        }
    invisible(NULL) 
    }
        
as.data.frame.effect <- function(x, row.names=NULL, optional=TRUE){
    if (is.null(x$se)) data.frame(x$x, effect=x$effect, fit=x$fit)
    else data.frame(x$x, effect=x$effect, fit=x$fit, se=x$se, lower=x$lower, upper=x$upper)
    }

plot.effect <- function(x, x.var=which.max(levels), 
    z.var=which.min(levels), multiline=is.null(x$se), rug=TRUE, xlab,
    ylab=x$response, colors=palette(), symbols=1:10, lines=1:10, cex=1.5, ylim,
    factor.names=TRUE, transform=list(link=I, inverse=I, at=NULL, n=5), ...){
    lrug <- function(x) {
                if (length(unique(x)) < 0.8 * length(x)) x <- jitter(x)
                grid.segments(x, unit(0, "npc"), x, unit(0.5, "lines"),
                    default.units="native")
                }
    ticks <- function(range, link, inverse, at, n) {
                        if (is.null(inverse)) inverse <- I
                        if (is.null(link)) link <- function(x) nlm(function(y) (inverse(y) - x)^2, 
                            mean(range))$estimate
                        if (is.null(n)) n <- 5
                        labels <- if (is.null(at)){
                            labels <- pretty(sapply(range, inverse), n=n+1)
                            }
                            else at
                        ticks <- sapply(labels, link)
                        list(at=ticks, labels=as.character(labels))
                        }
    require(lattice)
    ylab # force evaluation
    x.data <- x$data
    effect <- paste(sapply(x$variables, "[[", "name"), collapse="*")
    vars <- x$variables
    x <- as.data.frame(x)
    for (i in 1:length(vars)){
        if (!(vars[[i]]$is.factor)) next
        x[,i] <- factor(x[,i], levels=vars[[i]]$levels)
        }
    has.se <- !is.null(x$se)
    n.predictors <- ncol(x) - 2 - 3*has.se
    if (n.predictors == 1){
        range <- if (has.se) range(c(x$lower, x$upper)) else range(x$fit)
        ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
                                                range[2] + .025*(range[2] - range[1]))
        ticks <- ticks(ylim, link=transform$link, inverse=transform$inverse, 
            at=transform$at, n=transform$n)
        if (is.factor(x[,1])){
            levs <- levels(x[,1])
            print(xyplot(eval(parse(
                text=paste("fit ~ as.numeric(", names(x)[1], ")"))), 
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                panel=function(x, y, lower, upper, has.se, ...){
                    llines(x, y, lwd=2, col=colors[1], type='b', pch=19, cex=cex, ...)
                    if (has.se){
                        llines(x, lower, lty=2, col=colors[2])
                        llines(x, upper, lty=2, col=colors[2])
                        }
                    },
                ylim=ylim,
                ylab=ylab,
                xlab=if (missing(xlab)) names(x)[1] else xlab,
                scales=list(x=list(at=1:length(levs), labels=levs), 
                    y=list(at=ticks$at, labels=ticks$labels)),
                main=paste(effect, "effect plot"),
                lower=x$lower, upper=x$upper, has.se=has.se, data=x, ...))
            }        
        else {
            x.vals <- x.data[, names(x)[1]]
            print(xyplot(eval(parse(
                text=paste("fit ~", names(x)[1]))),
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                panel=function(x, y, x.vals, rug, lower, upper, has.se, ...){
                    llines(x, y, lwd=2, col=colors[1], ...)
                    if (rug) lrug(x.vals)
                    if (has.se){
                        llines(x, lower, lty=2, col=colors[2])
                        llines(x, upper, lty=2, col=colors[2])
                        }
                    },
                ylim=ylim,
                ylab=ylab,
                x.vals=x.vals, rug=rug,
                main=paste(effect, "effect plot"),
                lower=x$lower, upper=x$upper, has.se=has.se, data=x, 
                scales=list(y=list(at=ticks$at, labels=ticks$labels)), ...))
            }
        return(invisible())
        }
    predictors <- names(x)[1:n.predictors]
    levels <- sapply(apply(x[,predictors], 2, unique), length)
    if (x.var == z.var) z.var <- z.var + 1
    range <- if (has.se && (!multiline)) range(c(x$lower, x$upper)) else range(x$fit)
    ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
                                                range[2] + .025*(range[2] - range[1]))
    ticks <- ticks(ylim, link=transform$link, inverse=transform$inverse, 
        at=transform$at, n=transform$n)
    if (multiline){
        zvals <- unique(x[, z.var])
        if (length(zvals) > min(c(length(colors), length(lines), length(symbols))))
            stop(paste('Not enough colors, lines, or symbols to plot', length(zvals), 'lines'))
        if (is.factor(x[,x.var])){
            levs <- levels(x[,x.var])
            print(xyplot(eval(parse( 
                text=paste("fit ~ as.numeric(", predictors[x.var], ")",
                    if (n.predictors > 2) paste(" |", 
                    paste(predictors[-c(x.var, z.var)])), collapse="*"))),
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                panel=function(x, y, subscripts, z, ...){
                    for (i in 1:length(zvals)){
                        sub <- z[subscripts] == zvals[i]
                        llines(x[sub], y[sub], lwd=2, type='b', col=colors[i], 
                            pch=symbols[i], lty=lines[i], cex=cex, ...)
                        }
                    },
                ylim=ylim,
                ylab=ylab,
                xlab=if (missing(xlab)) predictors[x.var] else xlab,
                z=x[,z.var],
                scales=list(x=list(at=1:length(levs), labels=levs), 
                    y=list(at=ticks$at, labels=ticks$labels)),
                zvals=zvals,
                main=paste(effect, "effect plot"),
                key=list(title=predictors[z.var], cex.title=1, border=TRUE,
                    text=list(as.character(zvals)), 
                    lines=list(col=colors[1:length(zvals)], lty=lines[1:length(zvals)], lwd=2), 
                    points=list(pch=1:length(zvals))),
                data=x, ...))
            }    
        else{
        x.vals <- x.data[, names(x)[x.var]]
        print(xyplot(eval(parse( 
                text=paste("fit ~", predictors[x.var], 
                    if (n.predictors > 2) paste(" |", 
                    paste(predictors[-c(x.var, z.var)])), collapse="*"))),
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                panel=function(x, y, subscripts, x.vals, rug, z, ...){
                    if (rug) lrug(x.vals)
                    for (i in 1:length(zvals)){
                        sub <- z[subscripts] == zvals[i]
                        llines(x[sub], y[sub], lwd=2, type='l', col=colors[i], lty=lines[i], cex=cex, ...)
                        }
                    },
                ylim=ylim,
                ylab=ylab,
                xlab=if (missing(xlab)) predictors[x.var] else xlab,
                x.vals=x.vals, rug=rug,
                z=x[,z.var],
                zvals=zvals,
                main=paste(effect, "effect plot"),
                key=list(title=predictors[z.var], cex.title=1, border=TRUE,
                    text=list(as.character(zvals)), 
                    lines=list(col=colors[1:length(zvals)], lty=lines[1:length(zvals)], lwd=2)), 
                data=x, scales=list(y=list(at=ticks$at, labels=ticks$labels)), ...))
            }
        return(invisible())
        }
    if (is.factor(x[,x.var])){
        levs <- levels(x[,x.var])
        print(xyplot(eval(parse( 
            text=paste("fit ~ as.numeric(", predictors[x.var], ") |", 
                paste(predictors[-x.var], collapse="*")))),
            strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
            panel=function(x, y, subscripts, lower, upper, has.se, ...){
                llines(x, y, lwd=2, type='b', col=colors[1], pch=19, cex=cex, ...)
                if (has.se){
                    llines(x, lower[subscripts], lty=2, col=colors[2])
                    llines(x, upper[subscripts], lty=2, col=colors[2])
                    }
                },
            ylim=ylim,
            ylab=ylab,
            xlab=if (missing(xlab)) predictors[x.var] else xlab,
            scales=list(x=list(at=1:length(levs), labels=levs), 
                y=list(at=ticks$at, labels=ticks$labels)),
            main=paste(effect, "effect plot"),
            lower=x$lower, upper=x$upper, has.se=has.se, data=x, ...))
        }    
    else{
        x.vals <- x.data[, names(x)[x.var]]
        print(xyplot(eval(parse( 
            text=paste("fit ~", predictors[x.var], "|", 
                paste(predictors[-x.var], collapse="*")))),
            strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
            panel=function(x, y, subscripts, x.vals, rug, lower, upper, has.se, ...){
                llines(x, y, lwd=2, col=colors[1], ...)
                if (rug) lrug(x.vals)
                if (has.se){
                    llines(x, lower[subscripts], lty=2, col=colors[2])
                    llines(x, upper[subscripts], lty=2, col=colors[2])
                    }
                },
            ylim=ylim,
            ylab=ylab,
            xlab=if (missing(xlab)) predictors[x.var] else xlab,
            x.vals=x.vals, rug=rug,
            main=paste(effect, "effect plot"),
            lower=x$lower, upper=x$upper, has.se=has.se, data=x, 
            scales=list(y=list(at=ticks$at, labels=ticks$labels)), ...))
        }
    }

plot.effect.list <- function(x, selection, ...){
    if (!missing(selection)){
        plot(x[[selection]], ...)
        return(invisible())
        }
    effects <- gsub(":", "*", names(x))
    repeat {
        selection <- menu(effects)
        if (selection == 0) break
        else plot(x[[selection]], ...)
        }
    }
