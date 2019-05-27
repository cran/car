# added 2017-12-13 by J. Fox
# 2017-12-14: improved recovery of model data
#             removed faulty one-step approximations
# 2018-01-28: fix computation of Cook's D for lme models
# 2018-05-23: fixed bug when more than one grouping variable (reported by Maarten Jung)
# 2018-06-07: skip plot of "sigma^2" in GLMM if dispersion fixed to 1; improved labelling for covariance components
# 2018-11-04: tweak to dfbetas.influence.merMod() suggested by Ben Bolker.
# 2018-11-09: parallel version of influence.merMod()

# influence diagnostics for mixed models

globalVariables(".groups")

influence.merMod <- function(model, groups, data, maxfun=1000,
                             ncores=1, ...) {
    .vcov <- function(x) Matrix::as.matrix(vcov(x))
    if (is.infinite(ncores)) {
        ncores <- parallel::detectCores(logical=FALSE)
    }
    if (missing(data)){
        data <- getCall(model)$data
        data <- if (!is.null(data)) eval(data, parent.frame())
        else stop("model did not use the data argument")
    }
    if (missing(groups)) {
        groups <- ".case"
        data$.case <- rownames(data)
    }
    else if (length(groups) > 1){
        del.var <- paste0(groups, collapse=".")
        data[, del.var] <- apply(data, 1, paste0, collapse=".")
        groups <- del.var
    }
    unique.del <- unique(data[, groups])
    data[[".groups"]] <- data[, groups]
    par <- list(theta=lme4::getME(model, "theta"))
    if (inherits(model, "glmerMod")) par$fixef <- fixef(model)
    fixed <- fixef(model)
    Vs <- lme4::VarCorr(model)
    nms <- names(Vs)
    sep <- ":"
    if (length(nms) == 1) {
        nms <- ""
        sep <- ""
    }
    vc <- lme4::getME(model, "sigma")^2
    names(vc) <- "sigma^2"
    for (i in 1:length(Vs)){
        V <- Vs[[i]]
        c.names <- colnames(V)
        e.names <- outer(c.names, c.names, function(a, b) paste0("C[", a, ",", b, "]"))
        diag(e.names) <- paste0("v[", c.names, "]")
        v <- V[lower.tri(V, diag=TRUE)]
        names(v) <- paste0(nms[i], sep, e.names[lower.tri(e.names, diag=TRUE)])
        vc <- c(vc, v)
    }
    control <- if (inherits(model, "lmerMod")) lme4:: lmerControl(optCtrl=list(maxfun=maxfun))
    else if (inherits(model, "glmerMod")) lme4:: glmerControl(optCtrl=list(maxfun=maxfun))
    deleteGroup <- function(del){
        data$del <- del
        mod.1 <- suppressWarnings(update(model, data=data,
                                         subset=(.groups != del),
                                         start=par,
                                         control=control))
        opt <- mod.1@optinfo
        feval <- opt$feval
        converged <- opt$conv$opt == 0 && length(opt$warnings) == 0
        fixed.1 <- fixef(mod.1)
        Vs.1 <- lme4::VarCorr(mod.1)
        vc.0 <- lme4::getME(mod.1, "sigma")^2
        for (V in Vs.1){
            vc.0 <- c(vc.0, V[lower.tri(V, diag=TRUE)])
        }
        vc.1 <- vc.0
        vcov.1 <<- .vcov(mod.1)
        list(fixed.1=fixed.1, vc.1=vc.1, vcov.1=vcov.1, converged=converged, feval=feval)
    }
    result <- if(ncores >= 2){
        message("Note: using a cluster of ", ncores, " cores")
        cl <- parallel::makeCluster(ncores)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterEvalQ(cl, require("lme4"))
        parallel::clusterApply(cl, unique.del, deleteGroup)
    } else {
        lapply(unique.del, deleteGroup)
    }
    result <- combineLists(result)
    fixed.1 <- result$fixed.1
    rownames(fixed.1) <- unique.del
    colnames(fixed.1) <- names(fixed)
    vc.1 <- result$vc.1
    rownames(vc.1) <- unique.del
    colnames(vc.1) <- names(vc)
    feval <- as.vector(result$feval)
    converged <- as.vector(result$converged)
    vcov.1 <- result$vcov.1
    names(vcov.1) <- names(feval) <- names(converged) <- unique.del
    left <- "[-"
    right <- "]"
    if (groups == ".case") {
        groups <- "case"
    }
    nms <- c("fixed.effects", paste0("fixed.effects", left, groups, right),
             "var.cov.comps", paste0("var.cov.comps", left, groups, right),
             "vcov", paste0("vcov", left, groups, right),
             "groups", "deleted", "converged", "function.evals")
    result <- list(fixed, fixed.1, vc, vc.1,
                   .vcov(model), vcov.1, groups, unique.del, converged, feval)
    names(result) <- nms
    class(result) <- "influence.merMod"
    result
}


# influence.merMod <- function(model, groups, data, maxfun=1000, ...){
#     if (missing(data)){
#         data <- getCall(model)$data
#         data <- if (!is.null(data)) eval(data, parent.frame())
#         else stop("model did not use the data argument")
#     }
#     if (missing(groups)) {
#         groups <- ".case"
#         data$.case <- rownames(data)
#     }
#     else if (length(groups) > 1){
#         del.var <- paste0(groups, collapse=".")
#         data[, del.var] <- apply(data[, groups], 1, function (row) paste0(row, collapse="."))
#         groups <- del.var
#     }
#     unique.del <- unique(data[, groups])
#     data$.groups <- data[, groups]
#     par <- list(theta=lme4::getME(model, "theta"))
#     if (inherits(model, "glmerMod")) par$fixef <- fixef(model)
#     fixed <- fixef(model)
#     fixed.1 <- matrix(0, length(unique.del), length(fixed))
#     rownames(fixed.1) <- unique.del
#     colnames(fixed.1) <- names(fixed)
#     Vs <- lme4::VarCorr(model)
#     nms <- names(Vs)
#     sep <- ":"
#     if (length(nms) == 1) {
#         nms <- ""
#         sep <- ""
#     }
#     vc <- lme4::getME(model, "sigma")^2
#     names(vc) <- "sigma^2"
#     for (i in 1:length(Vs)){
#         V <- Vs[[i]]
#         c.names <- colnames(V)
#         e.names <- outer(c.names, c.names, function(a, b) paste0("C[", a, ",", b, "]"))
#         diag(e.names) <- paste0("v[", c.names, "]")
#         v <- V[lower.tri(V, diag=TRUE)]
#         names(v) <- paste0(nms[i], sep, e.names[lower.tri(e.names, diag=TRUE)])
#         vc <- c(vc, v)
#     }
#     vc.1 <- matrix(0, length(unique.del), length(vc))
#     rownames(vc.1) <- unique.del
#     colnames(vc.1) <- names(vc)
#     feval <- numeric(length(unique.del))
#     converged <- logical(length(unique.del))
#     vcov.1 <- vector(length(unique.del), mode="list")
#     names(vcov.1) <- names(feval) <- names(converged) <- unique.del
#     # control <- if (one.step){
#     #     if (inherits(model, "lmerMod")) lme4::lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B", maxit=1))
#     #     else if (inherits(model, "glmerMod")) lme4::glmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B", maxit=1))
#     # } else {
#     control <- if (inherits(model, "lmerMod")) lme4::lmerControl(optCtrl=list(maxfun=maxfun))
#         else if (inherits(model, "glmerMod")) lme4::glmerControl(optCtrl=list(maxfun=maxfun))
#     # }
#     if (is.null(control)) stop("model not an lmer or glmer object")
#     for (del in unique.del){
#         data$del <- del
#         mod.1 <- suppressWarnings(update(model, data=data, subset=.groups != del, start=par,
#                                          control=control))
# 
#         opt <- mod.1@optinfo
#         feval[del] <- opt$feval
#         converged[del] <- opt$conv$opt == 0 && length(opt$warnings) == 0
#         fixed.1[del, ] <- fixef(mod.1)
#         Vs.1 <- lme4::VarCorr(mod.1)
#         vc.0 <- lme4::getME(mod.1, "sigma")^2
#         for (V in Vs.1){
#             vc.0 <- c(vc.0, V[lower.tri(V, diag=TRUE)])
#         }
#         vc.1[del, ] <- vc.0
#         vcov.1[[del]] <- vcov(mod.1)
#     }
#     left <- "[-"
#     right <- "]"
#     if (groups == ".case") {
#         groups <- "case"
#     }
#     nms <- c("fixed.effects", paste0("fixed.effects", left, groups, right),
#              "var.cov.comps", paste0("var.cov.comps", left, groups, right),
#              "vcov", paste0("vcov", left, groups, right),
#              "groups", "deleted", "converged", "function.evals")
#     result <- list(fixed, fixed.1, vc, vc.1, vcov(model), vcov.1, groups, unique.del, converged, feval)
#     names(result) <- nms
#     class(result) <- "influence.merMod"
#     result
# }

dfbeta.influence.lme <- function(model, which=c("fixed", "var.cov"), ...){
    if (missing(which)) which <- "fixed"
    dfbeta.influence.merMod(model, which=which, ...)
}

dfbeta.influence.merMod <- function(model, which=c("fixed", "var.cov"), ...){
    which <- match.arg(which)
    b <- if (which == "fixed") model[[2]] else model[[4]]
    b0 <- if (which == "fixed") model[[1]] else model[[3]]
    b - matrix(b0, nrow=nrow(b), ncol=ncol(b), byrow=TRUE)
}

dfbetas.influence.lme <- function(model, ...){
    dfbetas.influence.merMod(model, ...)
}

dfbetas.influence.merMod <- function(model, ...){
    dfbeta(model)/t(sapply(model[[6]], function(x) sqrt(diag(as.matrix(x)))))
}

cooks.distance.influence.lme <- function(model, ...){
    db <- dfbeta(model)
    n <- nrow(db)
    p <- ncol(db)
    d <- numeric(n)
    vcovs <- model[[6]]
    sig.sq <- (exp(model[[4]][, ncol(model[[4]])]))^2
    for (i in 1:n){
        d[i] <- (db[i, ] %*% solve(vcovs[[i]]) %*% db[i, ])/(p*sig.sq[i])
    }
    d
}

cooks.distance.influence.merMod <- function(model, ...){
    db <- dfbeta(model)
    n <- nrow(db)
    p <- ncol(db)
    d <- numeric(n)
    vcovs <- model[[6]]
    sig.sq <- model[[4]][, 1]
    for (i in 1:n){
        d[i] <- (db[i, ] %*% solve(vcovs[[i]]) %*% db[i, ])/(p*sig.sq[i])
    }
    d
}

influence.lme <- function(model, groups, data, ncores=1, ...){
    if (is.infinite(ncores)) {
        ncores <- parallel::detectCores(logical=FALSE)
    }
    if (missing(data)) data <- model$data
    if (is.null(data)){
        data <- getCall(model)$data
        data <- if (!is.null(data)) eval(data, parent.frame())
        else stop("model did not use the data argument")
    }
    if (missing(groups)) {
        groups <- ".case"
        data$.case <- rownames(data)
    }
    else if (length(groups) > 1){
        del.var <- paste0(groups, collapse=".")
        data[, del.var] <- apply(data[, groups], 1, function (row) paste0(row, collapse="."))
        groups <- del.var
    }
    unique.del <- unique(data[, groups])
    data$.groups <- data[, groups]
    fixed <- fixef(model)
    fixed.1 <- matrix(0, length(unique.del), length(fixed))
    rownames(fixed.1) <- unique.del
    colnames(fixed.1) <- names(fixed)
    vc <- attr(model$apVar, "Pars")
    vc.1 <- matrix(0, length(unique.del), length(vc))
    rownames(vc.1) <- unique.del
    colnames(vc.1) <- names(vc)
    vcov.1 <- vector(length(unique.del), mode="list")
    names(vcov.1) <-  unique.del
    deleteGroup <- function(del){
        data$del <- del
        mod.1 <- suppressWarnings(update(model, data=data, subset=.groups != del,
                                         control=nlme::lmeControl(returnObject=TRUE)))
        fixed.1 <- fixef(mod.1)
        vc.0 <- attr(mod.1$apVar, "Pars")
        vc.1 <- if (!is.null(vc.0)) vc.0 else rep(as.numeric(NA), length(vc))
        vcov.1 <- vcov(mod.1)
        list(fixed.1=fixed.1, vc.1=vc.1, vcov.1=vcov.1)
    }
    result <- if(ncores >= 2){
        message("Note: using a cluster of ", ncores, " cores")
        cl <- parallel::makeCluster(ncores)
        on.exit(parallel::stopCluster(cl))
        parallel::clusterEvalQ(cl, require("nlme"))
        parallel::clusterApply(cl, unique.del, deleteGroup)
    } else {
        lapply(unique.del, deleteGroup)
    }
    result <- combineLists(result)
    left <- "[-"
    right <- "]"
    if (groups == ".case") {
        groups <- "case"
    }
    nms <- c("fixed.effects", paste0("fixed.effects", left, groups, right),
             "var.cov.comps", paste0("var.cov.comps", left, groups, right),
             "vcov", paste0("vcov", left, groups, right),
             "groups", "deleted")
    result <- list(fixed, fixed.1=result$fixed.1, vc, vc.1=result$vc.1, vcov(model), 
                   vcov.1=result$vcov.1, groups, unique.del)
    names(result) <- nms
    class(result) <- "influence.lme"
    result
}

# influence.lme <- function(model, groups, data, ...){
#     if (missing(data)) data <- model$data
#     if (is.null(data)){
#         data <- getCall(model)$data
#         data <- if (!is.null(data)) eval(data, parent.frame())
#         else stop("model did not use the data argument")
#     }
#     if (missing(groups)) {
#         groups <- ".case"
#         data$.case <- rownames(data)
#     }
#     else if (length(groups) > 1){
#         del.var <- paste0(groups, collapse=".")
#         data[, del.var] <- apply(data[, groups], 1, function (row) paste0(row, collapse="."))
#         groups <- del.var
#     }
#     unique.del <- unique(data[, groups])
#     data$.groups <- data[, groups]
#     fixed <- fixef(model)
#     fixed.1 <- matrix(0, length(unique.del), length(fixed))
#     rownames(fixed.1) <- unique.del
#     colnames(fixed.1) <- names(fixed)
#     vc <- attr(model$apVar, "Pars")
#     vc.1 <- matrix(0, length(unique.del), length(vc))
#     rownames(vc.1) <- unique.del
#     colnames(vc.1) <- names(vc)
#     vcov.1 <- vector(length(unique.del), mode="list")
#     names(vcov.1) <-  unique.del
#     for (del in unique.del){
#         data$del <- del
#         mod.1 <- suppressWarnings(update(model, data=data, subset=.groups != del,
#                                          control=nlme::lmeControl(returnObject=TRUE)))
#         fixed.1[del, ] <- fixef(mod.1)
#         vc.0 <- attr(mod.1$apVar, "Pars")
#         vc.1[del, ] <- if (!is.null(vc.0)) vc.0 else rep(NA, length(vc))
#         vcov.1[[del]] <- vcov(mod.1)
#     }
#     left <- "[-"
#     right <- "]"
#     if (groups == ".case") {
#         groups <- "case"
#     }
#     nms <- c("fixed.effects", paste0("fixed.effects", left, groups, right),
#              "var.cov.comps", paste0("var.cov.comps", left, groups, right),
#              "vcov", paste0("vcov", left, groups, right),
#              "groups", "deleted")
#     result <- list(fixed, fixed.1, vc, vc.1, vcov(model), vcov.1, groups, unique.del)
#     names(result) <- nms
#     class(result) <- "influence.lme"
#     result
# }

infIndexPlot.influence.lme <- function(model, vars=c("dfbeta", "dfbetas", "var.cov.comps", "cookd"), id=TRUE, grid=TRUE,
                                          main="Diagnostic Plots", ...){
    if (missing(vars)) vars <- c("dfbeta", "cookd")
    infIndexPlot.influence.merMod(model, vars=vars, id=id, grid=grid, main=main)
}

infIndexPlot.influence.merMod <- function(model, vars=c("dfbeta", "dfbetas", "var.cov.comps", "cookd"), id=TRUE, grid=TRUE,
         main="Diagnostic Plots", ...){
    id <- applyDefaults(id, defaults=list(method="y", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
    if (isFALSE(id)){
        id.n <- 0
        id.method <- "none"
        labels <- id.cex <- id.col <- id.location <- NULL
    }
    else{
        labels <- id$labels
        if (is.null(labels)) labels <- row.names(model[[2]])
        id.method <- id$method
        id.n <- if ("identify" %in% id.method) Inf else id$n
        id.cex <- id$cex
        id.col <- id$col
        id.location <- id$location
    }
    if (missing(vars)) vars <- c("dfbeta", "cookd")
    what <- pmatch(tolower(vars),
                   c("dfbeta", "dfbetas", "var.cov.comps", "cookd"))
    if(length(what) < 1) stop("Nothing to plot")
    X <- cbind(if (1 %in% what) dfbeta(model), if (2 %in% what) dfbetas(model),
               if (3 %in% what) dfbeta(model, "var.cov"), if (4 %in% what) cooks.distance(model))
    if (4 %in% what) colnames(X)[ncol(X)] <- "Cook's D"
    names <- colnames(X)
    # check for row.names, and use them if they are numeric.
    oldwarn <- options()$warn
    options(warn=-1)
    xaxis <- as.numeric(row.names(model[[2]]))
    options(warn=oldwarn)
    if (any (is.na(xaxis))) xaxis <- 1:length(xaxis)
    plotnum <- 0
    nplots <- ncol(X)
    if ("sigma^2" %in% names){
        if (all(X[, "sigma^2"] == 0)){ # check for fixed dispersion
            X <- X[, names != "sigma^2", drop=FALSE]
            names <- names[names != "sigma^2"]
            nplots <- nplots - 1
        }
    }
    op <- par(mfrow=c(nplots, 1), mar=c(1, 4, 0, 2) + .0,
              mgp=c(2, 1, 0), oma=c(6, 0, 6, 0))
    on.exit(par(op))
    for (j in 1:nplots){
        plotnum <- plotnum + 1
        plot(xaxis, X[, j], type="n", ylab=names[j], xlab="", xaxt="n", tck=0.1, ...)
        if(grid){
            grid(lty=1, equilogs=FALSE)
            box()}
        points(xaxis, X[, j], type="h", ...) #}
        points(xaxis, X[, j], type="p", ...)
        abline(h=0, lty=2 )
        axis(1, labels= ifelse(plotnum < nplots, FALSE, TRUE))
        showLabels(xaxis, X[, j], labels=labels,
                   method=id.method, n=id.n, cex=id.cex,
                   col=id.col, location=id.location)
    }
    mtext(side=3, outer=TRUE ,main, cex=1.2, line=1)
    mtext(side=1, outer=TRUE, paste0("Index(", model$groups, ")"), line=3)
    invisible()
}
