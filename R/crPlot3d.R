# added 2022-05-24 by J. Fox
# 2023-01-01: use rgl::*3d() function. Duncan Murdoch
# 2024-04-11: invisibly return coordinates. J. Fox

crPlot3d <- function(model, var1, var2, ...) {
  UseMethod("crPlot3d")
}

crPlot3d.lm <- function (model,
                         var1,
                         var2,
                         xlab = var1,
                         ylab = paste0("C+R(", eff$response, ")"),
                         zlab = var2,
                         axis.scales = TRUE,
                         axis.ticks = FALSE,
                         revolutions = 0,
                         bg.col = c("white", "black"),
                         axis.col = if (bg.col == "white") c("darkmagenta", "black", "darkcyan") 
                           else c("darkmagenta", "white", "darkcyan"),
                         surface.col = carPalette()[2:3],
                         surface.alpha = 0.5,
                         point.col = "yellow",
                         text.col = axis.col,
                         grid.col = if (bg.col ==  "white") "black" else "gray",
                         fogtype = c("exp2", "linear", "exp", "none"),
                         fill = TRUE,
                         grid = TRUE,
                         grid.lines = 26,
                         smoother = c("loess", "mgcv", "none"),
                         df.mgcv = NULL,
                         loess.args = NULL,
                         sphere.size = 1,
                         radius = 1,
                         threshold = 0.01,
                         speed = 1,
                         fov = 60,
                         ellipsoid = FALSE,
                         level = 0.5,
                         ellipsoid.alpha = 0.1,
                         id = FALSE,
                         mouseMode=c(none="none", left="polar", right="zoom", middle="fov", 
                                     wheel="pull"),
                         ...)
{
  smoother <- match.arg(smoother)
  
  if (!requireNamespace("rgl")) stop("rgl package missing")
  if (!requireNamespace("mgcv") && smoother == "mgcv") stop("mgcv package missing")
  if (!requireNamespace("effects")) stop("effects package missing")
  
  rgl::par3d(mouseMode=mouseMode)
  
  loess.args <- applyDefaults(
    loess.args,
    defaults = list(
      span = 2/3,
      family = if (inherits(model, "glm") && family(model)$family %in% 
                   c("binomial", "poisson", "quasibinomial", "quasipoisson")) 
        "gaussian" else "symmetric",
      degree = 1
    )
  )
  
  levels <- list(grid.lines, grid.lines)
  names(levels) <- c(var1, var2)
  eff <- effects::Effect(c(var1, var2),
                         model,
                         residuals = TRUE,
                         xlevels = levels)
  
  x <- eff$data[[var1]]
  z <- eff$data[[var2]]
  n <- length(x)
  
  id <- applyDefaults(
    id,
    defaults = list(
      method = "mahal",
      n = 2,
      labels = as.character(seq(along = x)),
      offset = ((100 / length(x)) ^ (1 / 3)) *
        0.02
    ),
    type = "id"
  )
  if (isFALSE(id)) {
    id.n <- 0
    id.method <- "none"
    labels <- NULL
  }
  else {
    labels <- rownames(eff$data)
    id.method <- id$method
    id.n <- if ("identify" == id.method) Inf else id$n
    offset <- id$offset
  }
  
  bg.col <- match.arg(bg.col)
  fogtype <- match.arg(fogtype)
  
  rgl::next3d()
  rgl::view3d(fov = fov)
  rgl::bg3d(color = bg.col, fogtype = fogtype)
  
  x.grid <- eff$variables[[var1]]$levels
  z.grid <- eff$variables[[var2]]$levels
  xz.grid <- expand.grid(x.grid, z.grid)
  
  iqr.x <- diff(quantile(x, c(.75, .25)))
  iqr.z <- diff(quantile(z, c(.75, .25)))
  
  which.grid <- rep(0, n)
  for (i in 1:n){
    which.grid[i] <- which.min(((x[i] - xz.grid[, 1])/iqr.x)^2 + 
                                 ((z[i] - xz.grid[, 2])/iqr.z)^2)
  }
  
  x <- xz.grid[which.grid, 1]
  z <- xz.grid[which.grid, 2]
  
  yhat <- as.vector(eff$fit)
  y <- yhat[which.grid] + eff$residuals
  
  minx <- min(x)
  maxx <- max(x)
  miny <- min(y)
  maxy <- max(y)
  minz <- min(z)
  maxz <- max(z)
  if (axis.scales) {
    lab.min.x <- nice(minx)
    lab.max.x <- nice(maxx)
    lab.min.y <- nice(miny)
    lab.max.y <- nice(maxy)
    lab.min.z <- nice(minz)
    lab.max.z <- nice(maxz)
    minx <- min(lab.min.x, minx)
    maxx <- max(lab.max.x, maxx)
    miny <- min(lab.min.y, miny)
    maxy <- max(lab.max.y, maxy)
    minz <- min(lab.min.z, minz)
    maxz <- max(lab.max.z, maxz)
    min.x <- (lab.min.x - minx) / (maxx - minx)
    max.x <- (lab.max.x - minx) / (maxx - minx)
    min.y <- (lab.min.y - miny) / (maxy - miny)
    max.y <- (lab.max.y - miny) / (maxy - miny)
    min.z <- (lab.min.z - minz) / (maxz - minz)
    max.z <- (lab.max.z - minz) / (maxz - minz)
    if (axis.ticks) {
      if (axis.scales) {
        x.labels <- seq(lab.min.x, lab.max.x, by = diff(range(lab.min.x, lab.max.x))/4)
        x.at <- seq(min.x, max.x, by = nice(diff(range(min.x, max.x))/4))
        rgl::text3d(x.at,-0.05, 0, x.labels, col = axis.col[1])
        z.labels <-seq(lab.min.z, lab.max.z, by = diff(range(lab.min.z, lab.max.z))/4)
        z.at <- seq(min.z, max.z, by = diff(range(min.z, max.z))/4)
        rgl::text3d(0,-0.1, z.at, z.labels, col = axis.col[3])
        y.labels <- seq(lab.min.y, lab.max.y, by = diff(range(lab.min.y, lab.max.y))/4)
        y.at <- seq(min.y, max.y, by = diff(range(min.y, max.y))/4)
        rgl::text3d(-0.05, y.at,-0.05, y.labels, col = axis.col[2])
      }
    }
    else {
      rgl::text3d(min.x,-0.05, 0, lab.min.x, col = axis.col[1])
      rgl::text3d(max.x,-0.05, 0, lab.max.x, col = axis.col[1])
      rgl::text3d(0,-0.1, min.z, lab.min.z, col = axis.col[3])
      rgl::text3d(0,-0.1, max.z, lab.max.z, col = axis.col[3])
      rgl::text3d(-0.05, min.y,-0.05, lab.min.y, col = axis.col[2])
      rgl::text3d(-0.05, max.y,-0.05, lab.max.y, col = axis.col[2])
    }
  }
  
  x <- (x - minx) / (maxx - minx)
  y <- (y - miny) / (maxy - miny)
  z <- (z - minz) / (maxz - minz)
  
  size <- sphere.size * ((100 / length(x)) ^ (1 / 3)) * 0.015
  radius <- radius / median(radius)
  if (size > threshold)
    rgl::spheres3d(x, y, z, color = point.col, radius = size*radius)
  else
    rgl::points3d(x, y, z, color = point.col)
  
  if (!axis.scales) axis.col[1] <- axis.col[3] <- axis.col[2]
  rgl::segments3d(c(0, 1), c(0, 0), c(0, 0), color = axis.col[1])
  rgl::segments3d(c(0, 0), c(0, 1), c(0, 0), color = axis.col[2])
  rgl::segments3d(c(0, 0), c(0, 0), c(0, 1), color = axis.col[3])
  rgl::text3d(1, 0, 0, xlab, adj = 1, color = axis.col[1])
  rgl::text3d(0, 1.05, 0, ylab, adj = 1, color = axis.col[2])
  rgl::text3d(0, 0, 1, zlab, adj = 1, color = axis.col[3])
  if (ellipsoid) {
    dfn <- 3
    dfd <- length(x) - 1
    ell.radius <- sqrt(dfn * qf(level, dfn, dfd))
    ellips <- ellipsoid(center = c(mean(x), mean(y), mean(z)), 
                        shape = cov(cbind(x, y, z)),
                        radius = ell.radius
    )
    if (fill) rgl::shade3d(ellips, col = surface.col[1], alpha = ellipsoid.alpha, lit = FALSE)
    if (grid) rgl::wire3d(ellips, col = surface.col[1], lit = FALSE)
  } 
  
  x.grid <- (x.grid - minx)/(maxx - minx)
  z.grid <- (z.grid - minz)/(maxz - minz)
  yhat <- (yhat - miny) / (maxy - miny)
  yhat <- matrix(yhat, grid.lines, grid.lines)
  
  rgl::surface3d(x = x.grid,
                 z = z.grid,
                 y = yhat,
                 color = surface.col[1],
                 alpha = surface.alpha,
                 lit = FALSE)
  if (grid)
    rgl::surface3d(
      x = x.grid,
      z = z.grid,
      y = yhat,
      color = if (fill) grid.col else surface.col[1],
      alpha = surface.alpha,
      lit = FALSE,
      front = "lines",
      back = "lines"
    )
  if (smoother != "none"){
    smooth <- if (smoother == "loess") {
      loess(y ~ x*z, span=loess.args$span, degree=loess.args$degree, family=loess.args$family)
    } else {
      if (is.null(df.mgcv)) mgcv::gam(y ~ s(x, z))
      else mgcv::gam(y ~ s(x, z, fx=TRUE, k=df.mgcv))
    }
    
    yhat <- predict(smooth, newdata=data.frame(expand.grid(x=x.grid, z=z.grid)))
    yhat <- matrix(yhat, grid.lines, grid.lines)
    
    rgl::surface3d(x = x.grid,
                   z = z.grid,
                   y = yhat,
                   color = surface.col[2],
                   alpha = surface.alpha,
                   lit = FALSE)
    if (grid)
      rgl::surface3d(
        x = x.grid,
        z = z.grid,
        y = yhat,
        color = if (fill) grid.col else surface.col[2],
        alpha = surface.alpha,
        lit = FALSE,
        front = "lines",
        back = "lines"
      )
  }
  
  if (id.method == "identify"){
    Identify3d(x, y, z, axis.scales=axis.scales, labels=labels, 
               col=surface.col[1], offset=offset)
  }
  else if (id.method != "none") {
    showLabels3d(
      x,
      y,
      z,
      labels,
      id.method = id.method,
      id.n = id.n,
      col = surface.col[1],
      offset = offset
    )
  }
  
  if (revolutions > 0) {
    for (i in 1:revolutions) {
      for (angle in seq(1, 360, length.out = 360 / speed))
        rgl::view3d(-angle, fov = fov)
    }
  }
  
  D <- data.frame(
    x = minx + (maxx - minx)*x,
    y = miny + (maxy - miny)*y,
    z = minz + (maxz - minz)*z
  )
  return(invisible(D))
}
