# last modified 25 Februrary 2010 by J. Fox
# rewritten 15 April 2010 S Weisberg
# 2013-02-07 S. Weisberg bug fix for use with 'scatterplot' with groups.
#   Added an argument to showLabels1 'all' that gives a list of two
#   elements for the original labels and subset indicator.  See
#   scatterplot.R for an example of its use.
#   If a list of cases to be labelled is supplied, id.n is needed only
#   if all n labels are to be printed.
# 2014-03-12 added new id.method "r" that labels using order(abs(y), decreasing=TRUE)
# 2016-05-16 added argument id.location = c("lr", "ab") for location of point labels
# 2017-01-08 added "avoid" to id.location arg. J. Fox
# 2017-01-08 removed ".id" from arg names for showLabels()
# 2017-01-10 special handling for method="none".
# 2017-02-13 fixed showLabels1() when location="avoid"
# 2017-03-25: don't supply names if indexes are the same as labels. J. Fox


showLabels <- function(x, y, labels=NULL, method="identify",
     n = length(x), cex=1, col=carPalette()[1], 
     location=c("lr", "ab", "avoid"), ...) {
  location <- match.arg(location)
  res <- NULL
  method <- if(is.list(method)) method else list(method)
  for (meth in method){
     if (length(meth) == 1 && is.character(meth) && meth == "none") next
     res <- c(res, showLabels1(x, y, labels, meth, n, cex,
              col, location, ...))
  }
  return(if(is.null(res)) invisible(res) else res)
  }

showLabels1 <- function(x, y, labels=NULL, id.method="identify",
	   id.n = length(x), id.cex=1, id.col=carPalette()[1], 
	   id.location="lr", all=NULL, ...) { 
# If labels are NULL, try to get the labels from x:
  if (is.null(labels)) labels <- names(x)
  if (is.null(labels)) labels <- paste(seq_along(x))
  if (is.null(id.col)) id.col <- carPalette()[1]
  if (is.null(id.location)) id.location <- "lr"
# logged-axes?
  log.x <- par("xlog")
  log.y <- par("ylog")
# id.method can be any of the following:
#    --- a list of row numbers
#    --- a list of labels
#    --- a vector of n numbers
#    --- a text string:  'identify', 'x', 'y', 'mahal', 'r'
  idmeth <- pmatch(id.method[1], c("x", "y", "mahal", "identify", "r"))
  if(!is.na(idmeth)) 
    idmeth <- c("x", "y", "mahal", "identify", "r")[idmeth]
# if idmeth is NA, then id.method must be <= n numbers or labels
  id.var <- NULL
  if(is.na(idmeth)){
    if(is.null(all)) 
      all <- list(labels=labels, subs=rep(TRUE, length(labels)))
    names(all$labels) <- all$labels
    if(length(id.method) >= length(x)){
      id.var <- id.method[which(all$subs)]
      id.n <- min(id.n, length(id.var))
      }
    else {
      id.var <- rep(0, length(x))
      names(id.var) <- labels
      inSubset <- all$labels[all$subs] %in% all$labels[id.method]
      id.var[inSubset] <- 1
      id.n <- sum(inSubset)
      }
  }
  else {
# use identify?
  if(idmeth == "identify"){
    	  result <- labels[identify(x, y, labels, n=length(x), 
    	                            cex=id.cex, col=id.col)]
    	  if(length(result) > 0) return(unique(result)) else return(NULL)
  }
# missing values need to be removed
	ismissing <- is.na(x) | is.na(y) | is.na(labels)
	if( any(ismissing) ) {
		x <- x[!ismissing]
		y <- y[!ismissing]
		labels <- labels[!ismissing]
	}
# other methods:
  id.var <- switch(id.method,
		x = if(log.x==TRUE)
          suppressWarnings(if(all(x) > 0)
				   abs(log(x) - mean(log(x))) else
           return(invisible(NULL)))  else
           abs(x - mean(x)),
		y = if(log.y==TRUE)
          suppressWarnings(if(all(y) > 0)
					 abs(log(y) - mean(log(y))) else
           return(invisible(NULL)))  else
           abs(y - mean(y)),
		r = if(log.y==TRUE)
					suppressWarnings(if(all(y) > 0)
					 abs(log(y)) else
					 return(invisible(NULL)))  else
					 abs(y),
    mahal = if(log.x == TRUE & log.y == TRUE) {
          suppressWarnings(if(all(x) > 0 & all(y) > 0)
					 rowSums( qr.Q(qr(cbind(1, log(x), log(y))))^2 ) else
           return(invisible(NULL))) } else {
            if(log.x == TRUE) {
             suppressWarnings(if(all(x) > 0 )
						 rowSums( qr.Q(qr(cbind(1, log(x), y)))^2 ) else
             return(invisible(NULL))) } else {
            if(log.y == TRUE) {
             suppressWarnings(if(all(y) > 0 )
							rowSums( qr.Q(qr(cbind(1, x, log(y))))^2 ) else
              return(invisible(NULL)))  } else {
              rowSums( qr.Q(qr(cbind(1, x, y)))^2 ) }}})
     }
# require id.n positive
  if(id.n <= 0L) return(invisible(NULL))
# criterion
  ind <-  order(id.var, decreasing=TRUE)[1L:min(length(id.var), id.n)]
# position, now depends on id.location (as of 5/16/2016)
  if (id.location != "avoid"){
      if(id.location == "lr"){
      mid <- mean(if(par("xlog")==TRUE) 10^(par("usr")[1:2]) else
                  par("usr")[1:2])
    	labpos <- c(4,2)[1+as.numeric(x > mid)]
      } else {
        mid <- mean(if(par("ylog")==TRUE) 10^(par("usr")[3:4]) else
          par("usr")[3:4])
        labpos <- c(3,1)[1+as.numeric(y > mid)]
      }
    # print
    	for (i in ind) {
    		text(x[i], y[i], labels[i], cex = id.cex, xpd = TRUE,
    			col = id.col, pos = labpos[i], offset = 0.25)}
  }
  else maptools::pointLabel(c(x[ind], x[ind]), c(y[ind], y[ind]),
                c(paste0(" ", labels[ind], " "), rep(" ", length(ind))),
                cex=id.cex, xpd=TRUE, col=id.col)
  if (any(as.character(ind) != labels[ind])) names(ind) <- labels[ind]
  result <- ind
  if (length(result) == 0) return(NULL) else return(result)
}










