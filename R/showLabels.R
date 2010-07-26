# last modified 25 Februrary 2010 by J. Fox
# rewritten 15 April 2010 S Weisberg

showLabels <- function(x, y, labels=NULL, id.method="identify",  
  id.n = length(x), id.cex=1, id.col=palette()[1],  ...) {
  if(id.n <= 0L) return(invisible(NULL))
  res <- NULL
  id.method <- if(is.list(id.method)) id.method else list(id.method)
  for (meth in id.method) 
     res <- c(res, showLabels1(x, y, labels, meth, id.n, id.cex, 
              id.col,  ...))
  return(if(is.null(res)) invisible(res) else res)
  }   

showLabels1 <- function(x, y, labels=NULL, id.method="identify",
	id.n = length(x), id.cex=1, id.col=palette()[1],  ...) {       
# If labels are NULL, try to get the labels from x:
  if (is.null(labels)) 
    labels <- names(x)
	if (is.null(labels))
		labels <- paste(seq_along(x))
# id.method can be a character string like "x" or "y", or it can be
# a vector like abs(rstudent(model)) or c(1, 2, 4) or a vector of labels
  use.built.in.method <- is.character(id.method) & length(id.method) == 1
  if(use.built.in.method==TRUE) 
	  match.arg(id.method, c("mahal", "x", "y", "identify")) 
# label color
	if (is.null(id.col))
		id.col <- palette()[1]
# Use identify?
  if(use.built.in.method==TRUE){  
   if(id.method == "identify") { 
    	result <- labels[identify(x, y, labels, n=length(x), cex=id.cex, 
                 col=id.col, ...)]
    	if(length(result) > 0) return(unique(result)) else return(NULL)
 	}}
# require id.n positive
  if(id.n <= 0L) return(invisible(NULL))  	
# logged-axes?
  log.x <- par("xlog")
  log.y <- par("ylog") 
# id.method might be a vector of labels or row numbers. 
  id.var <- NULL 
  if(use.built.in.method==FALSE) {
    if(length(id.method) == length(x)) id.var <- id.method else {
      id.var <- rep(0, length(x))
      names(id.var) <- labels
      id.var[id.method] <- 1
      id.n <- length(id.method)
      }}
# missing values need to be removed   
	ismissing <- is.na(x) | is.na(y) | is.na(labels) 
	if( any(ismissing) ) {
		x <- x[!ismissing]
		y <- y[!ismissing]
		labels <- labels[!ismissing]
		id.var <- id.var[!ismissing]
	}
# construct id.var if use.built.in.method==TRUE
   if(use.built.in.method==TRUE) {
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
# criterion
  ind <-  order(-id.var)[1L:id.n]
  mid <- mean(if(par("xlog")==TRUE) 10^(par("usr")[1:2]) else 
              par("usr")[1:2])
	labpos <- c(4,2)[1+as.numeric(x > mid)]
	for (i in ind) {
		text(x[i], y[i], labels[i], cex = id.cex, xpd = TRUE,
			col = id.col, pos = labpos[i], offset = 0.25, ...)} 
	result <- labels[ind]
	if (length(result) == 0) return(NULL) else return(result)
}






