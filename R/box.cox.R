# Box-Cox power transformations, with automatic start (J. Fox)

# last modified 27 April 04 by J. Fox

bc<-function(x,p) box.cox(x,p)

box.cox<-function(x,p, start=0){
    min<-min(x, na.rm=TRUE)
    s<-if (missing(start) & (min <= 0)) {
        IQR <- diff(quantile(x, c(.25,.75), na.rm=TRUE))
        if (IQR <= 10*.Machine$double.eps) stop("First and third quartile are equal.\nAutomatic start cannot be computed.")
        nice(-min +.05*IQR, "up")
        }
        else start
    if (missing(start) & s != 0) warning(paste("start = ", s, "added to data prior to transformation"))
    x<-x+s
    if (p==0) log(x)
        else (x^p-1)/p
    }
