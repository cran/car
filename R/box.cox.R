# Box-Cox power transformations, with automatic start (J. Fox)

# last modified 2 April 02 by J. Fox

bc<-function(x,p) box.cox(x,p)

box.cox<-function(x,p, start=0){
    # last modified 15 Dec 2000 by J. Fox
    min<-min(x, na.rm=TRUE)
    s<-if (missing(start) & (min <= 0)) nice(-min +.05*diff(quantile(x,c(.25,.75), na.rm=TRUE)), "up")
        else start
    if (missing(start) & s != 0) warning(paste("start = ", s, "added to data prior to transformation"))
    x<-x+s
    if (p==0) log(x)
        else (x^p-1)/p
    }
