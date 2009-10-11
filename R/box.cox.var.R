# constructed variable for Box-Cox transformation (J. Fox)

# last modified 2 April 02 by J. FOx

box.cox.var<-function(y) {
    geo.mean<-exp(mean(log(y),na.rm=TRUE))
    y*(log(y/geo.mean) - 1)
    }
