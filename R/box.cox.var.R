# constructed variable for Box-Cox transformation (J. Fox)

box.cox.var<-function(y) {
    geo.mean<-exp(mean(log(y),na.rm=T))
    y*(log(y/geo.mean) - 1)
    }
