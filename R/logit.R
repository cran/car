# logit transformation of proportion or percent (J. Fox)

# last modified 2 April 02

logit<-function(p, percents=max(p, na.rm=TRUE)>1, adjust){   
    if (percents) p<-p/100
    a<-if (missing(adjust)) {
        if (min(p, na.rm=TRUE)==0 | max(p, na.rm=TRUE)==1) .025 else 0
        }
        else adjust
    if (missing(adjust) & a != 0) warning(paste("Proportions remapped to (",
        a,",",1-a,")", sep=""))
    a<-1-2*a
    log((.50+a*(p-.50))/(1-(.50+a*(p-.50))))
    }
    
