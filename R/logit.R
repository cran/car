# logit transformation of proportion or percent (J. Fox)

logit<-function(p, percents=max(p, na.rm=T)>1, adjust){   
    if (percents) p<-p/100
    a<-if (missing(adjust)) {
        if (min(p, na.rm=T)==0 | max(p, na.rm=T)==1) .025 else 0
        }
        else adjust
    if (missing(adjust) & a != 0) warning(paste("Proportions remapped to (",
        a,",",1-a,")", sep=""))
    a<-1-2*a
    log((.50+a*(p-.50))/(1-(.50+a*(p-.50))))
    }
    
