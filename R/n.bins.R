# number of bins for histogram by Freedman-Diaconis rule (J. Fox)

n.bins<-function(x, rule=c("freedman.diaconis","sturges","scott","simple")){
    #last modified 16 Dec 2000 by J. Fox
    rule<-match.arg(rule)
    x<-x[!is.na(x)]
    n<-length(x)
    Q<-quantile(x, c(.25,.75))
    X<-range(x)
    result<-switch(rule,
        freedman.diaconis=ceiling(n^(1/3)*(X[2]-X[1])/(2*(Q[2]-Q[1]))),
        sturges=ceiling(log(n,2)+1),
        scott=ceiling((X[2]-X[1])*(n^(1/3))/(3.5*sqrt(var(x)))),
        simple=floor(if (n > 100) 10*log(n,10) else 2*sqrt(n))
        )
    names(result)<-NULL
    result
    }
 
