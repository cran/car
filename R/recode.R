# recode function (J. Fox)
# last modified 2 April 2002

recode<-function(var, recodes, as.factor.result){
    recode.list<-rev(strsplit(gsub(" ","",recodes),";")[[1]])
    is.fac<-is.factor(var)
    if (missing(as.factor.result)) as.factor.result <- is.fac
    if (is.fac) var<-as.character(var)
    result<-var
    if (is.numeric(var)) {
        lo<-min(var, na.rm=TRUE)
        hi<-max(var, na.rm=TRUE)
        }
    for (term in recode.list){
        if (0<length(grep(":", term))) {
            range<-strsplit(strsplit(term, "=")[[1]][1],":")
            low<-eval(parse(text=range[[1]][1]))
            high<-eval(parse(text=range[[1]][2]))
            target<-eval(parse(text=strsplit(term, "=")[[1]][2]))
            result[(var>=low)&(var<=high)]<-target
            }
        else if (0<length(grep("else", term))) {
            target<-eval(parse(text=strsplit(term, "=")[[1]][2]))
            result[1:length(var)]<-target
            }
        else {
            set<-eval(parse(text=strsplit(term, "=")[[1]][1]))
            target<-eval(parse(text=strsplit(term, "=")[[1]][2]))
            for (val in set){
                if (is.na(val)) result[is.na(var)]<-target
                    else result[var==val]<-target
                }
            }
        }
    if (as.factor.result) result<-as.factor(result)
        else if (!is.numeric(result)) {
            result.valid <- na.omit(result)
            if (length(result.valid) == length(grep("[0-9]",result.valid)))
                result <- as.numeric(result)
            }
    result
    }
 
