# change an argument to a function interactively (J. Fox)

Ask<-function(arg, fun, ...){ 
    fun<-fun
    repeat{   
        value<-readline(paste("Enter",deparse(substitute(arg)),": "))
        if (value == "") break()
        eval(parse(text=paste("fun(",deparse(substitute(arg)),"=",value,",...)")))
        }
    }
 
