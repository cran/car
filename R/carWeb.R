# March 20, 2017  CarWeb needs revision so that it will
#                 go to website3 if page = website
#                 view errata in browser if page=errata
#                 go to taskviews if page = taskview
#                 download a file if file=filename
#                 download a script file if script = chap-num
#  Add more to carWeb including cheat sheets
# 2/21/2018 deleted "ethics" from the deafulat for page


carWeb <-
function (page = c("webpage", "errata", "taskviews"), script, data, setup)
{
    rstudiocheat <- "https://www.rstudio.com/resources/cheatsheets/"
    ide.cheat <- "https://www.rstudio.com/wp-content/uploads/2016/01/rstudio-IDE-cheatsheet.pdf"
    data.page <- "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/data/"
    files <- c("Duncan.txt", "Hamlet.txt", "Prestige.txt", "Prestige-bugged.txt", "Datasets.xls")
    script.page <- "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/scripts/"
    ethics <- "http://www.amstat.org/asa/files/pdfs/EthicalGuidelines.pdf"
    page = match.arg(page)
    urls = c(webpage = "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/",
        errata = "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/errata.html",
        taskviews = "http://cran.r-project.org/web/views",
        ethics = ethics)
    url <- urls[page]
    if(!missing(data)) {
       dfile <- unlist(strsplit(data, ".", fixed=TRUE))
       if(length(dfile) > 1) dfile <- dfile[1:(length(dfile)-1)]
       dfile <- paste(c(dfile, "txt"), collapse="." )
       url <- paste(data.page, dfile, sep="")}
    if(!missing(script)) {
       sfile <- unlist(strsplit(script, ".", fixed=TRUE))
       if(length(sfile) > 1) sfile <- sfile[1:(length(sfile)-1)]
       sfile <- paste(c(sfile, "R"), collapse="." )
       url <- paste(script.page, sfile, sep="")}
    if(!missing(setup)){
      for(f in files) download.file(paste(data.page, f, sep=""),
                                    paste(getwd(), f, sep="/"))
      return(NULL)
    }
    browseURL(url)
}
