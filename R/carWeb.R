# March 20, 2017  CarWeb needs revision so that it will
#                 go to website3 if page = website
#                 view errata in browser if page=errata
#                 go to taskviews if page = taskview
#                 download a file if file=filename
#                 download a script file if script = chap-num
#  Add more to carWeb including cheat sheets
# 2/21/2018 deleted "ethics" from the deafulat for page
# 2018-04-25: J. Fox. Update website URLs; update setup files
# 2018-04-28: J. Fox. Check whether file exists before overwriting


carWeb <-
function (page = c("webpage", "errata", "taskviews"), script, data, setup)
{
    rstudiocheat <- "https://www.rstudio.com/resources/cheatsheets/"
    ide.cheat <- "https://www.rstudio.com/wp-content/uploads/2016/01/rstudio-IDE-cheatsheet.pdf"
    data.page <- "https://socialsciences.mcmaster.ca/jfox/Books/Companion/data/"
    setup.dir <- "https://socialsciences.mcmaster.ca/jfox/Books/Companion/setup/"
    files <- c("Duncan.txt", "Duncan.csv", "Duncan.xlsx", "Duncan.Rmd", "Hamlet.txt", "RMarkdownTest.Rmd", 
               "zipmod.R", "zipmodBugged.R", "zipmod-generic.R", paste0("chap-", 1:10, ".R"))
    script.page <- "https://socialsciences.mcmaster.ca/jfox/Books/Companion/scripts/"
    ethics <- "http://www.amstat.org/asa/files/pdfs/EthicalGuidelines.pdf"
    page = match.arg(page)
    urls = c(webpage = "https://socialsciences.mcmaster.ca/jfox/Books/Companion/",
        errata = "https://socialsciences.mcmaster.ca/jfox/Books/Companion/errata.html",
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
    if(!missing(setup) && isTRUE(setup)){
      downloaded <- character(0)
      for(f in files) {
        if (file.exists(f)){
          response <- askYesNo(paste0(f, " exists, replace?"), prompts=c("yes", "no", "cancel"), default=FALSE)
          if (is.na(response)) {
            if (length(downloaded) > 0) cat("\nFiles downloaded:", paste(downloaded, collapse=", "), "\n")
            return(invisible(response))
          }
        }
        else response <- TRUE
        if (isTRUE(response)) {
          download.file(paste(setup.dir, f, sep=""),
                                              paste(getwd(), f, sep="/"))
          downloaded <- c(downloaded, f)
        }
      }
      if (length(downloaded) > 0) cat("\nFiles downloaded:", paste(downloaded, collapse=", "), "\n")
      return(invisible(NULL))
    }
    browseURL(url)
}
