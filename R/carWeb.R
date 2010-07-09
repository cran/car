carWeb <-
function (page = c("webpage", "errata", "taskviews"), rfile, data)
{
    data.page <- "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/data/"
    rfile.page <- "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/rfile/"
    page = match.arg(page)
    urls = c(webpage = "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/",
        errata = "http://socserv.socsci.mcmaster.ca/jfox/Books/Companion/errata.pdf",
        taskviews = "http://cran.r-project.org/web/views")
    url = urls[page]
    if(!missing(rfile)) url <- paste(rfile.page, rfile, ".R",sep="")
    if(!missing(data)) url <- paste(data.page, data, ".txt", sep="")
    browseURL(url)
}