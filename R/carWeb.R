# March 20, 2017  CarWeb needs revision so that it will
#                 go to website3 if page = website
#                 view errata in browser if page=errata
#                 go to taskviews if page = taskview
#                 download a file if file=filename
#                 download a script file if script = chap-num
#  Add more to carWeb including cheat sheets
# 2/21/2018 deleted "ethics" from the default for page
# 2018-04-25: J. Fox. Update website URLs; update setup files
# 2018-04-28: J. Fox. Check whether file exists before overwriting
# 2024-09-19: J. Fox. Update URLs.
# 2024-09-23 Deprecated most of the carWeb features; updated .Rd page

carWeb <- function (...){ 
  browseURL(webpage <- "https://www.john-fox.ca/Companion/index.html")
  }

