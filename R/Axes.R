# Axes for transformations (J. Fox)

# last modified 2 April 02 by J. Fox

# function to find "nice" numbers

nice<-function(x, direction=c("round", "down", "up")){
    direction<-match.arg(direction)
    if (length(x)>1) return(sapply(x, nice, direction=direction))
    if (x==0) return(0)
    power.10<-floor(log(abs(x),10))
    lead.digit<-switch(direction,
        round=round(abs(x)/10^power.10),
        down=floor(abs(x)/10^power.10),
        up=ceiling(abs(x)/10^power.10))
    sign(x)*lead.digit*10^power.10
    }


# functions to add untransformed axis to right or top of a plot
#  for power or Box-Cox power transformations

power.axis<-function(power, base=exp(1), side=c("right", "above", "left", "below"), 
    at, grid=FALSE, grid.col=gray(.50), grid.lty=3,
    axis.title = "Untransformed Data", cex = 1, las=par("las")) {
    # last modified 20 Feb 2002 by J. Fox
    side<-if(is.numeric(side)) side 
        else which(match.arg(side)==c("below", "left", "above", "right"))
    axp<-if (side %% 2 == 1) par("xaxp") else par("yaxp")
    ticks<-nice(seq(from=axp[1], to=axp[2], length=axp[3]+1))
    ticks.x<- if (power !=0) nice(ticks[ticks>0]^(1/power)) 
        else nice(log(base)*exp(ticks))
    ticks.x <- if (missing(at)) ticks.x
        else at
    ticks.text <- as.character(ticks.x)
    ticks.trans<-if (power !=0) ticks.x^power else log(ticks.x, base)
    axis(side, labels = ticks.text, at = ticks.trans, las=las)
    if (grid & (side %% 2 == 0)) abline(h=ticks.trans, lty=grid.lty, col=grid.col)
    if (grid & (side %% 2 == 1)) abline(v=ticks.trans, lty=grid.lty, col=grid.col)
    mtext(axis.title, side = side, line = 3, cex = cex)
    }

box.cox.axis<-function(power, side=c("right", "above", "left", "below"), 
    at, grid=FALSE, grid.col=gray(.50), grid.lty=3,
    axis.title = "Untransformed Data", cex = 1, las=par("las")) {
    # last modified 20 Feb 2002 by J. Fox
    inverse.power<-function(x,p){
        if (p==0) exp(x)
        else (1+p*x)^(1/p)
        }
    side<-if(is.numeric(side)) side 
        else which(match.arg(side)==c("below", "left", "above", "right"))
    axp<-if (side %% 2 == 1) par("xaxp") else par("yaxp")
    ticks<-nice(seq(from=axp[1], to=axp[2], length=axp[3]+1))
    ticks.x<- if (power !=0) nice(inverse.power(ticks[ticks>0], power))
        else nice(inverse.power(ticks, 0))
    ticks.x <- if (missing(at)) ticks.x
        else at
    ticks.text <- as.character(ticks.x)
    ticks.trans<-box.cox(ticks.x,power)
    axis(side, labels = ticks.text, at = ticks.trans, las=las)
    if (grid & (side %% 2 == 0)) abline(h=ticks.trans, lty=grid.lty, col=grid.col)
    if (grid & (side %% 2 == 1)) abline(v=ticks.trans, lty=grid.lty, col=grid.col)
    mtext(axis.title, side = side, line = 3, cex = cex)
    }


# function to add a right or top probability axis to a plot of logits

prob.axis<-function(at, side=c("right", "above", "left", "below"),
    grid=FALSE, grid.lty=3, grid.col=gray(.50),
    axis.title = "Probability", interval = 0.1, cex = 1, las=par("las"))
{
    # last modified 20 Feb 2002 by J. Fox
    side<-if(is.numeric(side)) side 
        else which(match.arg(side)==c("below", "left", "above", "right"))
    logit<-if (side %% 2 == 1) par("usr")[c(1,2)] else par("usr")[c(3,4)]
    fact <- 10^( - (floor(log(interval, 10))))
    p.min <- nice(1/(1 + exp( - logit[1])), direction="down")
    p.max <- nice(1/(1 + exp( - logit[2])), direction="up")
    tick.min <- max(interval, (floor(fact * p.min))/fact)
    tick.max <- min(1 - interval, (ceiling(fact * p.max))/fact)
    ticks.p <- seq(tick.min, tick.max, interval)
    if(p.min <= 0.05) ticks.p <- c(0.05, ticks.p)
    if(p.min <= 0.01) ticks.p <- c(0.01, ticks.p)
    if(p.max >= 0.95) ticks.p <- c(ticks.p, 0.95)
    if(p.max >= 0.99) ticks.p <- c(ticks.p, 0.99)
    ticks.p<-if (missing(at)) ticks.p else at
    ticks.text <- as.character(ticks.p)
    ticks.logit <- log(ticks.p/(1 - ticks.p))
    axis(side, labels = ticks.text, at = ticks.logit, las=las)
    if (grid & (side %% 2 == 0)) abline(h=ticks.logit, lty=grid.lty, col=grid.col)
    if (grid & (side %% 2 == 1)) abline(v=ticks.logit, lty=grid.lty, col=grid.col)
    mtext(axis.title, side = side, line = 3, cex = cex)
}
