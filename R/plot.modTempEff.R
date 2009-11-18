`plot.modTempEff` <-
function(x, which=c("cold","heat"), add=FALSE, new=TRUE, var.bayes=FALSE, delta.rr=TRUE, level=0.95,...){
    if(length(x$ToTheat)<=0){
        if(length(x$fit.seas)<=0) {stop("the model does not include csdl() or seas() terms") 
          } else {
          if(add) lines(exp(x$fit.seas),...) else plot(exp(x$fit.seas), 
            xlab="Time", ylab="Fitted Values", type="l", ...)
          return(invisible(x))
        }
      }
    which<-match.arg(which, several.ok = TRUE)
    cold<-match("cold",which,nomatch = 0)>0
    heat<-match("heat",which,nomatch = 0)>0
    if(add) {
      new<-FALSE
      if((cold+heat)>=2) stop("add=TRUE works only with a *single* DL curve")
      }
    if(new) {x11();par(mfrow=c(1,cold+heat))}
    z<-qnorm((1+level)/2)
    xf<- -x$betaCold
    xc<-x$betaHeat
    if(var.bayes) {
        SE.f<-x$SE.c.bayes
        SE.c<-x$SE.h.bayes
        } else {
          SE.f<-x$SE.c
          SE.c<-x$SE.h
               }
    if(cold){
      etich<-"logRR for Cold"
      rrr<- cbind(xf,xf-z*SE.f,xf+z*SE.f)
      if(!add){
          if(delta.rr) {
              etich<-"% change in mortality for Cold"
              rrr<-100*(exp(rrr)-1)}
          matplot(0:(length(xf)-1),rrr,
          type="l",lty=c(1,2,2),col=1,ylab=etich, xlab="Lag (day)",...)
          abline(h=0,lty=1,col=gray(.8))} else {
        matlines(0:(length(xf)-1),rrr,col=1,...)
      } }
    if(heat){
      etich<-"logRR for Heat"
      rrr<-cbind(xc,xc-z*SE.c,xc+z*SE.c)
      if(!add){
        if(delta.rr) {
            rrr<-100*(exp(rrr)-1)
            etich<-"% change in mortality for Heat"
            }
        matplot(0:(length(xc)-1),rrr,type="l",
          lty=c(1,2,2),col=1,ylab=etich,xlab="Lag (day)",...)
        abline(h=0,col=gray(.8))} else {
          matlines(0:(length(xc)-1),rrr,col=1,...)}
      }
    } #end_funct

