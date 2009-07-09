`plot.modTempEff` <-
function(x, which=c("cold","heat"), new=TRUE, var.bayes=FALSE, add=FALSE, delta.rr=TRUE, level=0.95,...){
    which<-match.arg(which, several.ok = TRUE)
    cold<-match("cold",which,nomatch = 0)>0
    heat<-match("heat",which,nomatch = 0)>0
    if(add) new<-FALSE
    if(new) {x11();par(mfrow=c(1,cold+heat))}
    z<-qnorm((1+level)/2)
    xf<- -x$betaFreddo
    xc<-x$betaCaldo
    if(var.bayes) {
        SE.f<-x$SE.f.bayes
        SE.c<-x$SE.c.bayes
        } else {
          SE.f<-x$SE.f
          SE.c<-x$SE.c
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

