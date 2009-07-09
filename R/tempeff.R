`tempeff` <-
function(formula, z, data, tcontrol=temp.control(), pcontrol=p.control(),
    fcontrol=fit.control(), etastart=NULL, ndx.seas=0, ...){
#usa tcontrol=NULL per non stimare la temp (comunque mi sa che z deve essere fornita)
#
#utilizza te1() +pen.opt()+bspline()
#modifica 14/11/05 per f; modifica 01/05/07; 17/05/07
#News: eliminata l'opzione f per sciegliere la base (adesso una vere bspline())
#------------------
#y: la risposta
#X: la matrice delle esplicative standard (escluso la temp)
#z: la temp
#psi: starting value. Se length(psi)=2 viene stimato un modello "\_/"
#L: un vettore che definisce i lag per il freddo e il caldo rispettivamente
#heat.power: esponente del caldo
#ndx: un vettore per passare argomenti a bspline() per la costruzione delle basi del freddo e del caldo.
#   considera che ndx è il n.intervalli, quindi ndx-1 sono i nodi interni e ndx+deg è il rango della base
#etastart=un vettore da passare a te1() come valori iniziali (attenzione alla lunghezza..)
#control: una lista con opzioni per la stima..(vedi fit.control()
#...: argomenti da passare a te1(), quali penalty, ridge.formulas, seasonP . Ad esempio:
#       penalty: list(DL=TRUE,diff.varying="no",ridge=TRUE)
#       ridge.formulas=list(cold="xlag^2",heat="xlag^2")
#
#NB: i primi max(L) osservazioni vengono scartate (of course). Quindi la dev che viene stampata se visual=T
#   è diversa da quella del modello originale in quanto mancano le prime max(L) osservazioni.
#si dovrebbe aggiungere l'argomento middle.null che se FALSE stima two-breakpoints model without constraints on the slopes

#--required functions: bspline() lagged()
bspline<-function(x, ndx, xlr=NULL, deg=3, deriv=0){
#x: vettore di dati
#xlr: il vettore di c(xl,xr)
#ndx: n.intervalli in cui dividere il range
#deg: il grado della spline
#Restituisce ndx+deg basis functions per ndx-1 inner nodi
require(splines)
    if(is.null(xlr)){
    xl<-min(x)-.01*diff(range(x))
    xr<-max(x)+.01*diff(range(x))
       } else {
         if(length(xlr)!=2) stop("quando fornito, xlr deve avere due componenti")
        xl<-xlr[1]
        xr<-xlr[2]
     }
    dx<-(xr-xl)/ndx
    knots<-seq(xl-deg*dx,xr+deg*dx,by=dx)
    B<-splineDesign(knots,x,ord=deg+1,derivs=rep(deriv,length(x)))
    #B<-spline.des(knots,x,bdeg+1,0*x) #$design
    B #the B-spline base matrix
}#end_fn

lagged<-function(x,lag=1){#by T.Lumley
        if (lag==0) return(x)
            n<-length(x)
            c(rep(NA,lag),x[-( (n-lag+1):n)])
      }
#-----------------------------------

    etic=deparse(substitute(z))
    if (missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m<-match(c("formula", "data", "z"), names(mf), 0)
    mf <- mf[c(1, m)]
    #get_all_vars(y~x+temp(z),data=d) potrebbe essere anche utile.
#    m <- match(c("formula", "data", "subset", "weights", "na.action",
#        "etastart", "mustart", "offset"), names(mf), 0)
#    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame") #restituisce model.frame(formula = y ~ x, data = d, drop.unused.levels = TRUE)
    mf <- eval(mf, parent.frame()) #restituisce una dataframe
    y <- model.response(mf, "any")
    mt <- attr(mf, "terms")
    X<-model.matrix(mt, mf, contrasts)
    #id.temp<-which(is.na(match(names(mf),c(names(mf)[1],colnames(X)))))
    #z<-mf[,id.temp]
    z<-mf[,length(all.vars(formula))+1]
    
    toll<-fcontrol$toll
    visual<-fcontrol$visual
    it.max<-fcontrol$it.max
    GLM<-fcontrol$GLM
    maxit.glm<-fcontrol$maxit.inner

    if(!is.na(tcontrol[1])){
        only.seas<-FALSE
        psi<- tcontrol$psi
        heat.power<-tcontrol$heat.power
        L <-tcontrol$L
        ndx <- tcontrol$ndx

        yy<-y[-(1:max(L))]
        XX<-X[-(1:max(L)),]

        k<-length(psi)
        if(any(psi<=min(z)) || any(psi>=max(z))) stop("Invalid starting values for psi")
        if(length(psi)>2) stop("One or two breakpoints allowed")
        if(length(psi)==1) psi<-c(psi,psi)
        psi.new<-psi
        n<-length(yy)
        PSI1 <- matrix(rep(rep(psi.new[1],(L[1]+1)), rep(n, L[1]+1)), ncol = L[1]+1)
        PSI2 <- matrix(rep(rep(psi.new[2],(L[2]+1)), rep(n, L[2]+1)), ncol = L[2]+1)
        Xlag<-as.matrix(sapply(0:max(L),function(i) lagged(z,i)))[-(1:max(L)),]
        U1<-pmin(Xlag[,1:(L[1]+1)]-PSI1,0)
        U2<-pmax(Xlag[,1:(L[2]+1)]-PSI2,0)^heat.power

        A1<-bspline(0:L[1],ndx=ndx[1]) #B-spline basis for cold
        A2<-bspline(0:L[2],ndx=ndx[2]) #B-spline basis for heat
        colnames(A1)<-colnames(A2)<-NULL
        Xf<-U1%*%A1
        Xc<-U2%*%A2
        colnames(Xf)<-paste("Xf",1:ncol(Xf),sep="")
        colnames(Xc)<-paste("Xc",1:ncol(Xc),sep="")
    #if(!middle.null){
        #Umezzo<-pmax(Xlag[,1:(L[2]+1)]-PSI1,0)
        #Xmezzo<-Umezzo%*%A2
        #colnames(Xmezzo)<-paste("Xmm",1:ncol(Xmezzo),sep="")
        #XREG<-cbind(XREG,Xmezzo)
        #}
        } else {
        only.seas<-TRUE
        if(ndx.seas==0) stop("tcontrol=NULL, please provide 'ndx.seas>0'")
        o<-tempeff.fit(y,X,etastart=etastart,only.seas=only.seas,ndx.seas=ndx.seas,...)
        o$call<-match.call()
        return(o)
        }
    if(it.max==0){
    #se it.max=0 restituisce un modello con assegnati breakpoints==psi
        o<-tempeff.fit(yy,XX,Af=A1,Ac=A2,Xf,Xc,V=NULL,etastart=etastart,only.seas=only.seas,...)
        o$call<-match.call()
        return(o)
#        COV<-solve(crossprod(XREG,diag(obj$weights))%*%XREG)
#        colnames(COV)<-rownames(COV)
#        varFr<-A1%*%COV[id1, id1]%*%t(A1)
#        varCa<-A2%*%COV[id2, id2]%*%t(A2)
        #if(!middle.null){
            #varFr<-A1%*%COV[id1, id1]%*%t(A1)
            #varM<-A2%*%COV[idm, idm]%*%t(A2)
            #AC<-cbind(A2,A2) #questa è la matrice che serve per ottenere le slope dopo psi2 (in realtà è cbind(Am,A2))
            #betaCa<-AC%*%c(am,a2) #queste sono le slope dopo psi2
            #varCa<-AC%*%COV[c(idm,id2),c(idm,id2)]%*%t(AC) #questa è la var(slope dopo psi2)
            #sumcoefFr<- crossprod(rep(1,length(beta1)),beta1)
            #sumcoefM<- crossprod(rep(1,length(betam)),betam)
            #sumcoefCa<- crossprod(rep(1,length(betaCa)),betaCa)
            #VarSumcoefFr<- crossprod(rep(1,length(beta1)),varFr)%*%rep(1,length(beta1))
            #VarSumcoefM<- crossprod(rep(1,length(betam)),varM)%*%rep(1,length(betam))
            #VarSumcoefCa<- crossprod(rep(1,length(betaCa)),varCa)%*%rep(1,length(betaCa))
            #Fr<-cbind(Est=c(beta1,sumcoefFr),SE=c(sqrt(diag(varFr)),sqrt(VarSumcoefFr)))
            #Fr<-cbind(Fr,tvalue=Fr[,1]/Fr[,2])
            #rownames(Fr)<-c(paste("lag",0:L[1],""),"sum")
            #Me<-cbind(Est=c(betam,sumcoefM),SE=c(sqrt(diag(varM)),sqrt(VarSumcoefM)))
            #Me<-cbind(Me,tvalue=Me[,1]/Me[,2])
            #rownames(Me)<-c(paste("lag",0:L[2],""),"sum")
            #Ca<-cbind(Est=c(betaCa,sumcoefCa),SE=c(sqrt(diag(varCa)),sqrt(VarSumcoefCa)))
            #Ca<-cbind(Ca,tvalue=Ca[,1]/Ca[,2])
            #rownames(Ca)<-c(paste("lag",0:L[2],""),"sum")
            #out<-list(modello=obj,cold=Fr,middle=Me,heat=Ca)
            #return(out)
            #}
#        sumcoefFr<- crossprod(rep(1,length(beta1)),beta1)
#        sumcoefCa<- crossprod(rep(1,length(beta2)),beta2)
#        VarSumcoefFr<- crossprod(rep(1,length(beta1)),varFr)%*%rep(1,length(beta1))
#        VarSumcoefCa<- crossprod(rep(1,length(beta2)),varCa)%*%rep(1,length(beta2))
#        Fr<-cbind(Est=c(beta1,sumcoefFr),SE=c(sqrt(diag(varFr)),sqrt(VarSumcoefFr)))
#        Fr<-cbind(Fr,tvalue=Fr[,1]/Fr[,2])
#        rownames(Fr)<-c(paste("lag",0:L[1],""),"sum")
#        Ca<-cbind(Est=c(beta2,sumcoefCa),SE=c(sqrt(diag(varCa)),sqrt(VarSumcoefCa)))
#        Ca<-cbind(Ca,tvalue=Ca[,1]/Ca[,2])
#        rownames(Ca)<-c(paste("lag",0:L[2],""),"sum")
#        out<-list(modello=obj,cold=Fr,heat=Ca)
#        return(out)
        } #end_ if(itmax=0)
    dev00<-glm.fit(x=XX, y=yy, family=poisson())$dev
    XREG<-cbind(XX,Xf,Xc)
    obj<-glm.fit(XREG, yy, family=poisson())
    id1<-grep("Xf",names(coef(obj)))
    id2<-grep("Xc",names(coef(obj)))
    a1<-coef(obj)[id1]
    a2<-coef(obj)[id2]
    beta1<- A1%*%a1
    beta2<- A2%*%a2
    #if(!middle.null){
        #idm<-grep("Xmm",names(coef(obj)))
        #am<-coef(obj)[idm]
        #betam<- A2%*%am
        #}
    initial <- psi
    it <- 1
    epsilon <- 10
    #X<-X[-(1:max(L)),]
        while(abs(epsilon)>toll) { #start_while
            eta0<-obj$linear.predictors
            U1<-pmin(Xlag[,1:(L[1]+1)]-PSI1,0)
            U2<-pmax(Xlag[,1:(L[2]+1)]-PSI2,0)^heat.power
            Xf<-U1%*%A1
            Xc<-U2%*%A2
            colnames(Xf)<-paste("Xf",1:ncol(Xf),sep="")
            colnames(Xc)<-paste("Xc",1:ncol(Xc),sep="")
            V1<-ifelse(Xlag[,1:(L[1]+1)] <PSI1, -1, 0)
            V2<-ifelse(Xlag[,1:(L[2]+1)] >PSI2, -1, 0)
            if(heat.power!=1) V2<-2*pmax(Xlag[,1:(L[2]+1)]-PSI2,0)*V2
            V1<-t(t(V1)*as.vector(beta1))
            V2<-t(t(V2)*as.vector(beta2))
            V1<-rowSums(V1)
            V2<-rowSums(V2)
            V<- if(k==1) V1+V2 else cbind(V1,V2)
            XREG<-cbind(XX,Xf,Xc,V)
            dev.old<-obj$dev
            #----stima il modello..
            if(GLM){
              obj<-suppressWarnings(glm.fit(XREG, yy, family=poisson(),
                  etastart=eta0,control=glm.control(maxit=maxit.glm)))
              a1<-obj$coef[grep("Xf",names(obj$coef))]
              a2<-obj$coef[grep("Xc",names(obj$coef))]
              beta1<- A1%*%a1
              beta2<- A2%*%a2
              d<-obj$coef[grep("V",names(obj$coef))]
              } else {
                obj<-tempeff.fit(yy,XX,Af=A1,Ac=A2,Xf,Xc,V=V,gam.fit.it=maxit.glm,etastart=eta0,ndx.seas=ndx.seas,
                  penalty=pcontrol, ...)
                beta1<-obj$betaFreddo
                beta2<-obj$betaCaldo
                d<- obj$delta
                }
            #----
            #controlla bene nel caso di 1-2 psi
            psi.old <- psi
            psi <- psi.old + d
#            if(psi<=min(z) || psi>=max(z)) stop("psi fuori dal range")
            #PSI <- matrix(rep(psi.new, rep(nrow(x), ncol(x))), ncol = ncol(x))
            PSI1 <- matrix(rep(rep(psi[1],(L[1]+1)), rep(n, L[1]+1)), ncol = L[1]+1)
            PSI2 <- matrix(rep(rep(psi[2],(L[2]+1)), rep(n, L[2]+1)), ncol = L[2]+1)
            dev.new <- obj$dev #if(GLM) obj$dev else obj$dev
        if (visual) {
            flush.console()
            if (it == 1)
                cat(0, " ", formatC(dev00, 3, format = "f"),
                  "", "-----without variable", "\n")
            spp <- if (it < 10) "" else NULL
            cat(it, spp, "", formatC(c(dev.new,unique(psi)), 3, format = "f"), "\n")
            }
        #epsilon <- (dev.new - dev.old)/dev.old #
        epsilon <- abs(dev.new-dev.old)/(abs(dev.new)+0.1)
        it <- it + 1
        if (it > it.max)
            break
        } #end_while
    if (it > it.max) warning("max number of iterations attained", call. = FALSE)
    obj<-tempeff.fit(yy,XX,Af=A1,Ac=A2,Xf,Xc,V=V,etastart=eta0,ndx.seas=ndx.seas,penalty=pcontrol,...)
    var.psi<- obj$Ve[obj$id.d,obj$id.d]
    var.psi.bayes<-  obj$Vp[obj$id.d,obj$id.d]
    if(k==2) {
        var.psi<-diag(var.psi)
        var.psi.bayes<-diag(var.psi.bayes)
          }
    psi<-cbind(unique(psi),sqrt(var.psi),sqrt(var.psi.bayes))
    colnames(psi)<-c("Est","SE.freq","SE.bayes")
    obj$psi<-psi
    obj$call<-match.call()
    class(obj)<-"modTempEff"
    return(obj)
    } #end_function..

