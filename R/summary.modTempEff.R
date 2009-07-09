`summary.modTempEff` <-
function(object, spar=TRUE, digits = max(3, getOption("digits") - 3), ...){
      n<-length(object$fitted.values)
      dev<-object$dev
      tot.edf<-sum(object$edf)
      bic<-object$aic-2*tot.edf+log(n)*tot.edf
      ubre<- (object$dev + 2*tot.edf -1)/n
      edf.cold.tot<-sum(object$edf.cold)
      edf.heat.tot<-sum(object$edf.heat)
      is.ridge<-suppressWarnings(!is.na(object$call[["pcontrol"]][["ridge.formulas"]][[1]]))
      if(is.ridge) {
        ridgeC<-object$call[["pcontrol"]][["ridge.formulas"]][["cold"]]
        ridgeH<-object$call[["pcontrol"]][["ridge.formulas"]][["cold"]]}

    	seasP<-if(is.null(object$edf.seas)) { NA }
          else {paste("edf =", round(sum(object$edf.seas),2), "( rank =",object$rank.seas,")")}
      #print(object$call)
      cat("Model:\t")
      print(object$call[["formula"]])
      cat("Temperature:", object$call[["z"]], "\t")
	   if(!is.ridge) {cat("ridge penalty: ", NA, "\n")} else {
      cat("ridge penalty: ", "cold=", ridgeC, "heat=", ridgeH,"\n")
          }
	   cat("Seasonality: ", seasP,"\n")
	   cat("Fit summary:\n")
      cat("AIC =",object$aic," BIC =",bic," ubre =",ubre," dev =" ,object$dev, "\n")
      xx<-matrix(,2,6)
      rownames(xx)<-c("Cold","Heat")
      colnames(xx)<-c("Est","SE.freq","SE.bayes","rank","edf","L")
      xx["Cold",1:2]<-object$ToTfreddo[1:2]; xx["Cold",3]<-object$ToTfreddo.bayes[2]
      xx["Heat",1:2]<-object$ToTcaldo[1:2]; xx["Heat",3]<-object$ToTcaldo.bayes[2]
      xx["Cold","rank"]<-object$rank.cold
      xx["Heat","rank"]<-object$rank.heat
      xx["Cold","edf"]<-edf.cold.tot; xx["Heat","edf"]<-edf.heat.tot
      xx["Cold","L"]<-length(object$betaFreddo)-1
      xx["Heat","L"]<-length(object$betaCaldo)-1
      cat("Net effects of temperature: \n")
      print(xx,digits=digits)
      cat("Threshold: \n")
      psi<-object$psi
      rownames(psi)<-if(nrow(object$psi)==1) {""} else {c("psi1","psi2")}
      print(psi,digits=3)
      cat("V variable(s):\n")
      print(object$Tdelta)
      if(spar){
      cat("log Smoothing parameters: \n")
      colnames(object$sp.mio)<-""
      print(t(log(object$sp.mio)),digits=digits)
        }
      }

