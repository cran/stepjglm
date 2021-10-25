stepdisp = function(form,dadosd,d,alpha,modini,saidas){
  intercept=F
  nomes = attr(stats::terms(form), "term.labels")


  modelo.principal = as.formula(paste("y~",modini))
  nomes.principais = attr(stats::terms(modelo.principal), "term.labels")

  ncov = ncol(dadosd)-1
  n = nrow(dadosd)


  ###############
  if (length(nomes.principais) > 1){
    model.folga = nomes.principais[1]
    if (length(nomes.principais)>2)
      for(i in 2:(length(nomes.principais)-1))   model.folga = paste(model.folga,"+",nomes.principais[i])


    mod.termo = as.formula(paste("d~",model.folga))  ##
    mod.1 = as.formula("d~ 1")
    ob1d = glm(mod.1,family = Gamma (link=log),data=cbind(dadosd,d))
    obr = glm(mod.termo,family = Gamma (link=log),data=cbind(dadosd,d))

    de1= ob1d$deviance
    de2 = obr$deviance
    q = length(summary(ob1d)$coef)/4
    p = length(summary(obr)$coef)/4
    Chi= 0.5*(de1-de2)
    pvalue = 1- pchisq(Chi,p-q)

    if (saidas) cat("Initial testing dispersion model Pvalue: ",pvalue,"\n")

    ###############

    if (pvalue < alpha){
      modatual = as.formula(paste("d~",modini))
      ref = obr$aic
    }
    else{
      modatual = mod.1
      ref = ob1d$aic
      modini = "1"
    }
  } else{
    modatual = as.formula(paste("d~",modini))
    ref = 0
  }

  R2 = data.frame(ncol=5,nrow=20)
  content = 0;
  j = 1
  repeat{
    j = j+1
    obj = glm(modatual,family = Gamma (link=log),data=cbind(dadosd,d))
    refaux = ref
    ref = 999999999
    q = length(summary(obj)$coef)/4
    for (i in 1:ncov){
      modtest = as.formula(paste(modatual[2],"~",modatual[3],"+",colapse= nomes[i]))
      if (!intercept)
        if (j==2)  modtest = as.formula(paste(modatual[2],"~",modini,"+",colapse= nomes[i]))
      ob = glm(modtest,family = Gamma (link=log),data=cbind(dadosd,d))
      phi1 = fitted.values(ob)
      aic =  ob$aic
      if (aic == refaux) aic = ref
      if (aic < ref){
        ref = ob$aic
        obref = ob
        coluna = i

      }
    }
    de1= obj$deviance
    de2 = obref$deviance
    p = length(summary(obref)$coef)/4
    Chi= 0.5*(de1-de2)
    if (j==2) {
      R2[1,1] =  as.character(formula(obj))[3]
      R2[1,2] = round(obj$aic,4)
      R2[1,3] = round(de1,2)
      R2[1,4] = NA
      R2[1,5] = NA
    }
    R2[j,1] =  as.character(formula(obref))[3]
    R2[j,2] = round(ref,4)
    R2[j,3] = round(de2,2)
    R2[j,4] = round(Chi,2)
    R2[j,5] = round(1- pchisq(Chi,p-q),2)
    names(R2) = c("Model","AIC","D","X^2","P")

    finaliza = FALSE

    if (intercept){
      if    ( ((p!=q) && 1- pchisq(Chi,p-q) < alpha ) ){  # (0.5*(de1-de2) > qchisq(1-alpha,p-q))
        #if (j==2) modatual = as.formula(paste(modatual[2],"~","-1+",colapse= nomes[coluna]))
        #else
        modatual = as.formula(paste(modatual[2],"~",modatual[3],"+",colapse= nomes[coluna]))
        content = content+1

      } else finaliza = TRUE
    }  else {

      if    (  ((p!=q) && (0.5*(de1-de2) > qchisq(1-alpha,p-q))) ){ # || (j==2) ){
        if (j==2) modatual = as.formula(paste(modatual[2],"~",modini,"+",colapse= nomes[coluna]))
        else
          modatual = as.formula(paste(modatual[2],"~",modatual[3],"+",colapse= nomes[coluna]))
        content = content+1
      } else finaliza = TRUE

    }

    if (finaliza) {
      phi = fitted.values(obj)
      break
    }
  }
  if (!intercept)
    #R2[2,3] = R2[2,4] =  R2[2,5] = NA

    #if (intercept) stepdisp = list(R2=R2,modelo=obj,phi=phi)

    #else stepdisp = list(R2=R2[-c(1),],modelo=obj,phi=phi)
    stepdisp = list(R2=R2,modelo=obj,phi=phi)
}
