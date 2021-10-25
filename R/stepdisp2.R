stepdisp2 = function(form,dadosd,d,alpha,modini,lambda,saidas){
  intercept=F
  nomes = attr(stats::terms(form), "term.labels")
  #assign(".dadosd",dadosd,envir=.GlobalEnv)
  .dadosd = dadosd

  modelo.principal = as.formula(paste("y~",modini))
  nomes.principais = attr(stats::terms(modelo.principal), "term.labels")

  ncov = ncol(dadosd)-1
  n = nrow(dadosd)


  ###############
  if (length(nomes.principais) > 1){
    model.folga = nomes.principais[1]
    if (length(nomes.principais)>2)
      for(i in 2:(length(nomes.principais)-1))   model.folga = paste(model.folga,"+",nomes.principais[i])

    mod.1 = as.formula("d~ 1")
    ob1d = glm(mod.1,family = Gamma (link=log),data=cbind(.dadosd,d))
    mod.termo = as.formula(paste("d~",model.folga))  ##
    obr = glm(mod.termo,family = Gamma (link=log),data=cbind(.dadosd,d))


    de1= ob1d$deviance
    de2 = obr$deviance
    q = length(summary(ob1d)$coef)/4
    p = length(summary(obr)$coef)/4
    Chi= 0.5*(de1-de2)
    pvalue = 1- pchisq(Chi,p-q)

    if (saidas) cat("Initial testing dispersion model pvalue: ",pvalue,"\n")

    ###############

    if (pvalue < alpha){
      ref =
        modatual = as.formula(paste("d~",modini))
      r2 = calcula.rsq(obr,ob1d,rep(1,n),lambda,"Gamma")
      ref = r2 #1-((n-1)/(n-lambda*length(nomes.principais))) * (1-r2)
    }
    else{
      modatual = mod.1
      r2 = calcula.rsq(ob1d,ob1d,rep(1,n),lambda,"Gamma")
      ref = 0 # 1-((n-1)/(n-lambda)) * (1-r2)
      modini = "1"
    }
  } else{
    modatual = as.formula(paste("d~",modini))
    ref = 0
  }

  R = data.frame(ncol=7,nrow=20)
  content = 0;
  j = 1
  repeat{
    j = j+1
    obj = glm(modatual,family = Gamma (link=log),data=cbind(.dadosd,d))
    refaux = ref
    ref = -9999999
    q = length(summary(obj)$coef)/4
    for (i in 1:ncov){
      modtest = as.formula(paste(modatual[2],"~",modatual[3],"+",colapse= nomes[i]))
      if (!intercept)
        if (j==2)  modtest = as.formula(paste(modatual[2],"~",modini,"+",colapse= nomes[i]))
      ob = glm(modtest,family = Gamma (link=log),data=cbind(.dadosd,d))
      phi1 = fitted.values(ob)
      #aic =  ob$aic
      rm = length(obj$coef)
      R2 = calcula.rsq(ob,ob1d,rep(1,n),lambda,"Gamma")
      #R2 = 1-((n-1)/(n-lambda*rm)) * (1-r2)

      if (R2 == refaux) R2 = ref

      if (R2 > ref) {
        ref = R2
        obref = ob
        coluna = i
      }


    }

    de1 = obj$deviance
    de2 = obref$deviance
    p = length(summary(obref)$coef)/4
    Chi= 0.5*(de1-de2)

    if (j==2) {
      R[1,1] = as.character(formula(obj))[3]
      r2 = calcula.rsq(obj,ob1d,rep(1,n),lambda,"Gamma")
      R[1,2] = round(r2,4) #round(1-((n-1)/(n-lambda*length(obj$coef))) * (1-r2),4)
      R[1,3] = round(de1,2)
      R[1,4] = NA
      R[1,5] = NA
      r2 = calcula.rsq(obj,ob1d,rep(1,n),1,"Gamma")
      R[1,6] = round(r2,4) #round(1-((n-1)/(n-length(obj$coef))) * (1-r2),4)
      R[1,7] = ""
    }
    R[j,1] =  as.character(formula(obref))[3]
    r2 = calcula.rsq(obref,ob1d,rep(1,n),lambda,"Gamma")
    R[j,2] =  round(r2,4) #round(1-((n-1)/(n-lambda*length(obref$coef))) * (1-r2),4) #    R[j,2] = 1-((n-1)/(n-rm*sqrt(n))  * (1-RQ1))
    r2 = calcula.rsq(obref,ob1d,rep(1,n),1,"Gamma")
    R[j,6] =  round(r2,4) #round(1-((n-1)/(n-length(obref$coef))) * (1-r2),4) #    R[j,6] =1-((n-1)/(n-rm)  * (1-RQ1))
    R[j,7] = ""
    R[j,3] = round(de2,2)
    R[j,4] = round(Chi,2)
    R[j,5] = round(1- pchisq(Chi,p-q),2)
    names(R) = c("Model","R2(f(n))","D","X^2","P","R2(1)","final")
    finaliza = FALSE
    if (intercept){
      if ((R[j,5] < alpha ))
        modatual = as.formula(paste(modatual[2],"~",modatual[3],"+",colapse= nomes[coluna]))
      else finaliza = TRUE
    }  else {
      if ((R[j,5] < alpha )){ #|| (j==2)) {
        if (j==2) modatual = as.formula(paste(modatual[2],"~",modini,"+",colapse= nomes[coluna]))
        else
          modatual = as.formula(paste(modatual[2],"~",modatual[3],"+",colapse= nomes[coluna]))
      }else finaliza = TRUE
    }


    if (finaliza) {
      modelofim = obj
      d = resid(obj,type="deviance")^2 /(1-hatvalues(obj))
      phi = fitted.values(obj)
      R[j-1,7] = "X"


      for (i in 2: nrow(R)){
        if (R[i,2] < R[i-1,2]){
          if (R[i,5]< alpha) modd = as.formula(paste("d~",R[i,1]))
          else modd = as.formula(paste("d~",R[i-1,1]))
          modelofim = glm(modd, family = Gamma(link=log),data=cbind(.dadosd,d))
          d = resid(modelofim,type="deviance")^2 /(1-hatvalues(modelofim))
          phi = fitted.values(modelofim)
          R[,7] = ""
          if (R[i,5]< alpha) R[i,7] = "X"
          else R[i-1,7] = "X"
          break
        }
      }

      #  if (j>2)
      #    for (i in 3:j){
      #    if (R[i,2] < R[i-1,2]) {
      #    modd = as.formula(paste(form[2],"~",R[i-1,1]))
      #    modelofim = glm(modd, family = "gaussian",data=.dad2,weights=1/phi)
      #    d = resid(modelofim,type="deviance")^2 /(1-hatvalues(modelofim))
      #     R[,7] = ""
      #     R[i-1,7] = "X"
      # }
      #  }


      break
    }
  }
  if (!intercept)
    #R2[2,3] = R2[2,4] =  R2[2,5] = NA

    #if (intercept) stepdisp = list(R2=R2,modelo=obj,phi=phi)

    #else stepdisp = list(R2=R2[-c(1),],modelo=obj,phi=phi)
    stepdisp = list(R=R,modelo=obj,phi=phi)
}
