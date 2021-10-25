stepmean2 = function(form,dad2,Y,phi,alpha,familia,lambda,modini,moddisp,saidas){
  intercept = F

  nomes = attr(stats::terms(form), "term.labels")

  modelo.principal = as.formula(paste("y~",modini))
  nomes.principais = attr(stats::terms(modelo.principal), "term.labels")


  #assign(".dad2",dad2,envir=.GlobalEnv)
  #assign(".familia",familia, envir=.GlobalEnv)
  .dad2 = dad2
  .familia = familia

  #require(rsq)
  n = nrow(dad2)
  if (phi[1]==1) phi = rep(1,n)
  ncov = length(nomes)

  if (length(nomes.principais) > 1){
    model.folga = nomes.principais[1]
    if (length(nomes.principais)>2)
      for(i in 2:(length(nomes.principais)-1))   model.folga = paste(model.folga,"+",nomes.principais[i])


    mod.termo = as.formula(paste(form[2],"~",model.folga))  ##
    mod.1 = as.formula(paste(form[2],"~ 1"))
    obr = glm(mod.termo, family = .familia,data=.dad2,weights=1/phi)
    ob1 = glm(mod.1, family = .familia,data=.dad2,weights=1/phi)

    dest1 = resid(ob1,type="deviance")^2 /(1-hatvalues(ob1))
    de1 = sum(dest1)
    dest2 = resid(obr,type="deviance")^2 /(1-hatvalues(obr))
    de2 = sum(dest2)
    p = length(summary(obr)$coef)/4
    q = length(summary(ob1)$coef)/4
    F = ((de1-de2)/(de2/(n-p)))
    pvalue = 1- pf(F,p-q,n-p)
    d1 = (resid(obr,type="deviance")^2) /(1-hatvalues(obr))

    if (moddisp == "constante") npardisp = 1
    else {
      ajd = glm(moddisp,family="Gamma",data=cbind(.dad2,d1))
      npardisp = length(ajd$coef)}


    if (saidas) cat("Initial testing mean model Pvalue: ",pvalue,"\n")



    if (pvalue < alpha){
      modatual = as.formula(paste(form[2],"~",modini))
      #ref = calcula.rsq(obr,phi,lambda,familia)
      ref = calculaqde(obr$y,d1,phi,obr$family[1]$family,.dad2,n) * (length(obr$coef)+npardisp)
    }
    else{
      modatual = mod.1
      ref = 0
      modini = "1"
    }

  } else{
    modatual = as.formula(paste(form[2],"~",modini))
    ref = 0
  }


  R = data.frame(ncol=7,nrow=20)

  j = 1
  repeat{
    j = j+1
    obj = glm(modatual,family = .familia,data=.dad2,weights=1/phi)
    refaux = ref
    ref = 999999999
    q = length(summary(obj)$coef)/4
    for (i in 1:ncov){

      modtest = as.formula(paste(modatual[2],"~",modatual[3],"+",colapse= nomes[i]))
      if (!intercept)
        if (j==2)  modtest = as.formula(paste(modatual[2],"~",modini,"+",colapse= nomes[i]))
      ob = glm(modtest, family = .familia,data=.dad2,weights=1/phi)
      rm = length(ob$coef)
      d1 = (resid(ob,type="deviance")^2) /(1-hatvalues(ob))

      if (moddisp == "constante") npardisp = 1
      else {
        ajd = glm(moddisp,family="Gamma",data=cbind(.dad2,d1))
        npardisp = length(ajd$coef)}

      #RQ1 = summary(lm(formula(ob),data=dad2,weights=1/phi))$r.squared
      #R2 = calcula.rsq(ob,phi,lambda,familia)
      k = (length(ob$coef)+npardisp)
      AICE2K = calculaqde(ob$y,d1,phi,ob$family[1]$family,.dad2,n) + (2*k*n)/(n-k-1)

      # R2 = 1-((n-1)/(n-rm*sqrt(n))  * (1-RQ1))
      if (AICE2K == refaux) AICE2K = ref
      if (AICE2K < ref){
        ref = AICE2K
        obref = ob
        coluna = i

      }
    }


    dest1 = resid(obj,type="deviance")^2 /(1-hatvalues(obj))
    de1 = sum(dest1)
    dest2 = resid(obref,type="deviance")^2 /(1-hatvalues(obref))
    de2 = sum(dest2)
    p = length(summary(obref)$coef)/4
    F = ((de1-de2)/(de2/(n-p)))
    d1 = (resid(obj,type="deviance")^2) /(1-hatvalues(obj))
    if (j==2) {
      R[1,1] = as.character(formula(obj))[3]
      k = (length(obj$coef)+npardisp)
      R[1,2] = round(calculaqde(obj$y,d1,phi,obj$family[1]$family,.dad2,n) + (2*k*n)/(n-k-1) ,4)
      R[1,3] = round(de1,2)
      R[1,4] = 0
      R[1,5] = 0
      R[1,6] = round(calcula.rsq(obj,ob1,phi,1,familia),4)
      R[1,7] = ""
    }
    R[j,1] =  as.character(formula(obref))[3]

    # RQ1 = summary(lm(formula(obref),data=dad2),weights=1/phi)$r.squared
    #RQ1 = calcula.rsq(obref,phi,sqrt(n),"binomial")
    d1 = (resid(obref,type="deviance")^2) /(1-hatvalues(obref))
    k = (length(obref$coef)+npardisp)
    R[j,2] = round(calculaqde(obref$y,d1,phi,obref$family[1]$family,.dad2,n) + (2*k*n)/(n-k-1) ,4)#    R[j,2] = 1-((n-1)/(n-rm*sqrt(n))  * (1-RQ1))
    R[j,6] =  round(calcula.rsq(obref,ob1,phi,1,familia),4) #    R[j,6] =1-((n-1)/(n-rm)  * (1-RQ1))
    R[j,7] = ""
    # d1 = (resid(obref,type="deviance")^2) /(1-hatvalues(obref))
    # R[j,2] = calculaqde(obref$y,d1,phi,"gaussian",n) * (rm+lengthd)
    R[j,3] = round(de2,2)
    R[j,4] = round(F,2)
    if (p==q) R[j,5] = 0.9
    else
      R[j,5] = round(1- pf(F,p-q,n-p),4)
    names(R) = c("Model","EAIC","D","F","P","R2(1)","final")
    #if ( F > qf(alpha,p-q,n-p)){

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
      R[j-1,7] = "X"


      for (i in 2: nrow(R)){
        if (R[i,2] > R[i-1,2]){
          if (R[i,5]< alpha) modd = as.formula(paste(form[2],"~",R[i,1]))
          else modd = as.formula(paste(form[2],"~",R[i-1,1]))
          modelofim = glm(modd, family = .familia,data=.dad2,weights=1/phi)
          d = resid(modelofim,type="deviance")^2 /(1-hatvalues(modelofim))
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
  # if (!intercept) R[2,3] = R[2,4] = R[2,5] = NA
  # if (intercept) stepmean = list(R=R,modelo=obj,d=d)
  #  else stepmean = list(R=R[-c(1),],modelo=obj,d=d)
  list(R=R,modelo=modelofim,d=d)
}
