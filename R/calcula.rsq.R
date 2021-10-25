#
calcula.rsq = function(modelom,modelo1,phi,lambda,familye){
  
  yfitm <- modelom$fitted.values
  rm = length(modelom$coef)
  n = length(modelom$y)
  if (modelom$family$family == "poisson"){
    SSm = sum((1+phi^2)*(modelom$y-yfitm)^2)/sum((1+phi^2)*(modelom$y-mean(modelom$y))^2)
    r2mm = 1-((n-1)/(n-lambda*rm)) * (SSm)
  }
  else if (modelom$family$family == "gaussian") {
    
    
    # r2 = summary(lm(formula(modelom),data=modelom$data,weights=1/phi))$r.squared
    
    r2 = rsq(lm(formula(modelom),data=modelom$data,weights=1/phi),adj=F)
    r2mm = 1-((n-1)/(n-lambda*rm)) * (1-r2)
  }
  
  else if (modelom$family$family == "binomial"){  
    a = 1-2*modelom$y
    b1 = 1-2*yfitm
    #b2 = 1-2*mean(modelom$y)
    b2 = 1-2*modelo1$fitted.values
    dv1 = sum((1/16*phi^2)*(log((phi*b1+sqrt(1+phi^2*b1^2))/(phi*a+sqrt(1+phi^2*a^2)))+ (phi*b1+sqrt(1+phi^2*b1^2)) - (phi*a+sqrt(1+phi^2*a^2)) )^2)
    dv2 = sum((1/16*phi^2)*(log((phi*b2+sqrt(1+phi^2*b2^2))/(phi*a+sqrt(1+phi^2*a^2)))+ (phi*b2+sqrt(1+phi^2*b2^2)) - (phi*a+sqrt(1+phi^2*a^2)) )^2)
    SSm = dv1/dv2
    r2mm = 1-((n-1)/(n-lambda*rm)) * (SSm)
  }
  else if (modelom$family$family == "Gamma"){  
    a = modelom$y
    b1 = yfitm
    #b2 = mean(modelom$y)
    b2 = modelo1$fitted.values
    dv1 = sum((1/(16*phi^2))*( log(aux(phi,b1)/aux(phi,a)) + aux(phi,b1) - aux(phi,a))^2 )
    dv2 = sum((1/(16*phi^2))*( log(aux(phi,b2)/aux(phi,a)) + aux(phi,b2) - aux(phi,a))^2 )
    SSm = dv1/dv2
    r2mm = 1-((n-1)/(n-lambda*rm)) * (SSm)
  }
  
  
  calcula.rsq  = r2mm
  
}