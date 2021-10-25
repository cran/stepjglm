#' @title Variable selection in joint modeling of mean and dispersion
#' @name stepjglm
#'
#' @description A Procedure for selecting variables in JMMD (including mixture models) based on hypothesis testing and the quality of the model's fit.
#'
#' @param model an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. if \code{datafram} is a mixture data, \code{datafram} doesn't contain the principal mixture components.
#' @param alpha1 significance level for testing add new terms on the mean models.
#' @param alpha2 significance level for testing add new terms on the dispersion models.
#' @param datafram a data frame containing the data.
#' @param family a character string naming a family function or the result of a call to a family function. For \code{glm.fit} only the third option is supported. (See \code{family} for details of family functions). Describe the family function for the mean model (families implemented by package \code{stats}). For the dispersion model, the Gamma family whit log link is assumed.
#' @param lambda1 some function of the sample size to calculate the \eqn{\tilde{R}_m^{2}}  (See Pinto and Pereira (in press) and Zhang (2017) for more details). If equal to 1 (default), uses the standard correction for the \eqn{\tilde{R}_m^{2}}. If equal to "EAIC", uses the \eqn{EAIC} criterion.
#' @param lambda2 some function of the sample size to calculate the \eqn{\tilde{R}_d^{2}}  (See Pinto and Pereira (in press) and Zhang (2017) for more details). If equal to 1 (default), uses the standard correction for the \eqn{\tilde{R}_d^{2}}. If equal to "AIC", uses the corrected \eqn{AIC_c} criterion.
#' @param startmod if \code{datafram} is a mixture data, \code{startmod} is the principal mixture components, else, \code{startmod} must be equal to 1 (default).
#' @param interations if \code{TRUE} shows the outputs of iterations procedure step by step. The default is \code{FALSE}.
#'
#' @details  The function implements a method for selection of variables for both the mean and dispersion models in the JMMD introduced by Nelder and Lee (1991)
#' considering the \emph{Adjusted Quasi Extended Likelihood} introduced by Lee and Nelder (1998).
#' The method is a procedure for  selecting  variables,  based  on  hypothesis  testing  and  the  quality  of  the  model's  fit.
#' A criterion for checking the goodness of fit is used, in each iteration of the selection process,
#' as a filter for choosing the terms that will be evaluated by a hypothesis test. For more details on selection algorithms, see Pinto and Pereira (in press).
#'
#'
#' @return  \tabular{ll}{
#'    \code{ model.mean}  \tab   a \code{glm} object with the adjustments for the mean model. \cr
#'    \tab \cr
#'    \code{model.disp}  \tab   a \code{glm} object with the adjustments for the dispersion model. \cr
#'    \tab \cr
#'    \code{EAIC} \tab a numeric object containing the \emph{Extended Akaike Information Criterion}.  \cr
#'    \tab  For details, see Wang and Zhang (2009). \cr
#'    \tab \cr
#'    \code{EQD}  \tab a numeric object containing the \emph{Extended Quasi Deviance}. \cr
#'    \tab  For details, see Nelder and Lee (1991). \cr
#'    \tab \cr
#'    \code{R2m}  \tab a numeric object containing the standard correction for the \eqn{\tilde{R}_m^{2}}.  \cr
#'    \tab For details, see Pinto and Pereira (in press). \cr
#'    \tab \cr
#'    \code{R2d}  \tab a numeric object containing the standard correction for the \eqn{\tilde{R}_d^{2}}. \cr
#'    \tab For details, see Pinto and Pereira (in press). \cr
#' }
#' @author Leandro Alves Pereira, Edmilson Rodrigues Pinto.
#'
#' @seealso \code{\link[stats]{glm}}
#' @seealso \code{\link[stats]{summary.glm}}
#'
#' @import rsq
#' @import stats
#'
#' @usage  stepjglm(model,alpha1,alpha2,datafram,family,lambda1=1,lambda2=1,startmod=1,
#'                  interations=FALSE)
#' @examples
#'
#' # Application to the bread-making problem:
#'
#' data(bread_mixture)
#'
#' Form =
#' as.formula(y~ x1:x2+x1:x3+x2:x3+x1:x2:(x1-x2)+x1:x3:(x1-x3)+
#'             + x1:z1+x2:z1+x3:z1+x1:x2:z1
#'             + x1:x3:z1+x1:x2:(x1-x2):z1
#'             + x1:x3:(x1-x3):z1
#'             + x1:z2+x2:z2+x3:z2+x1:x2:z2
#'             + x1:x3:z2+x1:x2:(x1-x2):z2
#'             +x1:x3:(x1-x3):z2)
#'
#' object=stepjglm(Form,0.1,0.1,bread_mixture,gaussian,sqrt(90),"AIC","-1+x1+x2+x3")
#'
#' summary(object$modelo.mean)
#' summary(object$modelo.disp)
#'
#' object$EAIC  # Print the EAIC for the final model
#'
#'
#'
#' # Application to the injection molding data:
#'
#' form = as.formula(Y ~ A*M+A*N+A*O+B*M+B*N+B*O+C*M+C*N+C*O+D*M+D*N+D*O+
#'                       E*M+E*N+E*O+F*M+F*N+F*O+G*M+G*N+G*O)
#'
#' data(injection_molding)
#'
#' obj.dt = stepjglm(form, 0.05,0.05,injection_molding,gaussian,sqrt(nrow(injection_molding)),"AIC")
#'
#' summary(obj.dt$modelo.mean)
#' summary(obj.dt$modelo.disp)
#'
#' obj.dt$EAIC  # Print the EAIC for the final model
#' obj.dt$EQD   # Print the EQD for the final model
#' obj.dt$R2m   # Print the R2m for the final model
#' obj.dt$R2d   # Print the R2d for the final model
#'
#' @references Hu, B. and Shao, J. (2008). Generalized linear model selection using \eqn{R^2}. \emph{Journal of Statistical Planning and Inference}, 138, 3705-3712.
#'
#' @references Lee, Y., Nelder, J. A. (1998). Generalized linear models for analysis of quality improvement experiments. \emph{The Canadian Journal of Statistics}, v. 26, n. 1, pp. 95-105.
#'
#' @references Nelder, J. A., Lee, Y. (1991). Generalized linear models for the analysis of Taguchi-type experiments. \emph{Applied Stochastic Models and Data Analysis}, v. 7, pp. 107-120.
#'
#' @references Pinto, E. R., Pereira, L. A. (in press). On variable selection in joint modeling of mean and dispersion. \emph{Brazilian Journal of Probability and Statistics}. Preprint at   \url{https://arxiv.org/abs/2109.07978} (2021).
#'
#' @references Wang, D. and Zhang, Z. (2009). Variable selection in joint generalized linear models. \emph{Chinese Journal of Applied Probability and Statistics},  v. 25,  pp.245-256.
#'
#' @references Zhang, D. (2017). A coefficient of determination for generalized linear models. \emph{The American Statistician}, v. 71, 310-316.
#'
#' @export
stepjglm <- function(model,alpha1,alpha2,datafram,family,lambda1=1,lambda2=1,startmod=1,interations=FALSE){

  modelo = as.formula(model)
  dad = datafram
  familia = family
  modini = startmod
  saidas = interations

  n = nrow(dad)
  variaveis <- model.frame(paste("~",modelo[3]),data=dad) #variaveis explicativas
  resposta  <- model.frame(paste("~",modelo[2]),data=dad) # variavel resposta

  objetod = list(modelo = "constante")
  if (lambda1 == "EAIC") objetom = stepmean2(modelo,dad,resposta,1,alpha1,familia,lambda1,modini,"constante",saidas)
    else  objetom = stepmean(modelo,dad,resposta,1,alpha1,familia,lambda1,modini,saidas)

  if (saidas) print(objetom$R)
  if (saidas) print("---------------------------------------------------------------------------------------------------")
  R2 = sum(ifelse(objetom$R[,7]=="X",objetom$R[,2],0)) #Obtem o R2 do modelo final
  repeat{
    if (lambda2 == "AIC")  objetodaux = stepdisp(modelo,dad,objetom$d,alpha2,modini,saidas)
     else objetodaux = stepdisp2(modelo,dad,objetom$d,alpha2,modini,lambda2,saidas)
     if (saidas) print(objetodaux$R)
     if (saidas) print("-------------------------------------------------------------------------------------------------")
    if (lambda1 == "EAIC") objetomaux = stepmean2(modelo,dad,resposta,objetodaux$phi,alpha1,familia,lambda1,modini,objetodaux$modelo,saidas)
    else  objetomaux = stepmean(modelo,dad,resposta,objetodaux$phi,alpha1,familia,lambda1,modini,saidas)
    if (saidas) print(objetomaux$R)
    if (saidas) print("-------------------------------------------------------------------------------------------------")
    R2.novo = sum(ifelse(objetomaux$R[,7]=="X",objetomaux$R[,2],0)) #Obtem o R2 do modelo final

    if  (lambda1 == "EAIC"){
    if (R2.novo < R2 ){
      R2 = R2.novo    #  Se o R2 aumentou, entao continua
      objetom = objetomaux
      objetod = objetodaux
    }
    else break
    }# Se o R2 caiu, para o loop
    else{
      if (R2.novo > R2 ){
        R2 = R2.novo    #  Se o R2 aumentou, entao continua
        objetom = objetomaux
        objetod = objetodaux
      }
      else break
    }

 }
  modelo.mean = objetom$modelo
  modelo.disp = objetod$modelo
  k = (length(modelo.mean$coef)+ length(modelo.disp$coef))
  d1 = (resid( modelo.mean,type="deviance")^2) /(1-hatvalues( modelo.mean))
  phi = fitted.values(modelo.disp)

  R2m = calcula.rsq(modelo.mean,ob1,phi,1,modelo.mean$family[1]$family)

  d = resid(modelo.mean,type="deviance")^2 /(1-hatvalues(modelo.mean))
  ob1d = glm("d~ 1",family = Gamma (link=log),data=cbind(dad,d))
  R2d = calcula.rsq(modelo.disp,ob1d,rep(1,n),1,"Gamma")

  EAIC = calculaqde(modelo.mean$y,d1,phi,modelo.mean$family[1]$family,dad,n) + (2*k*n)/(n-k-1)
  EQD = calculaqde(modelo.mean$y,d1,phi,modelo.mean$family[1]$family,dad,n)

  stepjglm = list(modelo.mean = objetom$modelo, modelo.disp = objetod$modelo, EAIC = EAIC, EQD= EQD, R2m = R2m, R2d = R2d)


}

