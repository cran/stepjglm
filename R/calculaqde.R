calculaqde <- function(y,desvio,ph,familia,dat,n1){  ##calcula o quase desvio estendido 
  
  if (familia == "gaussian") vy = seq(1,n1)
  else if (familia == "poisson") vy = y
  else if (familia == "binomial"){
    ##yu = (y+0.5)/(tot+1)
    vy = y*(1-y) ## yu*(1-yu)  ### se y=1 da erro!!!!
  }
  else if (familia == "Gamma") vy = y^2
  else if (familia == "inv gauss") vy = y^3
  
  calculaqde <- sum(desvio/ph + log(2*pi*ph*vy)  )  #Verificar isso
}