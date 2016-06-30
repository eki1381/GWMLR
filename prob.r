prob <- function(x.design.1,beta.1.temp){
  nom <- matrix(NA,nrow(x.design.1),ncol(beta.1.temp))
  denom <- matrix(NA,nrow(x.design.1),1)
  for(i in 1:nrow(x.design.1)){
    denom.temp.1 <- 0
    for(j in 1:ncol(beta.1.temp)){
      denom.temp.2 <- 0
      nom.temp <- 0
      for(k in 1:ncol(x.design.1)){
        denom.temp.2 <- denom.temp.2 + (x.design.1[i,k]*beta.1.temp[k,j])
        nom.temp <- nom.temp + (x.design.1[i,k]*beta.1.temp[k,j])
      }
      denom.temp.2.exp <- exp(denom.temp.2)
      denom.temp.1 <- denom.temp.1 + denom.temp.2.exp
      nom[i,j] <- exp(nom.temp)
    }
    denom[i,1] <- 1 + denom.temp.1
  }
  
  p.temp <- matrix(NA,nrow(x.design.1),ncol(beta.1.temp))
  for(i in 1:nrow(x.design.1)){
    for(j in 1:ncol(beta.1.temp)){
      p.temp[i,j] <- nom[i,j]/denom[i,1]
    }
  }
  
  return(p.temp)
}
