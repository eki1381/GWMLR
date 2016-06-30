llike <- function(y.design.2,x.design.1,beta.1.temp,N,J,K){
  res <- 0
  for(i in 1:N){
    eq.1.2 <- 0
    eq.2.1 <- 0
    for(j in 1:(J-1)){
      eq.1.1 <- 0
      for(k in 1:(K+1)){
        eq.1.1 <- eq.1.1 + (x.design.1[i,k]*beta.1.temp[k,j])
      }
      eq.1.2 <- eq.1.2 + (y.design.2[i,j]*eq.1.1)
      eq.2.1 <- eq.2.1 + exp(eq.1.1)
    }
    res <- res + (eq.1.2 - log(1+eq.2.1))
  }
  return(res)
}
