kernel <- function(dm,b,kernel = "gaussian",i){
  n <- nrow(dm)
  w <- matrix(0,n,n)
  if(kernel == "gaussian"){
    for(i.prime in 1:n){
      w[i.prime,i.prime] <- exp(-1*(dm[i.prime,i]/b)^2)
    }
  }else if(kernel == "exponential"){
    for(i.prime in 1:n){
      w[i.prime,i.prime] <- exp(-1*(dm[i.prime,i]/b))
    }
  }
  return(w)
}
