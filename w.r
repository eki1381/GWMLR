w <- function(p.temp,J,N){
  list.w <- c()
  for(j in 1:(J-1)){
    for(jprime in 1:(J-1)){
      w.temp <- matrix(0,N,N)
      if(j == jprime){
        for(i in 1:N){
          w.temp[i,i] <- p.temp[i,j]*(1-p.temp[i,j])
        }
      }else{
        for(i in 1:N){
          w.temp[i,i] <- (-1)*p.temp[i,j]*p.temp[i,jprime]
        }
      }
      list.w <- c(list.w,w.temp)
    }
  }
  ar5 <- array(list.w, dim=c(N,N*(J-1),J-1))
  w <- do.call(rbind,lapply(seq(dim(ar5)[3]), function(i) ar5[,,i]))
  return(w)
}
