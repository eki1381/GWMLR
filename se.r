se <- function(varcov){
  se <- rep(NA,ncol(varcov))
  for(i in 1:nrow(varcov)){
    se[i] = sqrt(varcov[i,i])
  }
  return(se)
}
