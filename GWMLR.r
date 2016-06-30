GWMLR <- function(y,x,dm,b,kernel = "gaussian"){
  library(Matrix)
  library(MASS)
  
  y.design.1 <- model.matrix(~-1 + .,data = y)
  x.design.1 <- model.matrix(~.,data = x)
  N <- nrow(y.design.1)
  K <- ncol(x.design.1) - 1
  J <- ncol(y.design.1)
  y.design.2 <- as.matrix(y.design.1[,2:J])
  y.design.3 <- as.vector(y.design.2)
  list <- rep(list(x.design.1),J-1)
  x.design.2 <- as.matrix(bdiag(list))
  
  coef.matrix <- matrix(NA,1,(K+1)*(J-1))
  se.matrix <- matrix(NA,1,(K+1)*(J-1))
  
  for(i in 1:N){
    beta.1.temp <- matrix(0,K+1,J-1)
    beta.1 <- as.vector(beta.1.temp)
    beta.2.temp <- matrix(-Inf,K+1,J-1)
    beta.2 <- as.vector(beta.2.temp)
    diff.beta <- sqrt(sum((beta.1-beta.2)^2))
    
    llike.1 <- llike(y.design.2,x.design.1,beta.1.temp,N,J,K)
    llike.2 <- llike(y.design.2,x.design.1,beta.2.temp,N,J,K)
    diff.llike <- 1e9
    
    weight <- kernel(dm,b,i,kernel = kernel)
    
    iterations <- 1
    while((iterations <= 50) & (diff.beta > 1e-6) & (diff.llike > 1e-7)){
      iterations <- iterations + 1
      p.temp <- prob(x.design.1,beta.1.temp)
      p <- as.vector(p.temp)
      
      w <- gww(p.temp,J,N,weight)
      
      weight.rep <- rep(diag(weight),J-1)
      
      der1.llike <- t(x.design.2)%*%(weight.rep*(y.design.3-p))
      der2.llike <- t(x.design.2)%*%w%*%x.design.2
      
      beta.2 <- beta.1
      beta.1 <- beta.2 +(chol2inv(chol(der2.llike))%*%der1.llike)
      
      beta.1.temp <- matrix(beta.1,K+1,J-1)
      beta.2.temp <- matrix(beta.2,K+1,J-1)
      diff.beta <- sqrt(sum((beta.1-beta.2)^2))
      
      llike.2 <- llike.1
      llike.1 <- llike(y.design.2,x.design.1,beta.1.temp,N,J,K)
      diff.llike <- abs(llike.1 - llike.2)
    }
    llike.1 <- as.numeric(llike.1)
    aic <- as.numeric(2*(K+1)*(J-1)-2*(llike.1))
    se <- se(solve(der2.llike))

    coef.matrix <- rbind(coef.matrix,t(beta.1.temp))
    se.matrix <- rbind(se.matrix,se)
  }
  res <- cbind(coef.matrix[-1,],se.matrix[-1,])
  return(res)
}
