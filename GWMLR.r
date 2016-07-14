GWMLR <- function(y.design.1,x.design.1,dm,b,kernel = "gaussian",dp){
  library(Matrix)
  library(sp)
  
  # y.design.1 <- model.matrix(~-1 + .,data = y)
  # x.design.1 <- model.matrix(~.,data = x)
  N <- nrow(y.design.1)
  K <- ncol(x.design.1) - 1
  J <- ncol(y.design.1)
  y.design.2 <- as.matrix(y.design.1[,2:J])
  y.design.3 <- as.vector(y.design.2)
  list <- rep(list(x.design.1),J-1)
  x.design.2 <- as.matrix(bdiag(list))
  
  #coef.matrix <- matrix(NA,(K+1)*(J-1),1)
  se.matrix <- matrix(NA,1,(K+1)*(J-1))
  tval.matrix <- matrix(NA,1,(K+1)*(J-1))
  
  weight <- kernel(dm,b,kernel = kernel)
  llike.cont <- c()
  coef.matrices <- array()
  for(i in 1:N){
    beta.1.temp <- matrix(0,K+1,J-1)
    beta.1 <- as.vector(beta.1.temp)
    beta.2.temp <- matrix(-Inf,K+1,J-1)
    beta.2 <- as.vector(beta.2.temp)
    diff.beta <- sqrt(sum((beta.1-beta.2)^2))
    
    llike.1 <- llike.weight(y.design.2,x.design.1,beta.1.temp,N,J,K,i,weight)
    llike.2 <- llike.weight(y.design.2,x.design.1,beta.2.temp,N,J,K,i,weight)
    diff.llike <- 1e9
    
    weight.rep <- rep(weight[i,],J-1)
    weight.temp <- as.matrix(diag(weight[i,]))
    
    list2 <- rep(list(weight.temp),J-1)
    weighted <- as.matrix(bdiag(list2))
    x.design.2.weighted <- weighted%*%x.design.2
    
    iterations <- 1
    while((iterations <= 20) & (diff.beta > 1e-5) & (diff.llike > 1e-5)){
      iterations <- iterations + 1
      p.temp <- prob(x.design.1,beta.1.temp)
      p <- as.vector(p.temp)
      
      w <- gww(p.temp,J,N,weight)
      
      der1.llike <- t(x.design.2.weighted)%*%(y.design.3-p)
      der2.llike <- t(x.design.2)%*%w%*%x.design.2.weighted
      
      beta.2 <- beta.1
      beta.1 <- beta.2 + (solve(der2.llike)%*%der1.llike)
      
      beta.1.temp <- matrix(beta.1,K+1,J-1)
      beta.2.temp <- matrix(beta.2,K+1,J-1)
      diff.beta <- sqrt(sum((beta.1-beta.2)^2))
      
      llike.2 <- llike.1
      llike.1 <- llike.weight(y.design.2,x.design.1,beta.1.temp,N,J,K,i,weight)
      diff.llike <- abs(llike.1 - llike.2)
    }
    llike.1 <- as.numeric(llike.1)
    llike.cont <- c(llike.cont,llike.1)
    
    varcov <- solve(der2.llike)
    se <- sqrt(diag(varcov))
    tval <- as.vector(t(beta.1.temp))/se
    coef.matrix <- matrix(as.vector(t(beta.1.temp)),(K+1)*(J-1),1)
    coef.matrix <- cbind(coef.matrix,se,tval)
    if(i == 1){
      coef.matrices <- array(coef.matrix,c((K+1)*(J-1),3,i))
    }else{
      coef.matrices <- array(c(coef.matrices,coef.matrix),c((K+1)*(J-1),3,i))
    }
    se.matrix <- rbind(se.matrix,se)
    tval.matrix <- rbind(tval.matrix,tval)
  }
  llike <- as.numeric(max(llike.cont))[1]
  aic <- as.numeric(2*(K+1)*(J-1)-2*(llike))
  aicc <- aic + ((2*(K+1)*(J-1)*((K+1)*(J-1))+1)/(N-((K+1)*(J-1))-1))
  resi.dev <- -2*llike
  # coef.matrix <- coef.matrix[-1,]
  # coef.matrix <- t(coef.matrix)
  se.matrix <- se.matrix[-1,]
  tval.matrix <- tval.matrix[-1,]
  #rownames(coef.matrix) <- c(as.character(1:N))
  colname <- c()
  sename <- c()
  tvname <- c()
  for(j in 1:(J-1)){
    for(k in 1:(K+1)){
      varname <- paste(colnames(x.design.1)[k],"=",j)
      colname <- c(colname,varname)
      # setname <- paste(colnames(x.design.1)[k],"=",j,"_SE")
      # sename <- c(sename,setname)
      # tvtname <- paste(colnames(x.design.1)[k],"=",j,"_TV")
      # tvname <- c(tvname,tvtname)
    }
  }
  rownames(coef.matrices) <- colname
  colnames(coef.matrices) <- c("Coefficients","Std Error","t-Value")
  # colnames(se.matrix) <- sename
  # colnames(tval.matrix) <- tvname
  #result <- cbind(coef.matrix,se.matrix,tval.matrix)
  
  #sumstat <- matrix(NA,(J-1)*(K+1),5)
  #for(j in 1:(J-1)){
  # for(k in 1:(K+1)){
  #    sumstat[j*k,1] <- min(coef.matrix[,j*k])
  #    sumstat[j*k,2] <- quantile(coef.matrix[,j*k],0.25)
  #    sumstat[j*k,3] <- median(coef.matrix[,j*k])
  #    sumstat[j*k,4] <- quantile(coef.matrix[,j*k],0.75)
  #    sumstat[j*k,5] <- max(coef.matrix[,j*k])
  # }
  #}
  #colnames(sumstat) <- c("Min","1st","Median","3rd","Max")
  #rownames(sumstat) <- colnames(coef.matrix)
  
  list(coefficients = coef.matrices,sumstat = sumstat,llike = llike,aic = aic,aicc = aicc,resi.dev = resi.dev,kernel = kernel,bandwidth = b,N = N,dp = dp)
}
