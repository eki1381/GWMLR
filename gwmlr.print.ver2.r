gwmlr <- function(x, ...)UseMethod("gwmlr")

gwmlr.default <- function(y,x,kernel,b,dp, ...){
  require("Matrix")
  require("sp")
  dm <- dist.m(dp)
  est <- list()
  est$local <- GWMLR(y = y,x = x,dm = dm,b = b,kernel = kernel,dp = dp)
  est$global <- multinomial(y = y,x = x)
  est$call <- match.call()
  class(est) <- "gwmlr"
  est
}

print.gwmlr <- function(x, ...){
  cat("***********************************************************************\n")
  cat("**************Geographically Weighted Logistic Regression**************\n")
  cat("**************              ASGARD R Package             **************\n")
  cat("**************                    by                     **************\n")
  cat("**************           Yohanes Eki Apriliawan          **************\n")
  cat("***********************************************************************\n")
  cat("\n")
  cat("(Dependent Variable class reference has been reset to zero)")
  cat("\n")
  cat("Model\t\t=\t Logistic Regression")
  cat("\n")
  cat("Kernel\t\t=\t",x$local$kernel)
  cat("\n")
  cat("Bandwidth\t=\t",x$local$bandwidth)
  cat("\n")
  cat("Log-Likelihood\t=\t",x$local$llike)
  cat("\n")
  cat("AIC\t\t=\t",x$local$aic)
  cat("\n")
  cat("Residual Dev.\t=\t",x$local$resi.dev)
  cat("\n")
  cat("***********************************************************************\n")
  cat("\n")
  for(i in 1:x$local$N){
    cat("Obs = ",i,"\tX-Coordinate =",as.numeric(x$local$dp[i,1]),"\tY-Coordinate =",as.numeric(x$local$dp[i,2]))
    cat("\n\n")
    print(x$local$coefficients[,,i])
    cat("\n")
    cat("***********************************************************************\n")
  }
}

gwmlr.formula <- function(formula, data=list(),kernel = "gaussian",b, ...)
{
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- as.data.frame(model.response(mf))
  y.design.1 <- model.matrix(~-1 + .,data = y)
  coord <- coordinates(data)
  est <- gwmlr.default(y = y.design.1,x = x,kernel = kernel,b = b,dp = coord, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}
