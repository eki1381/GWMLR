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
  cat("*************************************************************\n")
  cat("*********Geographically Weighted Logistic Regression*********\n")
  cat("*********              ASGARD R Package             *********\n")
  cat("*********                    by                     *********\n")
  cat("*********           Yohanes Eki Apriliawan          *********\n")
  cat("*************************************************************\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("(Dependent Variable class reference has been reset to zero)")
  cat("\n")
  cat
  cat("\nGlobal Model Coefficients:\n")
  print(x$global$coefficients)
  cat("\n")
  cat("AIC:",x$global$aic)
  cat("\n")
  cat("AICc:",x$global$aicc)
  cat("\n")
  cat("Log Likelihood:",x$global$likelihood)
  cat("\n")
  cat("Residual Deviance:",x$global$resi.dev)
  cat("\n")
  cat("\nLocal Model Summary Statistics:\n")
  print(x$local$sumstat)
  cat("\nAIC:")
  print(x$local$aic)
  cat("AICc:")
  print(x$local$aicc)
  cat("Log likelihood:")
  print(x$local$llike)
  cat("Residual Deviance:")
  print(x$local$resi.dev)
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
