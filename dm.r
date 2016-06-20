#Euclidean Distance Function
gwr.dm <- function(dp)
{
  #Step 1 : Count the length of data points
  length <- length(dp)/2
  
  #Step 2 : Make an empty distance matrix
  dm <- matrix(,nrow = length,ncol = length)
  
  #Step 3 : Count each distance of data points using Euclidean Distance and loop, then fill distance matrix
  for(i in 1:length)
  {
    for(j in 1:length){
      dm[i,j] <- sqrt((dp[i]-dp[j])^2+(dp[length+i]-dp[length+j])^2)
    }
  }
  
  #Step 4 : Return the distance matrix
  dm
}
