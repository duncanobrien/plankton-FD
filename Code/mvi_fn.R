### Multivariate Variance Index function ###

# replicated from Brock, W. A., and S. R. Carpenter. 2006. 
# Variance as a leading indicator of regime shift in ecosystem services. 
# Ecology and Society 11(2): 9

#@df = matrix of species abundances, names across columns, time across rows
#@window = number of time points in rolling window

multi.var.index <- function(df,window){
  
  out <- matrix(ncol = 2,nrow=dim(df)[1]-window+1)
  for(i in 1:(dim(df)[1]-window+1)){
    #tmp <- cov(df[i:(i+1),])
    tmp <- cov(df[i:(i+window-1),])
    out[i,1] <- sqrt(eigen(tmp,symmetric = T)$values[1])
    out[i,2] <- (i+window-1)
    
  }
  colnames(out) <- c("mvi","maxt")
  return(out)
}
