## SmoothFI Function ##

SmoothFI <- function(n,timedat,FI){
#Computes <FI>, a “m” point mean of FI smoothing results to focus on trends and not fluctuations.
data.smooth <- cbind(timedat, FI)
A <- data.smooth[order(data.smooth[,"V1"], decreasing = T),]
FIs <- A[,2]
timeseries <- A[,1]

window <- 0 
I <- c()
L <- c()
A <- matrix()
B <- matrix()
FIm <- matrix()

for (i in seq(from = 1, to = nrow(FI), by =  n)){ #loop through FIs

  lmin <- min(i,nrow(FIs))
  lmax <- min(i+n-1,nrow(FIs))
  NP <- nrow(FIs[lmin:lmax,1])
  
if (NP == n){
  window=window+1
  
I[window] <- mean(FIs[lmin:lmax])
L[window] <- ceiling(mean(timeseries[lmin:lmax]))
}
}

B <- cbind(t(L), t(I)) 
FIm <- B[order(B[,"V1"]),]
return(FIm)
}


#if plot==1
#%<FI>
 # subplot(2,1,1),plot(t,FI,’-*’); %plot original FI title(‘Original FI’); axis([t(1) t(end) 0 8]); subplot(2,1,2),plot(FIm(:,1),FIm(:,2),’-*’); %plot <I> title(‘Smoothed FI(< FI >)’); axis([t(1) t(end) 0 8]); end
#%%