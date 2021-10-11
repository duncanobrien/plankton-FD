
CalcMeanFish <- function(timedat,data,hwin,winspace,sost){
#This code was developed to calculate the mean value of the Fisher 
#information by computing FI over various tightening levels to find the 
# point above which perfect order is experienced. And then calculating a 
#mean of the values for each time window.
#AvgFish: Average Fisher information over TL from strict to relaxed
#PO_TL: Tightening Level the system is in perfect order. It is connected to pervasiveness (1-PO_TL)
#mY: Mean Year for each time window
TL <- NULL 
FisherMat <- list()#?
FMat <- FisherMat
for (i in 1:100){
    Tol <- 100-(i-1)
    
    fishertemp <- GFisher(timedat,data,sost,hwin,winspace,Tol)
    FMat[[i]] <- fishertemp$FI
    n <- length(fishertemp$FI) # num of Fisher information Calculations

#Calculates average FI for each time window by seeking the tightening level at which all 
#resulting fisher calculation for each widow is at its %max = 8, perfect order. 
#The average FI is calculated from the FisherMat which stores all of the Fisher results from 
#each time window less the windows where all of the Fisher results are 8.
      if (sum(FMat[[i]])<8*n){ #Tightening level with all Fishers = 8, perfect order % display(‘Yes’);
      FisherMat[[i]] <- FMat[[i]]
      TL <- Tol
      } #Final value is lowest tightening level before all Fisher results are 8 (max FI)
        }

if (isempty(TL)==1){
display("The system is in perfect order for this size of states at all tightening levels. Please feel free to make an adjustment accordingly!!")
PO_TL <- 100} else{
  PO_TL <- TL-1}

mY <-fishertemp$midt_win

#Provides an output when the system is completely orderly at all tightening levels. 
#This results in no collection of Fisher values before the system is a maximum level. 
#Accordingly, FisherMat is empty indicates either true perfect order or a need to adjust 
#the size of states (smaller)
if (isempty(FisherMat)==1){
AvgFish <- 8*matrix(1, nrow= length(fishertemp$FI), ncol = length(fishertemp$FI))} else{
  #Calculates mean FI for each time window
  AvgFish <- rowMeans(sapply(FisherMat, unlist))
}

#Provides an ouput when the system is completely orderly at all tighten levels. 
#Indicates either true perfect order or a need to adjust the size of states (smaller)

#Transposing (aesthetic)
AvgFish <- t(AvgFish)
mY <- t(fishertemp$midt_win)
list.out <- list("AvgFish" = AvgFish,"mY" = mY,"Min_TL" = 100-TL)
return(list.out)
}


#Simple Plots #
  #subplot(2,1,1), plot(timedat,data) %plots data
#%subplot(2,1,1), plot(data(:,1),data(:,2)) %plots data for ode title(‘Time series data’);
#subplot(2,1,2), plot(mY,AvgFish)
#axis([mY(1) mY(end) 0 9])
#title(‘Fisher information over time’); xlabel(‘Time’); ylabel(‘Fisher information’);