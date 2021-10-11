###########################################################################################
### diff.perm.ccf function ###
###########################################################################################
#@ts = vector of functional diversity values
#@timedat = vector of dates/time index
#@span = number of time points to cross correlate and/or spatial permute across
#@iter = number of permutations
#@comp.ts = vector of values to cross correlate against
#@comp.ts.timedat = vector of dates/time index for the comparative time series
#@diff = should both ts and comp.ts be differenced
#@lag = how many time points should the differencing take place across
#@pre.diff = is the comp.ts pre differenced?
#@perm.method = which pseudo-randmomisation technique should be used:
              # "spatial" (shuffling within temporal span),
              # "arima" (residuals of ARIMA model),
              # "red.noise" (red noise process using data mean, variance and autocorrelation coef) 
              # or "replacement" (sample with replacement)
#@scale = should both ts and comp.ts be centered to mean 0 and scaled to unit variance
#@normalise = should both ts and comp.ts be normalised by ranking observations by their distribution
#@identical.t = for cross correlation, should the permuted correlation coefs be taken from the same optimal 
              #lag as the observed (identical.t = T) or from the optimal lag of the permuted ts (identical.t = F)


diff.perm.ccf <- function(ts, timedat, span = 12*5, iter = 999,
                          comp.ts,comp.ts.timedat = NULL,diff = T, 
                          lag=12,pre.diff = F,
                          perm.method = c("spatial","arima","red.noise","replacement"),
                          scale = F,normalise = F,identical.t = F){
  
  require(fgpt) # spatial permutation
  require(foreach) # parallelable 'for' function
  require(doParallel) # helper functions for foreach
  require(parallel) # makeCluster function
  require(dplyr) # dplyr, ggplot etc.
  require(forecast) # ARIMA function
  require(zoo) # worker functions for timeseries 
  set.seed(123) 
  
  ## Necessary custom functions ##
  comb <- function(...) {
    mapply('rbind', ..., SIMPLIFY=FALSE)
  } #combination function for 'foreach'
  
  red.noise.ts <- function(ts,lag,length){
    x  <- rep(NA, length)
    ts <- rnorm(length,mean = mean(ts,na.rm=T),sd=sd(ts,na.rm=T))
    x[1] <- ts[1]
    for(i in seq(from = 2, to = length, by =1)){
      x[i] <- lag*x[i-1] + ((1-lag^2)^0.5)*ts[i]
    }
    return(x) #create rednoise process from a time series' mean, variance and autocorrelation coef
  }
  
  
  ## Create and register cluster for parallelisation
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  doParallel::registerDoParallel(cl)
  
  perm.meth <- match.arg(perm.method)
  
  if(any(is.na(ts))){
    print("Warning: ts has missing values. NAs have been interpolated")
    ts <- zoo::na.approx(ts,na.rm=F)
  }
  if(any(is.na(comp.ts))){
    print("Warning: comp.ts has missing values. NAs have been interpolated")
    comp.ts <- zoo::na.approx(comp.ts,na.rm=F)
  }
  
  if(isTRUE(diff)){
    if(isFALSE(pre.diff)){
      ts.diff <- diff(ts, lag = lag)
      comp.diff <- diff(comp.ts, lag = lag)
      
    }else{
      ts.diff <- diff(ts, lag = lag)
      comp.diff <- comp.ts
    }
  }else{
    ts.diff <- ts
    comp.diff <- comp.ts
  }
  ##normalise derivatives via quantiles
  if(isTRUE(normalise)){
    ts.diff.ecdf <- ecdf(ts.diff)
    ts.diff <- ts.diff.ecdf(ts.diff)
    comp.diff.ecdf <- ecdf(comp.diff)
    comp.diff <- comp.diff.ecdf(comp.diff)
  }else{
    ts.diff <- ts.diff
    comp.diff <- comp.diff
  }
  
  
  ##scale derivatives to mean zero and equal SD
  if(isTRUE(scale)){
    ts.diff <- c(scale(ts.diff))
    comp.diff <- c(scale(comp.diff))
  }else{
    ts.diff <- ts.diff
    comp.diff <- comp.diff
  }
  
  #match timedat if not 1:1 (& drop comp.ts dates X in FD data)
  if(!is.null(comp.ts.timedat)){
    #ts.diff <- na.omit(ts.diff[(comp.ts.timedat-lag)])
    ts.diff <- na.omit(ts.diff[(comp.ts.timedat)])
  }
  
  #calculate observed correlation
  ccf.tmp <- ccf(ts.diff,comp.diff,lag.max = span,plot = F)
  ccf.tmp <- data.frame("lag" = ccf.tmp$lag,"acf" = ccf.tmp$acf)
  
  
  ccf.obs <- data.frame("tmin" = ccf.tmp$lag[ccf.tmp$acf == min(ccf.tmp$acf,na.rm=TRUE)])
  ccf.obs$rmin <-  ccf.tmp$acf[ccf.tmp$lag == ccf.obs$tmin]
  ccf.obs$tmax <- ccf.tmp$lag[ccf.tmp$acf == max(ccf.tmp$acf,na.rm=TRUE)]
  ccf.obs$rmax <-  ccf.tmp$acf[ccf.tmp$lag == ccf.obs$tmax]
  ccf.obs$r0 <- ccf.tmp$acf[ccf.tmp$lag == 0]
  ccf.obs$t.absmax <- ccf.tmp$lag[abs(ccf.tmp$acf) == max(abs(ccf.tmp$acf),na.rm=TRUE)]
  ccf.obs$abs.rmax <- ccf.tmp$acf[ccf.tmp$lag == ccf.obs$t.absmax]
  
  ## spatial autocorrelation based permutation
  if(perm.meth == "spatial"){
    ## Create permutation order
    xy <- cbind(timedat,timedat) #temporal 'grid references' for spatial permutation
    
    perm.df <- fgpt::fgperm(xy,marks=ts, scale=span,  iter=iter, ratio=1, FUN=fgpt::fyshuffle, 
                            bootstrap = F, add.obs=F, as.matrix=F) #create semi-random order of indices
    
    ##loop through each permutation to cross correlate and extract optimal lag
    tmp <- foreach::foreach(i = c(1:iter),.combine = "comb",.multicombine = T, .packages = c("mgcv","gratia","dplyr")) %dopar%{
      
      if(isTRUE(diff)){
        ts.diff.perm <- diff(ts[perm.df[[i]]], lag = lag)
      }else{
        ts.diff.perm <- ts[perm.df[[i]]]
      }
      
      ##normalise derivatives via quantiles
      if(isTRUE(normalise)){
        ts.diff.perm.ecdf <- ecdf(ts.diff.perm)
        ts.diff.perm <- ts.diff.perm.ecdf(ts.diff.perm)
      }else{
        ts.diff.perm <- ts.diff.perm
      }
      ##scale derivatives to mean zero and equal SD
      if(isTRUE(scale)){
        ts.diff.perm <- c(scale(ts.diff.perm))
      }else{
        ts.diff.perm <- ts.diff.perm
      }
      
      
      #match timedat if not 1:1
      if(!is.null(comp.ts.timedat)){
        ts.diff.perm <- na.omit(ts.diff.perm[(comp.ts.timedat)])
      }
      
      #calculate observed correlation
      ccf.perm.tmp <- ccf(ts.diff.perm,comp.diff,lag.max = span,plot = F)
      ccf.perm.tmp <- data.frame("lag" = ccf.perm.tmp$lag,"acf" = ccf.perm.tmp$acf)
      
      
      if(isTRUE(identical.t)){
        ccf.perm.obs <- data.frame("tmin" =ccf.obs$tmin) #match observed tmin/tmax for comparison
        ccf.perm.obs$rmin <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmin]
        ccf.perm.obs$tmax <- ccf.obs$tmax
        ccf.perm.obs$rmax <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmax]
        ccf.perm.obs$r0 <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == 0]
        ccf.perm.obs$t.absmax <- ccf.obs$t.absmax
        ccf.perm.obs$abs.rmax <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$t.absmax]
      }else{
        ccf.perm.obs <- data.frame("tmin" = ccf.perm.tmp$lag[ccf.perm.tmp$acf == min(ccf.perm.tmp$acf,na.rm=TRUE)])
        ccf.perm.obs$rmin <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmin]
        ccf.perm.obs$tmax <- ccf.perm.tmp$lag[ccf.perm.tmp$acf == max(ccf.perm.tmp$acf,na.rm=TRUE)]
        ccf.perm.obs$rmax <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmax]
        ccf.perm.obs$r0 <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == 0]
        ccf.perm.obs$t.absmax <- ccf.perm.tmp$lag[abs(ccf.perm.tmp$acf) == max(abs(ccf.perm.tmp$acf),na.rm=TRUE)]
        ccf.perm.obs$abs.rmax <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$t.absmax]
      }
      
      return(ccf.perm.obs)
    }
  }
  
  ## ARIMA based permutation
  if(perm.meth == "arima"){
    ## Create permutation order
    fit <- forecast::auto.arima(stats::as.ts(ts),allowdrift =TRUE,seasonal = T)
    perm.df <- data.frame(date = timedat[,1],
                          pred = fitted(fit),
                          upr = fitted(fit) + 1.96*sqrt(fit$sigma2),
                          lwr = fitted(fit) - 1.96*sqrt(fit$sigma2))
    for (perm in c(1:iter)) { # for each permutation randomly sample from confidence interval at each time point
      new_col_name <- paste0("perm_", perm) # new column and colname for each permutation
      perm.df <- perm.df %>% 
        rowwise()%>%
        mutate(!!sym(new_col_name) := sample(x = seq(from = lwr, to = upr,by = 0.0001),size = 1))%>%
        ungroup()
    }
    
    ##loop through each permutation to cross correlate and extract optimal lag
    tmp <- foreach::foreach(i = 1:iter,.combine = "comb",.multicombine = T, .packages = c("dplyr")) %dopar%{     
      if(isTRUE(monthly)){ # first difference (seasonal)
        ts.diff.perm <- as.numeric(base::diff(as.ts(perm.df[,paste("perm",i,sep = "_")]), lag = 12)) #as.ts required due to finicky tibble interactions with diff
      }else{ #first difference (yearly)
        ts.diff.perm <- as.numeric(diff(as.ts(perm.df[,paste("perm",i,sep = "_")]), lag = 1))
      }
      
      ##normalise derivatives via quantiles
      if(isTRUE(normalise)){
        ts.diff.perm.ecdf <- ecdf(ts.diff.perm)
        ts.diff.perm <- ts.diff.perm.ecdf(ts.diff.perm)
      }else{
        ts.diff.perm <- ts.diff.perm
      }
      ##scale derivatives to mean zero and equal SD
      if(isTRUE(scale)){
        ts.diff.perm <- c(scale(ts.diff.perm))
      }else{
        ts.diff.perm <- ts.diff.perm
      }
      
      #calculate observed correlation
      if(monthly == T){
        ccf.perm.tmp <- ccf(ts.diff.perm,comp.diff,lag.max = 12*5,plot = F)
        ccf.perm.tmp <- data.frame("lag" = ccf.perm.tmp$lag,"acf" = ccf.perm.tmp$acf)
        
      }else{
        ccf.perm.tmp <- ccf(ts.diff.perm,comp.diff,lag.max = 5,plot = F)
        ccf.perm.tmp <- data.frame("lag" = ccf.perm.tmp$lag,"acf" = ccf.perm.tmp$acf)
      }
      
      ccf.perm.obs <- data.frame("tmin" =ccf.obs$tmin) #match observed tmin/tmax for comparison
      ccf.perm.obs$rmin <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmin]
      ccf.perm.obs$tmax <- ccf.obs$tmax
      ccf.perm.obs$rmax <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmax]
      ccf.perm.obs$r0 <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == 0]
      ccf.perm.obs$t.absmax <- ccf.obs$t.absmax
      ccf.perm.obs$abs.rmax <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$t.absmax]
      
      return(ccf.perm.obs)
    }
  }
  
  ## rednoise based permutation
  if(perm.meth == "red.noise"){
    ## Create permutation order
    # fit ar1 model to ts and extract autocorrelation coefficient
    if(span>=12){
      d1.ar1 <- stats::arima(stats::as.ts(ts), order = c(1, 0, 0),
                             seasonal = list(order = c(1, 0, 0), period = 12),
                             optim.control = list(maxit = 1000),method="ML")$coef[1]
    }else{
      d1.ar1 <- stats::arima(stats::as.ts(ts), order = c(1, 0, 0),optim.control = list(maxit = 1000),method="ML")$coef[1]
    }
    
    perm.df <- matrix(NA,nrow=length(timedat),ncol=iter)
    
    for (r in seq_len(iter)) {
      perm.df[,r] <- red.noise.ts(ts,lag=d1.ar1,length=length(timedat)) #permuted autocorrelated surrogates
    }  
    
    ##loop through each permutation to cross correlate and extract optimal lag
    tmp <- foreach::foreach(i = 1:iter,.combine = "comb",.multicombine = T, .packages = c("dplyr")) %dopar%{
      if(isTRUE(diff)){
        ts.diff.perm <- diff(perm.df[,i], lag = lag)
      }else{
        ts.diff.perm <- perm.df[,i]
      }
      ##normalise derivatives via quantiles
      if(isTRUE(normalise)){
        ts.diff.perm.ecdf <- ecdf(ts.diff.perm)
        ts.diff.perm <- ts.diff.perm.ecdf(ts.diff.perm)
      }else{
        ts.diff.perm <- ts.diff.perm
      }
      ##scale derivatives to mean zero and equal SD
      if(isTRUE(scale)){
        ts.diff.perm <- c(scale(ts.diff.perm))
      }else{
        ts.diff.perm <- ts.diff.perm
      }
      
      #match timedat if not 1:1
      if(!is.null(comp.ts.timedat)){
        ts.diff.perm <- na.omit(ts.diff.perm[(comp.ts.timedat)])
      }
      
      #calculate observed correlation
      ccf.perm.tmp <- ccf(ts.diff.perm,comp.diff,lag.max = span,plot = F)
      ccf.perm.tmp <- data.frame("lag" = ccf.perm.tmp$lag,"acf" = ccf.perm.tmp$acf)
      
      if(isTRUE(identical.t)){
        ccf.perm.obs <- data.frame("tmin" =ccf.obs$tmin) #match observed tmin/tmax for comparison
        ccf.perm.obs$rmin <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmin]
        ccf.perm.obs$tmax <- ccf.obs$tmax
        ccf.perm.obs$rmax <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmax]
        ccf.perm.obs$r0 <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == 0]
        ccf.perm.obs$t.absmax <- ccf.obs$t.absmax
        ccf.perm.obs$abs.rmax <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$t.absmax]
      }else{
        ccf.perm.obs <- data.frame("tmin" = ccf.perm.tmp$lag[ccf.perm.tmp$acf == min(ccf.perm.tmp$acf,na.rm=TRUE)])
        ccf.perm.obs$rmin <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmin]
        ccf.perm.obs$tmax <- ccf.perm.tmp$lag[ccf.perm.tmp$acf == max(ccf.perm.tmp$acf,na.rm=TRUE)]
        ccf.perm.obs$rmax <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmax]
        ccf.perm.obs$r0 <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == 0]
        ccf.perm.obs$t.absmax <- ccf.perm.tmp$lag[abs(ccf.perm.tmp$acf) == max(abs(ccf.perm.tmp$acf),na.rm=TRUE)]
        ccf.perm.obs$abs.rmax <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$t.absmax]
      }
      
      return(ccf.perm.obs)
    }
  }
  
  ## resample with replacement
  if(perm.meth == "replacement"){
    ## Create permutation order
    perm.df <- data.frame(date = timedat,
                          ts=ts)
    for (perm in c(1:iter)) { # for each permutation randomly sample from ts with replacement at each time point
      new_col_name <- paste0("perm_", perm) # new column and colname for each permutation
      perm.df <- perm.df %>% 
        mutate(!!sym(new_col_name) :=sample(x = ts,replace = TRUE,size = length(ts)))
    }
    tmp <- foreach::foreach(i = 1:iter,.combine = "comb",.multicombine = T, .packages = c("dplyr")) %dopar%{     
      
      if(isTRUE(diff)){ 
        ts.diff.perm <- diff(perm.df[,paste("perm",i,sep = "_")], lag = lag)
      }else{ 
        ts.diff.perm <- perm.df[,paste("perm",i,sep = "_")]
      }
      
      ##normalise derivatives via quantiles
      if(isTRUE(normalise)){
        ts.diff.perm.ecdf <- ecdf(ts.diff.perm)
        ts.diff.perm <- ts.diff.perm.ecdf(ts.diff.perm)
      }else{
        ts.diff.perm <- ts.diff.perm
      }
      ##scale derivatives to mean zero and equal SD
      if(isTRUE(scale)){
        ts.diff.perm <- c(scale(ts.diff.perm))
      }else{
        ts.diff.perm <- ts.diff.perm
      }
      
      
      #match timedat if not 1:1
      if(!is.null(comp.ts.timedat)){
        ts.diff.perm <- na.omit(ts.diff.perm[(comp.ts.timedat)])
      }
      
      #calculate observed correlation
      ccf.perm.tmp <- ccf(ts.diff.perm,comp.diff,lag.max = span,plot = F)
      ccf.perm.tmp <- data.frame("lag" = ccf.perm.tmp$lag,"acf" = ccf.perm.tmp$acf)
      
      if(isTRUE(identical.t)){
        ccf.perm.obs <- data.frame("tmin" =ccf.obs$tmin) #match observed tmin/tmax for comparison
        ccf.perm.obs$rmin <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmin]
        ccf.perm.obs$tmax <- ccf.obs$tmax
        ccf.perm.obs$rmax <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmax]
        ccf.perm.obs$r0 <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == 0]
        ccf.perm.obs$t.absmax <- ccf.obs$t.absmax
        ccf.perm.obs$abs.rmax <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$t.absmax]
      }else{
        ccf.perm.obs <- data.frame("tmin" = ccf.perm.tmp$lag[ccf.perm.tmp$acf == min(ccf.perm.tmp$acf,na.rm=TRUE)])
        ccf.perm.obs$rmin <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmin]
        ccf.perm.obs$tmax <- ccf.perm.tmp$lag[ccf.perm.tmp$acf == max(ccf.perm.tmp$acf,na.rm=TRUE)]
        ccf.perm.obs$rmax <-  ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$tmax]
        ccf.perm.obs$r0 <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == 0]
        ccf.perm.obs$t.absmax <- ccf.perm.tmp$lag[abs(ccf.perm.tmp$acf) == max(abs(ccf.perm.tmp$acf),na.rm=TRUE)]
        ccf.perm.obs$abs.rmax <- ccf.perm.tmp$acf[ccf.perm.tmp$lag == ccf.perm.obs$t.absmax]
      }
      
      
      return(ccf.perm.obs)
    }
  }
  parallel::stopCluster(cl)
  
  out <- matrix(NA,nrow = 7,ncol=5) #using permuted metrics, estimate permuted median and quantile of observed data
  
  out[1,] <- c("rmin.ccf",ecdf(tmp$rmin)(ccf.obs$rmin),ccf.obs$rmin,median(tmp$rmin),ccf.obs$rmin - median(tmp$rmin))
  out[2,] <- c("t.rmin.ccf",ecdf(tmp$tmin)(ccf.obs$tmin),ccf.obs$tmin,median(tmp$tmin),ccf.obs$tmin - median(tmp$tmin))
  out[3,] <- c("rmax.ccf",ecdf(tmp$rmax)(ccf.obs$rmax),ccf.obs$rmax,median(tmp$rmax),ccf.obs$rmax - median(tmp$rmax))
  out[4,] <- c("t.rmax.ccf",ecdf(tmp$tmax)(ccf.obs$tmax),ccf.obs$tmax,median(tmp$tmax),ccf.obs$tmax - median(tmp$tmax))
  out[5,] <- c("r0.ccf",ecdf(tmp$r0)(ccf.obs$r0),ccf.obs$r0,median(tmp$r0),ccf.obs$r0 - median(tmp$r0))
  out[6,] <- c("absmax.ccf",ecdf(tmp$abs.rmax)(ccf.obs$abs.rmax),ccf.obs$abs.rmax,median(tmp$abs.rmax),ccf.obs$abs.rmax - median(tmp$abs.rmax))
  out[7,] <- c("t.absmax.ccf",ecdf(tmp$t.absmax)(ccf.obs$t.absmax),ccf.obs$t.absmax,median(tmp$t.absmax),ccf.obs$t.absmax - median(tmp$t.absmax))
  
  
  
  out <- data.frame(out) %>% 
    stats::setNames(c("measure","quantile","obs.value","median.perm.value","obs.difference")) %>%
    mutate(across(quantile:obs.difference, ~as.numeric(as.character(.)))) #ensure numeric data not character
  
  out.ls <- list("perm.dens" = tmp, "summary" = out)
  return(out.ls)
}
