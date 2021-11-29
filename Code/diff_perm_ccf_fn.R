###########################################################################################
### diff.perm.ccf function ###
###########################################################################################
#@dat = matrix of data inputs. dat[,1] must be a numeric date vector,
# dat[,2] is the hypothesised causal variable,
# dat[,3] is the response variable of interest.
#@iter = number of permutations
#@span = number of lag time points to cross map across
#@return.raw = whether each cross skill at each lag should be returned
#@detrend.method = "none", "diff" (differenced), "lm" (standardised residuals of a lm between causal variable and time)
#@lag = if 'diff' detrending method selected, how many lags to difference over
#@perm.method = which pseudo-randmomisation technique should be used:
# "spatial" (shuffling of observed data within temporal span),
# "arima" (sampled from predictions of an ARIMA model),
# "red.noise" (red noise process using data mean, variance and autocorrelation coef) 
# or "replacement" (sampled from observed data with replacement)
#@identical.t = for cross correlation, should the permuted correlation coefs be taken from the same optimal 
#lag as the observed (identical.t = T) or from the optimal lag of the permuted ts (identical.t = F)


diff.perm.ccf <- function(dat, span = 12*5, iter = 999,
                          lag=0, detrend.method = c("none","diff","lm"),
                          perm.method = c("spatial","arima","red.noise","replacement"),
                          identical.t = F){
  
  require(fgpt) # spatial permutation
  require(foreach) # parallelable 'for' function
  require(doParallel) # helper functions for foreach
  require(parallel) # makeCluster function
  require(dplyr) # dplyr, ggplot etc.
  require(forecast) # auto.arima function
  require(zoo) # worker functions for timeseries 
  set.seed(123) 
  
  ## Necessary custom functions ##
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
  
  if(any(is.na(dat[,2]))){
    warning("predictor variable has missing values. NAs have been interpolated")
    dat[,2] <- zoo::na.approx(dat[,2],na.rm=F)
  }
  if(any(is.na(dat[,3]))){
    warning("response variable has missing values. NAs have been interpolated")
    dat[,3] <- zoo::na.approx(dat[,3],na.rm=F)
  }
  
  if(detrend.method == "diff"){
    sub.dat <- dat[(1+lag):nrow(dat),]
    sub.dat[,2] <- base::diff(dat[,2], lag = lag)
    sub.dat[,3] <- base::diff(dat[,3], lag = lag)
  }else if(detrend.method == "lm"){
    sub.dat <- dat
    sub.dat[,2] <- residuals(lm(sub.dat[,2] ~ as.numeric(sub.dat[,1]),na.action=na.exclude))- 
      c(decompose(ts(sub.dat[,2],frequency = 12),type = "additive")$seasonal)
    sub.dat[,3] <- residuals(lm(sub.dat[,3] ~ as.numeric(sub.dat[,1]),na.action=na.exclude))- 
      c(decompose(ts(sub.dat[,3],frequency = 12),type = "additive")$seasonal)
  }else{
    sub.dat <- dat
  }
  sub.dat <- na.omit(sub.dat) %>% mutate(across(everything(),~as.numeric(.))) %>% #drop lead/lag NAs and ensure all columns are numeric
    mutate(across(everything(),~scale(.)))
  
  #calculate observed correlation
  ccf.tmp <- ccf(c(sub.dat[,2]),c(sub.dat[,3]),lag.max = span,plot = F)
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
    xy <- cbind(dat[,1],dat[,1]) #temporal 'grid references' for spatial permutation
    
    perm.df <- fgpt::fgperm(xy,marks=dat[,2], scale=span,  iter=iter, ratio=1, FUN=fgpt::fyshuffle, 
                            bootstrap = F, add.obs=F, as.matrix=F) #create semi-random order of indices
    
    ##loop through each permutation to cross correlate and extract optimal lag
    tmp <- foreach::foreach(i = c(1:iter),.combine = "rbind",.multicombine = T, .packages = c("dplyr")) %dopar%{
      
      if(detrend.method == "diff"){
        perm.dat <- dat[(1+lag):nrow(dat),]
        perm.dat[,2] <- base::diff(dat[perm.df[[i]],2], lag = lag)
        perm.dat[,3] <- base::diff(dat[,3], lag = lag)
      }else if(detrend.method == "lm"){
        perm.dat <- dat
        perm.dat[,2] <- residuals(lm(perm.dat[perm.df[[i]],2] ~ as.numeric(perm.dat[,1]),na.action=na.exclude)) - 
          c(decompose(ts(perm.dat[perm.df[[i]],2],frequency = 12),type = "additive")$seasonal) #extract linear model residuals and subtract seasonal component
        perm.dat[,3] <- residuals(lm(perm.dat[,3] ~ as.numeric(perm.dat[,1]),na.action=na.exclude)) - 
          c(decompose(ts(perm.dat[,3],frequency = 12),type = "additive")$seasonal) #extract linear model residuals and subtract seasonal component
      }else{
        perm.dat <- dat
        perm.dat[,2] <- dat[perm.df[[i]],2]
      }
      perm.dat <- na.omit(perm.dat) %>% mutate(across(everything(),~as.numeric(.))) %>% #drop lead/lag NAs and ensure all columns are numeric
        mutate(across(everything(),~scale(.)))
      
      #calculate observed correlation
      ccf.perm.tmp <- ccf(c(perm.dat[,2]),c(perm.dat[,3]),lag.max = span,plot = F)
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
    fit <- forecast::auto.arima(stats::as.ts(dat[,2]),allowdrift =TRUE,seasonal = T)
    perm.df <- data.frame("date" = as.numeric(dat[,1]),
                          "pred" = fitted(fit),
                          "upr" = fitted(fit) + 1.96*sqrt(fit$sigma2),
                          "lwr" = fitted(fit) - 1.96*sqrt(fit$sigma2)) %>%
      stats::setNames(c("date","pred","upr","lwr"))
    
    for (perm in c(1:iter)) { # for each permutation randomly sample from confidence interval at each time point
      new_col_name <- paste0("perm_", perm) # new column and colname for each permutation
      perm.df <- perm.df %>% 
        rowwise()%>%
        mutate(!!sym(new_col_name) := sample(x = seq(from = lwr, to = upr,by = 0.0001),size = 1))%>%
        ungroup()%>%
        as.data.frame()
    }
    
    ##loop through each permutation to cross correlate and extract optimal lag
    tmp <- foreach::foreach(i = 1:iter,.combine = "rbind",.multicombine = T, .packages = c("dplyr")) %dopar%{     
      
      if(detrend.method == "diff"){
        perm.dat <- dat[(1+lag):nrow(dat),]
        perm.dat[,2] <- base::diff(perm.df[,paste("perm",i,sep = "_")], lag = lag)
        perm.dat[,3] <- base::diff(dat[,3], lag = lag)
      }else if(detrend.method == "lm"){
        perm.dat <- dat
        perm.dat[,2] <- residuals(lm(perm.df[,paste("perm",i,sep = "_")] ~ as.numeric(perm.dat[,1]),na.action=na.exclude)) - 
          c(decompose(ts(perm.df[,paste("perm",i,sep = "_")],frequency = 12),type = "additive")$seasonal)
        perm.dat[,3] <- residuals(lm(perm.dat[,3] ~ as.numeric(perm.dat[,1]),na.action=na.exclude)) - 
          c(decompose(ts(perm.dat[,3],frequency = 12),type = "additive")$seasonal)
      }else{
        perm.dat <- dat
        perm.dat[,2] <- perm.df[,paste("perm",i,sep = "_")]
      }
      perm.dat <- na.omit(perm.dat) %>% mutate(across(everything(),~as.numeric(.))) %>% #drop lead/lag NAs and ensure all columns are numeric
            mutate(across(everything(),~scale(.)))
      
      #calculate observed correlation
      ccf.perm.tmp <- ccf(c(perm.dat[,2]),c(perm.dat[,3]),lag.max = span,plot = F)
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
  
  ## rednoise based permutation
  if(perm.meth == "red.noise"){
    ## Create permutation order
    # fit ar1 model to ts and extract autocorrelation coefficient
    if(span>=12){
      d1.ar1 <- stats::arima(stats::as.ts(dat[,2]), order = c(1, 0, 0),
                             seasonal = list(order = c(1, 0, 0), period = 12),
                             optim.control = list(maxit = 1000),method="ML")$coef[1]
    }else{
      d1.ar1 <- stats::arima(stats::as.ts(dat[,2]), order = c(1, 0, 0),optim.control = list(maxit = 1000),method="ML")$coef[1]
    }
    
    perm.df <- matrix(NA,nrow=length(dat[,3]),ncol=iter)
    
    for (r in seq_len(iter)) {
      perm.df[,r] <- red.noise.ts(dat[,2],lag=d1.ar1,length=nrow(dat)) #permuted autocorrelated surrogates
    }  
    
    ##loop through each permutation to cross correlate and extract optimal lag
    tmp <- foreach::foreach(i = 1:iter,.combine = "rbind",.multicombine = T, .packages = c("dplyr")) %dopar%{
      
      if(detrend.method == "diff"){
        perm.dat <- dat[(1+lag):nrow(dat),]
        perm.dat[,2] <- base::diff(perm.df[,i], lag = lag)
        perm.dat[,3] <- base::diff(dat[,3], lag = lag)
      }else if(detrend.method == "lm"){
        perm.dat <- dat
        perm.dat[,2] <- residuals(lm(perm.df[,i] ~ as.numeric(perm.dat[,1]),na.action=na.exclude)) - 
          c(decompose(ts(perm.df[,i],frequency = 12),type = "additive")$seasonal)
        perm.dat[,3] <- residuals(lm(perm.dat[,3] ~ as.numeric(perm.dat[,1]),na.action=na.exclude)) - 
          c(decompose(ts(perm.dat[,3],frequency = 12),type = "additive")$seasonal)
      }else{
        perm.dat <- dat
        perm.dat[,2] <- perm.df[,i]
      }
      perm.dat <- na.omit(perm.dat) %>% mutate(across(everything(),~as.numeric(.))) %>% #drop lead/lag NAs and ensure all columns are numeric
        mutate(across(everything(),~scale(.)))
      
      #calculate observed correlation
      ccf.perm.tmp <- ccf(c(perm.dat[,2]),c(perm.dat[,3]),lag.max = span,plot = F)
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
    perm.df <- data.frame(date = dat[,1],
                          ts=dat[,2])
    for (perm in c(1:iter)) { # for each permutation randomly sample from ts with replacement at each time point
      new_col_name <- paste0("perm_", perm) # new column and colname for each permutation
      perm.df <- perm.df %>% 
        mutate(!!sym(new_col_name) :=sample(x = dat[,2],replace = TRUE,size = nrow(dat)))
    }
    
    tmp <- foreach::foreach(i = 1:iter,.combine = "rbind",.multicombine = T, .packages = c("dplyr")) %dopar%{     
      
      if(detrend.method == "diff"){
        perm.dat <- dat[(1+lag):nrow(dat),]
        perm.dat[,2] <- base::diff(perm.df[,paste("perm",i,sep = "_")], lag = lag)
        perm.dat[,3] <- base::diff(dat[,3], lag = lag)
      }else if(detrend.method == "lm"){
        perm.dat <- dat
        perm.dat[,2] <- residuals(lm(perm.df[,paste("perm",i,sep = "_")] ~ as.numeric(perm.dat[,1]),na.action=na.exclude))- 
          c(decompose(ts(perm.df[,paste("perm",i,sep = "_")],frequency = 12),type = "additive")$seasonal)
        perm.dat[,3] <- residuals(lm(perm.dat[,3] ~ as.numeric(perm.dat[,1]),na.action=na.exclude))- 
          c(decompose(ts(perm.dat[,3],frequency = 12),type = "additive")$seasonal)
      }else{
        perm.dat <- dat
        perm.dat[,2] <- perm.df[,paste("perm",i,sep = "_")]
      }
      perm.dat <- na.omit(perm.dat) %>% mutate(across(everything(),~as.numeric(.))) %>%#drop lead/lag NAs and ensure all columns are numeric
      mutate(across(everything(),~scale(.)))
      
      #calculate permuted correlation
      ccf.perm.tmp <- ccf(c(perm.dat[,2]),c(perm.dat[,3]),lag.max = span,plot = F)
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
