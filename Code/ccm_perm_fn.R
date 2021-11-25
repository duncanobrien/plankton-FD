###########################################################################################
### ccm.perm function ###
###########################################################################################
#@dat = matrix of data inputs. dat[,1] must be a numeric date vector,
# dat[,2] is the hypothesised causal variable,
# dat[,3] is the response variable of interest.
#@iter = number of permutations
#@span = number of lag time points to cross map across
#@return.raw = whether each cross skill at each lag should be returned
#@detrend.method = "none", "diff" (differenced), "lm" (standardised residuals of a lm between causal variable and time)
#@lag = if 'diff' detrending method selected, how many lags to difference over

ccm.perm <- function(dat,iter =999, span = 12,return.raw = F, 
                     detrend.method = c("none","diff","lm"), lag = 0){
  
  require(rEDM) # convergent cross mapping functions
  require(parallel) # parallel
  require(zoo) # na.approx function
  require(data.table) # data wrangling functions
  require(dplyr) # data wrangling functions
  require(foreach) # foreach and %dopar% functions
  
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
  
  if(any(is.na(dat[,2]))){
    dat[,2] <- zoo::na.approx(dat[,2],maxgap =3,na.rm=F)
  } # certain FD metrics may have NAs for single data points
  
  if(any(is.na(dat[,3]))){
    dat[,3] <- zoo::na.approx(dat[,3],maxgap =3,na.rm=F)
  } # certain state metrics may have NAs for single data points
  
  if(detrend.method == "diff"){
    sub.dat <- dat[(1+lag):nrow(dat),]
    sub.dat[,2] <- base::diff(dat[,2], lag = lag)
    sub.dat[,3] <- base::diff(dat[,3], lag = lag)
  }else if(detrend.method == "lm"){
    sub.dat <- dat
    sub.dat[,2] <- residuals(lm(sub.dat[,2] ~ as.numeric(sub.dat[,1]),na.action=na.exclude)) - 
      c(decompose(ts(sub.dat[,2],frequency = 12),type="additive")$seasonal)
    sub.dat[,3] <- residuals(lm(sub.dat[,3] ~ as.numeric(sub.dat[,1]),na.action=na.exclude)) - 
      c(decompose(ts(sub.dat[,3],frequency = 12),type="additive")$seasonal)
  }else{
    sub.dat <- dat
  }
  sub.dat <- na.omit(sub.dat) %>% mutate(across(everything(),~as.numeric(.))) %>% #drop lead/lag NAs and ensure all columns are numeric
    mutate(across(everything(),~scale(.)))
  
  ### Observed CCM ###
  simplex_FD <- rEDM::simplex(sub.dat[, c(1,2)], # identify optimal embedding dimension for FD i.e. when maximum info (rho) contained
                              E = 2:10, # range of possible embdedding dimensions
                              lib = c(1, nrow(sub.dat)), pred = c(1, nrow(sub.dat)),stats_only = T)%>%
    filter(as.numeric(rho) %in% max(as.numeric(rho))) %>% slice(1) %>% dplyr::select(E)
  
  simplex_state <- rEDM::simplex(sub.dat[,c(1,3)] , # identify optimal embedding dimension for system state
                                 E = 2:10, 
                                 lib = c(1, nrow(sub.dat)), pred = c(1, nrow(sub.dat)),stats_only = T)%>%
    filter(as.numeric(rho) %in% max(as.numeric(rho))) %>% slice(1) %>% dplyr::select(E)
  
  if(simplex_FD$E != simplex_state$E){ #select shared/maximum embedding dimension 
    simplex_use <- max(c(simplex_FD$E,simplex_state$E))
  }else{
    simplex_use <- simplex_FD$E
  }
  obs.params <- expand.grid(lib_column = colnames(sub.dat)[c(2,3)], target_column = colnames(sub.dat)[c(3,2)], 
                            tp = (-1*span):span) %>% # create frame of all possible combinations
    filter(lib_column != target_column)%>% # remove x->x and y->y matches
    slice(seq(1,n(),by=2)) # remove (reversed) duplicates 
  
  obs.ccm <- lapply(seq_len(nrow(obs.params)), function(kk){ # perform convergent cross mapping across lags
    rEDM::ccm(sub.dat[,c(2,3)], E = simplex_use, lib_sizes = floor(nrow(sub.dat)*0.75),  # single library size containing maximum information (as much of timeseries as possible)
              random_libs = FALSE, lib_column = paste(obs.params$lib_column[kk]), # lib_column = hypothesised causal variable
              target_column = paste(obs.params$target_column[kk]), # target_column = hypothesised response variable
              exclusion_radius = 0,tp = obs.params$tp[kk], silent = TRUE,	num_sample=1,RNGseed = 123)}) %>% # tp = lags to iterate over, num_sample for no repeats and RNGseed for reproducibility
    data.table::rbindlist(use.names=FALSE)%>%
    rename("x_y" = 2, "y_x" = 3) # rename columns to generic x_y/y_x for downstream wrangling
  
  obsx <- obs.ccm %>% filter(x_y == max(x_y) | tp == 0) %>% dplyr::select(-y_x) # extract observed cross skill at best lag and lag0 for x->y
  obsy <- obs.ccm %>% filter(y_x == max(y_x) | tp == 0) %>% dplyr::select(-x_y)# extract observed cross skill at best lag and lag0 for y->x
  
  obs.out <- data.frame(x_y.r0 = obsx$x_y[obsx$tp == 0],
                        x_y.lag = obsx$tp[obsx$x_y == max(obsx$x_y,na.rm=TRUE)])
  obs.out$x_y.skill <- obsx$x_y[obsx$tp == obs.out$x_y.lag]
  obs.out$y_x.r0  <- obsy$y_x[obsy$tp == 0]
  obs.out$y_x.lag  <- obsy$tp[obsy$y_x == max(obsy$y_x,na.rm=TRUE)]
  obs.out$y_x.skill <- obsy$y_x[obsy$tp == obs.out$y_x.lag] # create comparable out data frame
  
  
  ### Permuted CCM ###
  d1.ar1 <- stats::arima(stats::as.ts(dat[,2]), order = c(1, 0, 0),
                         seasonal = list(order = c(1, 0, 0), period = 12),
                         optim.control = list(maxit = 1000),method="ML")$coef[1]
  # fit ar1 model to hypothesised causal variable and extract autocorrelation coefficient
  
  perm.df <- data.frame(dat[,3]) %>% rename(!!sym(paste(colnames(dat)[3])) := 1)
  # initiate permutation data frame with response variable
  for (r in seq_len(iter)) {
    new_col_name <- paste0("perm_", r) # new column and colname for each permutation
    perm.df <- perm.df %>% 
      mutate(!!sym(new_col_name) :=  red.noise.ts(dat[,2],lag=d1.ar1,length=nrow(dat))) # permuted autocorrelated surrogates
  }  
  
  if(detrend.method == "diff"){
    perm.df <- perm.df %>%
      mutate(across(everything(), function(x){x - dplyr::lag(x,n=lag)})) %>%
      na.omit() %>% mutate(across(everything(),~scale(.)))
  }else if(detrend.method == "lm"){
    perm.df <- perm.df %>%
      mutate(across(everything(),function(.x){
        as.numeric(residuals(lm(.x ~ as.numeric(dat[,1]),na.action=na.exclude)) - 
                     c(decompose(ts(.x,frequency = 12),type="additive")$seasonal))})) %>%#standardised residuals
      #mutate(across(everything(),~rstandard(lm(.x ~ as.numeric(dat[,1]),na.action=na.exclude))))%>%
      na.omit() %>% mutate(across(everything(),~scale(.)))
  }else{
    perm.df <- perm.df %>% na.omit() %>% mutate(across(everything(),~scale(.)))
  }
  
  perm.tmp <- foreach::foreach(r = 1:iter,.combine = "rbind",.multicombine = F, .packages = c("dplyr","rEDM")) %dopar%{
    perm.params <- expand.grid(lib_column = c(paste(c("perm",r),collapse = "_"),colnames(sub.dat)[3]), 
                               target_column = c(colnames(sub.dat)[3],paste(c("perm",r),collapse = "_")), 
                               tp = (-1*span):span) %>% # create frame of all possible combinations
      filter(lib_column != target_column)%>% # remove x->x and y->y matches
      slice(seq(1,n(),by=2)) # remove (reversed) duplicates 
    
    perm.ccm <- lapply(seq_len(nrow(perm.params)), function(kk){ # perform convergent cross mapping across lags for each permutation (iter times)
      rEDM::ccm(perm.df[,c(paste(c("perm",r),collapse = "_"),colnames(sub.dat)[3])], 
                E = simplex_use, lib_sizes = floor(nrow(sub.dat)*0.75), 
                random_libs = FALSE, lib_column = paste(perm.params$lib_column[kk]), 
                target_column = paste(perm.params$target_column[kk]), exclusion_radius = 0,
                tp = perm.params$tp[kk], silent = TRUE, num_sample=1,RNGseed = 123)}) %>%
      data.table::rbindlist(use.names=FALSE)%>%
      rename("x_y" = 2, "y_x" = 3) 
    
    permx <- perm.ccm %>% filter(x_y == max(x_y) | tp == 0) %>% select(-y_x)
    permy <- perm.ccm %>% filter(y_x == max(y_x) | tp == 0) %>% select(-x_y)
    
    perm.out <- data.frame(x_y.r0 = permx$x_y[permx$tp == 0],
                           x_y.lag = permx$tp[permx$x_y == max(permx$x_y,na.rm=TRUE)])
    perm.out$x_y.skill <- permx$x_y[permx$tp == perm.out$x_y.lag]
    perm.out$y_x.r0  <- permy$y_x[permy$tp == 0]
    perm.out$y_x.lag  <- permy$tp[permy$y_x == max(permy$y_x,na.rm=TRUE)]
    perm.out$y_x.skill <- permy$y_x[permy$tp == perm.out$y_x.lag]
    
    return(perm.out)
  }
  parallel::stopCluster(cl)
  
  ### Combine information ###
  
  out <- matrix(NA,nrow = 3,ncol=7) # using permuted metrics, estimate permuted median and quantile of observed data for optimal lag and lag0
  
  out[1,] <- c("r0.skill",ecdf(perm.tmp$x_y.r0)(obs.out$x_y.r0),obs.out$x_y.r0,median(perm.tmp$x_y.r0),ecdf(perm.tmp$y_x.r0)(obs.out$y_x.r0),obs.out$y_x.r0,median(perm.tmp$y_x.r0))
  out[2,] <- c("max.skill",ecdf(perm.tmp$x_y.skill)(obs.out$x_y.skill),obs.out$x_y.skill,median(perm.tmp$x_y.skill),ecdf(perm.tmp$y_x.skill)(obs.out$y_x.skill),obs.out$y_x.skill,median(perm.tmp$y_x.skill))
  out[3,] <- c("t.max.skill",ecdf(perm.tmp$x_y.lag)(obs.out$x_y.lag),obs.out$x_y.lag,median(perm.tmp$x_y.lag),ecdf(perm.tmp$y_x.lag)(obs.out$y_x.lag),obs.out$y_x.lag,median(perm.tmp$y_x.lag))
  
  out <- data.frame(out) %>% 
    stats::setNames(c("measure","x_y.quantile", "x_y.obs_value","x_y.median_perm_value","y_x.quantile",
                      "y_x.obs_value","y_x.median_perm_value")) %>%
    mutate(across(x_y.quantile:y_x.median_perm_value, ~as.numeric(as.character(.)))) # ensure numeric data not character
  
  if(isTRUE(return.raw)){ # return observed cross skill values with lag if requested
    out.ls <- list("perm.dens" = perm.tmp, "summary" = out,"raw.obs" = obs.ccm[,c(2,3,6)])
  }else{
    out.ls <- list("perm.dens" = perm.tmp, "summary" = out)
  }
  return(out.ls)
}



