###########################################################################################
### cross.granger function ###
###########################################################################################
# @ts = hypothesised causative variable
# @comp.ts = response variable
# @span = number of lags to test causality over
# @method - "var"/"raw". "var" utilises vector autoregression residuals (package::vars)
            #whereas "raw" utilises the raw variables (package::lmtest)

cross.granger <- function(ts,comp.ts,span,method = c("var","raw"),covariates = NULL){
  
  if(is.null(covariates)){
  sub.dat <- data.frame("ts" = ts, "comp.ts" = comp.ts) # construct data frame. Data frame required in case input ts/comp.ts are scaled
  #frm.x_y <- formula(paste("comp.ts ~ ", "ts"))
  #frm.y_x <- formula(paste("ts ~ ", "comp.ts"))
  }else{
    sub.dat <-  data.frame("ts" = ts, "comp.ts" = comp.ts) %>%
      cbind(covariates)
  #  frm.x_y <- formula(paste("comp.ts ~ ", paste0(c("ts", colnames(covariates)), collapse = " + ")))
  #  frm.y_x <- formula(paste("ts ~ ", paste0(c("comp.ts", colnames(covariates)), collapse = " + ")))
   }

  if(method == "var"){ #utilise vars::VAR and vars::causality
    gc.df <- data.frame(statistic = c("x_y.F","x_y.P","x_y.sig","x_y.aic","y_x.F","y_x.P","y_x.sig","y_x.aic")) #initiate df
    
    for(kk in c(1:span)){
      gc_col_name <- paste0("lag_", kk) # new column and colname for each lag
    
      if(is.null(covariates)){ 
        obsx <- vars::VAR(sub.dat[,c("ts","comp.ts")], p = kk, type = "both",season = 12) # fit VAR model forward
        obsy<-  vars::VAR(sub.dat[,c("comp.ts","ts")], p = kk, type = "both",season = 12) # fit VAR model reverse
        }else{
          obsx <- vars::VAR(sub.dat[,c("ts","comp.ts")], p = kk, type = "both",season = 12,exogen = as.matrix(covariates)) # fit VAR model forward
          obsy<-  vars::VAR(sub.dat[,c("comp.ts","ts")], p = kk, type = "both",season = 12,exogen = as.matrix(covariates)) # fit VAR model reverse
        }
      gc.obsx <- vars::causality(obsx,cause = "ts")$Granger # perform forward (null = x not Granger cause y)
      gc.obsy <- vars::causality(obsy,cause = "comp.ts")$Granger # perform in reverse

      gc.out <- rbind(gc.obsx$statistic,gc.obsx$p.value,ifelse(gc.obsx$p.value <= 0.025,"*",""),AIC(obsx),
                      gc.obsy$statistic,gc.obsy$p.value,ifelse(gc.obsy$p.value <= 0.025,"*",""),AIC(obsy)) #return Wald test F stat, p-value and AIC
      gc.df <- gc.df %>% 
        mutate(!!sym(gc_col_name) := gc.out) # bind to returned data frame with labelled column
    }
  }else{
    
    gc.df <- data.frame(statistic = c("x_y.F","x_y.P","x_y.sig","y_x.F","y_x.P","y_x.sig")) #initiate df
    
    for(kk in c(1:span)){
      gc_col_name <- paste0("lag_", kk) # new column and colname for each lag
      gc.obsx <- lmtest::grangertest(x=sub.dat[,"ts"],y=sub.dat[,"comp.ts"],order=kk) # perform forward (null = x not Granger cause y)
      gc.obsy <- lmtest::grangertest(x=sub.dat[,"comp.ts"],y=sub.dat[,"ts"],order=kk) # perform in reverse
      #gc.obsx <- lmtest::grangertest(formula = frm.x_y ,data = list(sub.dat),order=kk) # perform forward (null = x not Granger cause y)
      #gc.obsy <- lmtest::grangertest(formula = frm.y_x,data = list(sub.dat),order=kk) # perform in reverse
      
      gc.out <- rbind(gc.obsx$F[2],gc.obsx$`Pr(>F)`[2],ifelse(gc.obsx$`Pr(>F)`[2] <= 0.025,"*",""),
                      gc.obsy$F[2],gc.obsy$`Pr(>F)`[2],ifelse(gc.obsy$`Pr(>F)`[2] <= 0.025,"*","")) #return Wald test F stat and p-value
      gc.df <- gc.df %>% 
        mutate(!!sym(gc_col_name) := gc.out) #bind to returned df with labelled column
    }
  }
  return(gc.df)
}
