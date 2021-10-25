###########################################################################################
### cross.granger function ###
###########################################################################################
# @ts = hypothesised causative variable
# @comp.ts = response variable
# @span = number of lags to test causality over
# @method = "var"/"raw". "var" utilises vector autoregression residuals (package::vars)
            #whereas "raw" utilises the raw variables (package::lmtest)
# @covariates = a matrix of covariates
# @ic = information criterion to select optimum lag. 
        # "AIC" = Akaike information criterion, "HQ" = Hannan–Quinn information criterion,
        # "SC" = Bayesian information criterion, "FPE" = Final Prediction Error criterion

cross.granger.ic <- function(ts,comp.ts,span,method = c("var","TY","TYrestrict"),covariates = NULL,ic = c("AIC", "HQ", "SC", "FPE")){
  
  require(vars)
  require(aod)
  
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
  
  if(is.null(covariates)){ 
      lag.info.fwd <- vars::VARselect(sub.dat[,c("ts","comp.ts")],lag.max=span, type = "both",season = 12) #identify optimal lag
      lag.info.rev <- vars::VARselect(sub.dat[,c("comp.ts","ts")],lag.max=span, type = "both",season = 12) #identify optimal lag
      
      if(ic %in% c("AIC")){
        kk.fwd <- lag.info.fwd$selection[1]
        kk.rev <- lag.info.rev$selection[1]
      }
      if(ic %in% c("HQ")){
          kk.fwd <- lag.info.fwd$selection[2]
          kk.rev <- lag.info.rev$selection[2]
        }
      if(ic %in% c("SC")){
        kk.fwd <- lag.info.fwd$selection[3]
        kk.rev <- lag.info.rev$selection[3]
      }
      if(ic %in% c("FPE")){
        kk.fwd <- lag.info.fwd$selection[4]
        kk.rev <- lag.info.rev$selection[4]
      }
      
      obsx <- vars::VAR(sub.dat[,c("ts","comp.ts")], p = kk.fwd+1,type = "both",season = 12) # fit VAR model forward
      obsy<-  vars::VAR(sub.dat[,c("comp.ts","ts")], p = kk.rev+1, type = "both",season = 12) # fit VAR model reverse
    
    }else{
      lag.info.fwd <- vars::VARselect(sub.dat[,c("ts","comp.ts")],lag.max=span, type = "both",season = 12,exogen = as.matrix(covariates)) #identify optimal lag
      lag.info.rev <- vars::VARselect(sub.dat[,c("comp.ts","ts")],lag.max=span, type = "both",season = 12,exogen = as.matrix(covariates)) #identify optimal lag
      
      if(ic %in% c("AIC")){
        kk.fwd <- lag.info.fwd$selection[1]
        kk.rev <- lag.info.rev$selection[1]
      }
      if(ic %in% c("HQ")){
        kk.fwd <- lag.info.fwd$selection[2]
        kk.rev <- lag.info.rev$selection[2]
      }
      if(ic %in% c("SC")){
        kk.fwd <- lag.info.fwd$selection[3]
        kk.rev <- lag.info.rev$selection[3]
      }
      if(ic %in% c("FPE")){
        kk.fwd <- lag.info.fwd$selection[4]
        kk.rev <- lag.info.rev$selection[4]
      }
      obsx <- suppressWarnings(vars::VAR(sub.dat[,c("ts","comp.ts")], p = kk.fwd+1,type = "both",season = 12,exogen = as.matrix(covariates))) # fit VAR model forward
      obsy<-  suppressWarnings(vars::VAR(sub.dat[,c("comp.ts","ts")], p = kk.rev+1, type = "both",season = 12,exogen = as.matrix(covariates))) # fit VAR model reverse
    
      } #warnings suppressed as covariates inputted without colnames kicks warning
  
  if(method == "var"){ #utilise vars::VAR and vars::causality
    gc.obsx <- vars::causality(obsx,cause = "ts",boot = T,boot.runs = 1000)$Granger # perform forward (null = x not Granger cause y)
    gc.obsy <- vars::causality(obsy,cause = "comp.ts",boot = T,boot.runs = 1000)$Granger # perform in reverse
    
    gc.df <- rbind(data.frame("measure" = "x_y","lag" = obsx$p, "F.value" = gc.obsx$statistic,"P.value" = gc.obsx$p.value,
                              "sig" = ifelse(gc.obsx$p.value <= 0.05,"*","")),
                   data.frame("measure" = "y_x","lag" = obsy$p, "F.value" = gc.obsy$statistic,"P.value" = gc.obsy$p.value,
                              "sig" = ifelse(gc.obsy$p.value <= 0.05,"*","")))%>% #return Wald test F stat, p-value and significance
             mutate(causality.direc = ifelse(grepl("^x_y",measure),"forward","reverse")) # classify directions
  }
  if(method == "TY"){ #TODA-YAMAMOTO IMPLEMENTATION 
    gc.obsx <- suppressWarnings(aod::wald.test(b=coef(obsx$varresult[[1]]), Sigma=vcov(obsx$varresult[[1]]), Terms=seq(from = 2,to = (2*kk.fwd), by=2)))
    gc.obsy <- suppressWarnings(aod::wald.test(b=coef(obsx$varresult[[2]]), Sigma=vcov(obsx$varresult[[2]]), Terms= seq(from = 1,to = (2*kk.fwd), by=2)))
    
    gc.df <- rbind(data.frame("measure" = "x_y","lag" = kk.fwd, "F.value" =  gc.obsx$result$chi2[1],"P.value" = gc.obsx$result$chi2[3],
                              "sig" = ifelse(gc.obsx$result$chi2[3] <= 0.05,"*","")),
                   data.frame("measure" = "y_x","lag" = kk.fwd, "F.value" = gc.obsy$result$chi2[1],"P.value" =gc.obsy$result$chi2[3],
                              "sig" = ifelse(gc.obsy$result$chi2[3] <= 0.05,"*","")))%>% #return Wald test F stat, p-value and significance
      mutate(causality.direc = ifelse(grepl("^x_y",measure),"forward","reverse")) # classify directions
  }
  if(method == "TYrestrict"){ #Toad-Yamamoto implementation restricted to single lag of interest 
    gc.obsx <- suppressWarnings(aod::wald.test(b=coef(obsx$varresult[[1]]), Sigma=vcov(obsx$varresult[[1]]), Terms=c(kk.fwd)))
    gc.obsy <- suppressWarnings(aod::wald.test(b=coef(obsx$varresult[[2]]), Sigma=vcov(obsx$varresult[[2]]), Terms= c(kk.fwd)))
    
    gc.df <- rbind(data.frame("measure" = "x_y","lag" = kk.fwd, "F.value" =  gc.obsx$result$chi2[1],"P.value" = gc.obsx$result$chi2[3],
                              "sig" = ifelse(gc.obsx$result$chi2[3] <= 0.05,"*","")),
                   data.frame("measure" = "y_x","lag" = kk.fwd, "F.value" = gc.obsy$result$chi2[1],"P.value" =gc.obsy$result$chi2[3],
                              "sig" = ifelse(gc.obsy$result$chi2[3] <= 0.05,"*","")))%>% #return Wald test F stat, p-value and significance
      mutate(causality.direc = ifelse(grepl("^x_y",measure),"forward","reverse")) # classify directions
  }
  
  return(gc.df)
}






cross.granger.unrestricted <- function(ts,comp.ts,span,method = c("var","raw"),covariates = NULL,ic=c("AIC","SC")){

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
    gc.df <- data.frame(statistic = c("x_y.lag","x_y.F","x_y.P","x_y.sig","x_y.aic","x_y.BPsig","y_x.lag","y_x.F","y_x.P","y_x.sig","y_x.aic","y_x.BPsig")) #initiate df

    for(kk in c(1:span)){
      gc_col_name <- paste0("run_", kk) # new column and colname for each lag

      if(is.null(covariates)){
        obsx <- vars::VAR(sub.dat[,c("ts","comp.ts")], p = kk, type = "both",season = 12) # fit VAR model forward
        obsy<-  vars::VAR(sub.dat[,c("comp.ts","ts")], p = kk, type = "both",season = 12) # fit VAR model reverse
        }else{
          obsx <- suppressWarnings(vars::VAR(sub.dat[,c("ts","comp.ts")], p = kk, type = "both",season = 12,exogen = as.matrix(covariates))) # fit VAR model forward
          obsy<-  suppressWarnings(vars::VAR(sub.dat[,c("comp.ts","ts")], p = kk, type = "both",season = 12,exogen = as.matrix(covariates))) # fit VAR model reverse
        } #warnings supressed as covariates inputted without colnames kicks warning
      
      hetero.testx <- lmtest::bptest(obsx$varresult$ts) # perform Breusch-Pagan test for heteroscedasticity
      hetero.testy <- lmtest::bptest(obsy$varresult$comp.ts)
      
      gc.obsx <- vars::causality(obsx,cause = "ts")$Granger # perform forward (null = x not Granger cause y)
      gc.obsy <- vars::causality(obsy,cause = "comp.ts")$Granger # perform in reverse

      if(ic == "AIC"){
      gc.out <- rbind(obsx$p,gc.obsx$statistic,gc.obsx$p.value,ifelse(gc.obsx$p.value <= 0.025,"*",""),AIC(obsx),hetero.testx$p.value,
                      obsy$p,gc.obsy$statistic,gc.obsy$p.value,ifelse(gc.obsy$p.value <= 0.025,"*",""),AIC(obsy),hetero.testy$p.value) #return Wald test F stat, p-value and AIC
      }else{
        gc.out <- rbind(obsx$p,gc.obsx$statistic,gc.obsx$p.value,ifelse(gc.obsx$p.value <= 0.025,"*",""),BIC(obsx),hetero.testx$p.value,
                        obsy$p,gc.obsy$statistic,gc.obsy$p.value,ifelse(gc.obsy$p.value <= 0.025,"*",""),BIC(obsy),hetero.testy$p.value) #return Wald test F stat, p-value and BIC
      }
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
