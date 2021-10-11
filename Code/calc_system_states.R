### Calculate System State ###

## Preamble ##
require(tidyverse) # dplyr, ggplot etc.

source("Code/mvi_fn.R")
fisher.scripts <- list.files(path ="Code/fisher_information", pattern="*.R",full.names = T) 
          # identify Fisher information scripts to be sourced
purrr::walk(fisher.scripts, source) # source silently

##########################################################################################
## Read in Plankton Abundance Data ##
##########################################################################################
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_plankton_data.R")
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_environmental_data.R")
#phyo, zoo + combined raw concentrations at week, month and year groupings

source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Zurich/Data/zurich_plankton_data.R")
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Zurich/Data/LZ_environmental_data.R")

source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Madison/Data/madison_plankton_data.R")
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Madison/Data/madison_environmental_data.R")

source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Windermere/Data/windermere_plankton_data.R")
source("/Users/duncanobrien/Desktop/Academia/PhD/Data/Windermere/Data/windermere_environmental_data.R")

##########################################################################################
## Kinneret 'Fisher Information' and 'Multivariate Index of Variance' ##
##########################################################################################
sd.kinmth <- apply(plank.kin.combo.mth[,4:79], MARGIN = 2, FUN = sd)
kinFImth <- GFisher(seq(1,dim(plank.kin.combo.mth)[1],1),plank.kin.combo.mth[,4:79],
                 sost =  matrix(sd.kinmth,1,dim(plank.kin.combo.mth[,4:79])[2]), 
                 hwin = 12, winspace = 1 , TL =90)
kinFImth.dat <- data.frame(FI = kinFImth$FI,
                        date = plank.kin.combo.mth$Date[apply(kinFImth$t_win, MARGIN = 2, FUN = max)],
                        maxt =apply(kinFImth$t_win, MARGIN = 2, FUN = max))

kinmth.mvi <- data.frame(multi.var.index(df=plank.kin.combo.mth[,4:79],window = 12))%>%
  mutate(date = plank.kin.combo.mth$Date[maxt])

sd.kinyr <- apply(plank.kin.combo.yr[,2:77], MARGIN = 2, FUN = sd)
kinFIyr <- GFisher(seq(1,dim(plank.kin.combo.yr)[1],1),plank.kin.combo.yr[,2:77],
                    sost =  matrix(sd.kinyr,1,dim(plank.kin.combo.yr[,2:77])[2]), 
                    hwin = 5, winspace = 1 , TL =90)
kinFIyr.dat <- data.frame(FI = kinFIyr$FI,
                          date = plank.kin.combo.yr$Date[apply(kinFIyr$t_win, MARGIN = 2, FUN = max)],
                          maxt =apply(kinFIyr$t_win, MARGIN = 2, FUN = max))

kinyr.mvi <- data.frame(multi.var.index(df=plank.kin.combo.yr[,2:77],window = 5))%>%
  mutate(date = plank.kin.combo.yr$Date[maxt])

##########################################################################################
## Lower Zurich 'Fisher Information' and 'Multivariate Index of Variance' ##
##########################################################################################
sd.LZmth <- apply(plank.LZ.combo.mth[,4:205], MARGIN = 2, FUN = sd)
LZFImth <- GFisher(seq(1,dim(plank.LZ.combo.mth)[1],1),plank.LZ.combo.mth[,4:205],
                sost =  matrix(sd.LZmth,1,dim(plank.LZ.combo.mth[,4:205])[2]), 
                hwin = 12, winspace = 1 , TL =95)
LZFImth.dat <- data.frame(FI = LZFImth$FI,date = plank.LZ.combo.mth$date[apply(LZFImth$t_win, MARGIN = 2, FUN = max)],
                       maxt = apply(LZFImth$t_win, MARGIN = 2, FUN = max))

LZmth.mvi <- data.frame(multi.var.index(df=plank.LZ.combo.mth[,4:205],window = 12))%>%
  mutate(date = plank.LZ.combo.mth$date[maxt])

sd.LZyr <- apply(plank.LZ.combo.yr[,2:203], MARGIN = 2, FUN = sd)
LZFIyr <- GFisher(seq(1,dim(plank.LZ.combo.yr)[1],1),plank.LZ.combo.yr[,2:203],
                   sost =  matrix(sd.LZyr,1,dim(plank.LZ.combo.yr[,2:203])[2]), 
                   hwin = 5, winspace = 1 , TL =85)
LZFIyr.dat <- data.frame(FI = LZFIyr$FI,
                         date = plank.LZ.combo.yr$date[apply(LZFIyr$t_win, MARGIN = 2, FUN = max)],
                         maxt =apply(LZFIyr$t_win, MARGIN = 2, FUN = max))
LZyr.mvi <- data.frame(multi.var.index(df=plank.LZ.combo.yr[,2:203],window = 5))%>%
  mutate(date = plank.LZ.combo.yr$date[maxt])

##########################################################################################
## Mendota 'Fisher Information' and 'Multivariate Index of Variance' ##
##########################################################################################
sd.madmth <- apply(plank.mad.combo.mth[,4:224], MARGIN = 2, FUN = sd)
madFImth <- GFisher(seq(1,dim(plank.mad.combo.mth)[1],1),plank.mad.combo.mth[,4:224],
                 sost =  matrix(sd.madmth,1,dim(plank.mad.combo.mth[,4:224])[2]), 
                 hwin = 12, winspace = 1 , TL =95)
madFImth.dat <- data.frame(FI = madFImth$FI,date = plank.mad.combo.mth$date[apply(madFImth$t_win, MARGIN = 2, FUN = max)],
                        maxt = apply(madFImth$t_win, MARGIN = 2, FUN = max))

madmth.mvi <- data.frame(multi.var.index(df=plank.mad.combo.mth[,4:224],window = 12))%>%
  mutate(date = plank.mad.combo.mth$date[maxt])

sd.madyr <- apply(plank.mad.combo.yr[,2:222], MARGIN = 2, FUN = sd)
madFIyr <- GFisher(seq(1,dim(plank.mad.combo.yr)[1],1),plank.mad.combo.yr[,2:222],
                    sost =  matrix(sd.madyr,1,dim(plank.mad.combo.yr[,2:222])[2]), 
                    hwin = 5, winspace = 1 , TL =85)
madFIyr.dat <- data.frame(FI = madFIyr$FI,
                          date = plank.mad.combo.yr$date[apply(madFIyr$t_win, MARGIN = 2, FUN = max)],
                          maxt =apply(madFIyr$t_win, MARGIN = 2, FUN = max))

madyr.mvi <- data.frame(multi.var.index(df=plank.mad.combo.yr[,2:222],window = 5))%>%
  mutate(date = plank.mad.combo.yr$date[maxt])

##########################################################################################
## Windermere 'Fisher Information' and 'Multivariate Index of Variance' ##
##########################################################################################
sd.windmth <- apply(phyto_env.windmthdata[,c(4:20,22:25)], MARGIN = 2, FUN = sd)
windFImth <- GFisher(seq(1,dim(phyto_env.windmthdata)[1],1),phyto_env.windmthdata[,c(4:20,22:25)],
                  sost =  matrix(sd.windmth,1,dim(phyto_env.windmthdata[,c(4:20,22:25)])[2]), 
                  hwin = 12, winspace = 1 , TL =95)
windFImth.dat <- data.frame(FI = windFImth$FI,date =phyto_env.windmthdata$Date[apply(windFImth$t_win, MARGIN = 2, FUN = max)],
                         maxt =apply(windFImth$t_win, MARGIN = 2, FUN = max) )

windmth.mvi <- data.frame(multi.var.index(df=phyto_env.windmthdata[,c(4:20,22:25)],window = 12))%>%
  mutate(date = phyto_env.windmthdata$Date[maxt])

sd.windyr <- apply(phyto_env.windyrdata[,c(2:18,20:23)], MARGIN = 2, FUN = sd)
windFIyr <- GFisher(seq(1,dim(phyto_env.windyrdata)[1],1),phyto_env.windyrdata[,c(2:18,20:23)],
                     sost =  matrix(sd.windyr,1,dim(phyto_env.windyrdata[,c(2:18,20:23)])[2]), 
                     hwin = 5, winspace = 1 , TL =80)
windFIyr.dat <- data.frame(FI = windFIyr$FI,
                           date = phyto_env.windyrdata$Date[apply(windFIyr$t_win, MARGIN = 2, FUN = max)],
                           maxt =apply(windFIyr$t_win, MARGIN = 2, FUN = max))
windyr.mvi <- data.frame(multi.var.index(df=phyto_env.windyrdata[,c(2:18,20:23)],window = 5))%>%
  mutate(date = phyto_env.windyrdata$Date[maxt])

##########################################################################################
## Combine and save out all system states ##
##########################################################################################
all.system.states <- list(kin.mth = data.frame("data.source" = "Kinneret",
                                               "res" = "Month",
                                               "date" = plank_env.data.mth$Date,
                                               "density" = rowSums(plank_env.data.mth[,4:79]),
                                               "community" = prcomp(scale(plank_env.data.mth[,4:79]))$x[,1],
                                               "zp.ratio" = (rowSums(plank.kin.combo.mth[,49:79])/rowSums(plank.kin.combo.mth[,4:48])))%>%
                            left_join(kinFImth.dat,by="date") %>% left_join(kinmth.mvi,by=c("date","maxt")),
                          kin.yr = data.frame("data.source" = "Kinneret",
                                              "res" = "Year",
                                              "date" = plank_env.data.yr$Date,
                                              "density" = rowSums(plank.kin.combo.yr[,2:77]),
                                              "community" = prcomp(scale(plank_env.data.yr[,2:77]))$x[,1],
                                              "zp.ratio" = (rowSums(plank.kin.combo.yr[,47:77])/rowSums(plank.kin.combo.yr[,2:46])))%>%
                            left_join(kinFIyr.dat,by="date") %>% left_join(kinyr.mvi,by=c("date","maxt")),
                          mad.mth = data.frame("data.source" = "Mendota",
                                               "res" = "Month",
                                               "date" =plank.mad.combo.mth$date,
                                               "density" = rowSums(plank.mad.combo.mth[,4:224]),
                                               "community" = prcomp(scale(plank.mad.combo.mth[,4:224]))$x[,1],
                                               "zp.ratio" = rowSums(plank.mad.combo.mth[,201:224])/rowSums(plank.mad.combo.mth[,4:200]))%>%
                            left_join(madFImth.dat,by="date") %>% left_join(madmth.mvi,by=c("date","maxt")),
                          mad.yr = data.frame("data.source" = "Mendota",
                                              "res" = "Year",
                                              "date" = plank_env.madyrdata$date,
                                              "density" = rowSums(plank_env.madyrdata[,2:222]),
                                              "community" = prcomp(scale(plank_env.madyrdata[,2:222]))$x[,1],
                                              "zp.ratio" = rowSums(plank.mad.combo.yr[,199:222])/rowSums(plank.mad.combo.yr[,2:198]))%>%
                            left_join(madFIyr.dat,by="date") %>% left_join(madyr.mvi,by=c("date","maxt")),
                          LZ.mth = data.frame("data.source" = "Lower Zurich",
                                              "res" = "Month",
                                              "date" =plank.LZ.combo.mth$date,
                                              "density" = rowSums(plank.LZ.combo.mth[,4:205]),
                                              "community" = prcomp(scale(plank.LZ.combo.mth[,4:205]))$x[,1],
                                              "zp.ratio" = rowSums(plank.LZ.combo.mth[,179:205])/rowSums(plank.LZ.combo.mth[,4:178]))%>%
                            left_join(LZFImth.dat,by="date") %>% left_join(LZmth.mvi,by=c("date","maxt")),
                          LZ.yr = data.frame("data.source" = "Lower Zurich",
                                             "res" = "Year",
                                             "date" = plank_env.LZyrdata$date,
                                             "density" = rowSums(plank_env.LZyrdata[,2:203]),
                                             "community" = prcomp(scale(plank_env.LZyrdata[,2:203]))$x[,1],
                                             "zp.ratio" = rowSums(plank.LZ.combo.yr[1:33,177:203])/rowSums(plank.LZ.combo.yr[1:33,2:176]))%>%
                            left_join(LZFIyr.dat,by="date") %>% left_join(LZyr.mvi,by=c("date","maxt")),
                          wind.mth = data.frame("data.source" = "Windermere",
                                                "res" = "Month",
                                                "date" = phyto_env.windmthdata$Date,
                                                "density" = rowSums(phyto_env.windmthdata[,c(4:20,22:25)]),
                                                "community" = prcomp(scale(phyto_env.windmthdata[,c(4:20,22:25)]))$x[,1],
                                                "zp.ratio" = rowSums(phyto_env.windmthdata[,22:25])/rowSums(phyto_env.windmthdata[,4:20]))%>%
                            left_join(windFImth.dat,by="date") %>% left_join(windmth.mvi,by=c("date","maxt")),
                          wind.yr = data.frame("data.source" = "Windermere",
                                               "res" = "Year",
                                               "date" = phyto_env.windyrdata$Date,
                                               "density" = rowSums(phyto_env.windyrdata[,c(2:18,20:23)]),
                                               "community" = prcomp(scale(phyto_env.windyrdata[,c(2:18,20:23)]))$x[,1],
                                               "zp.ratio" = rowSums(phyto_env.windyrdata[,20:23])/rowSums(phyto_env.windyrdata[,2:18]))%>%
                            left_join(windFIyr.dat,by="date") %>% left_join(windyr.mvi,by=c("date","maxt")))
save(all.system.states,file = "Data/all.system.states.RData")

