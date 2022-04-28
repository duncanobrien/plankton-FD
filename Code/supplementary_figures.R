### Supplementary Figures ###

## Preamble 
require(tidyverse) # dplyr, ggplot etc.
require(ggh4x) # facet_nested function
require(ggpubr) # multi-panel plots

load(file = "Results/ccf/raw_data/raw.ccf.RData") #load raw permutation cross correlation data
# RData required to reduce file size compared to .csv
summary.ccf <- read.csv("Results/ccf/raw_data/summary.ccf.csv") #load summary cross correlation data

load(file = "Results/ccm/raw_data/ccm_raw.RData") #load raw permutation convergent cross mapping data
summary.ccm <- read.csv(file ="Results/ccm/raw_data/ccm_summary.csv") #load summary convergent cross mapping data

###########################################################################
## Cross correlation supplementary violin plots ##
###########################################################################
pdf(file="Results/ccf/supplementary_figure3.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth1,aes(x = state.metric, y =  r0, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf[summary.ccf$measure %in% "r0.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf[summary.ccf$measure %in% "r0.ccf",], 
            aes(x = state.metric, y = 0.35,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  facet_grid(system~troph)+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF t0 between detrended FD and system state")
dev.off()

pdf(file="Results/ccf/supplementary_figure4.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth1,aes(x = state.metric, y =  abs.rmax, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf[summary.ccf$measure %in% "absmax.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf[summary.ccf$measure %in% "absmax.ccf",], 
            aes(x = state.metric, y = 0.45,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = summary.ccf[summary.ccf$measure %in% "t.absmax.ccf",], 
            aes(x = state.metric, y = 0.4,label = obs.value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  #facet_wrap(~troph,nrow=2,strip.position = "right")+
  facet_grid(system~troph)+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between detrended FD and system state")
dev.off()

###########################################################################
## Cross correlation lag density plot ##
###########################################################################

lag.density.df <- summary.ccf %>%
  dplyr::select(!c(quantile,median.perm.value,obs.difference,res))%>%
  filter(measure %in% c("absmax.ccf","t.absmax.ccf"))%>%
  pivot_wider(names_from = measure,values_from = obs.value)%>%
  ungroup()%>%
  mutate(t.absmax.ccf=ifelse(is.numeric(absmax.ccf) & !is.na(t.absmax.ccf),t.absmax.ccf,
                             dplyr::lead(t.absmax.ccf)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
  na.omit()

pdf(file="Results/ccf/supplementary_figure5.pdf",
    width=6, height = 5)
ggplot(lag.density.df,aes(x=t.absmax.ccf)) + 
  xlab("Lag") + ylab("Density")+
  geom_density(aes(linetype = "All"),col = "#90ADC6",fill="#90ADC6",alpha = 0.4,bw = 4)+
  geom_density(data = filter(lag.density.df,sig == "*"),aes(linetype = "Significant"),col = "#FAD02C",fill="#FAD02C",alpha = 0.4,bw = 4) + 
  geom_vline(xintercept =quantile(lag.density.df$t.absmax.ccf,probs = c(0.1,0.80)), linetype = "longdash",col="#90ADC6")+
  geom_vline(xintercept =quantile(filter(lag.density.df,sig == "*")$t.absmax.ccf,probs = c(0.1,0.80)), linetype = "longdash",col="#FAD02C")+
  theme_bw() +  scale_linetype_manual(values = c(1,1),
                                      guide = guide_legend(title = "Cross-correlation\ngroup",override.aes = list(size = 1,col=c("#90ADC6","#FAD02C"),fill = c("grey","#FAD02C"))))
dev.off()

###########################################################################
## Cross correlation statistic tables ##
###########################################################################
ccf.lag0.tab <- summary.ccf %>%
  filter(measure %in% "r0.ccf")%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(median.cor = median(obs.value),
            cor.se = sd(obs.value)/n(),
            nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(median.cor:cor.se,~round(.x,digits=3)))%>%
  rowwise()%>%
  mutate(median.cor = paste(c(median.cor,cor.se),collapse = " ± ")) %>%
  select(-c(cor.se))
write.csv(ccf.lag0.tab,file ="Results/ccf/ccf_tables/supplementary_table1.csv",row.names = F)

ccf.lagx.tab <- summary.ccf %>%
  dplyr::select(!c(quantile,median.perm.value,obs.difference,res))%>%
  filter(measure %in% c("absmax.ccf","t.absmax.ccf"))%>%
  pivot_wider(names_from = measure,values_from = obs.value)%>%
  ungroup()%>%
  mutate(t.absmax.ccf=ifelse(is.numeric(absmax.ccf) & !is.na(t.absmax.ccf),t.absmax.ccf,
                             dplyr::lead(t.absmax.ccf)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
  na.omit() %>%#drop duplicate rows
  group_by(troph,FD.metric,state.metric) %>%
  summarise(median.cor = median(absmax.ccf),
            cor.se = sd(absmax.ccf)/n(), 
            median.lag = median(t.absmax.ccf),
            lag.se = sd(t.absmax.ccf)/n(), 
            nsig=sum(sig %in% "*"),
            prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(median.cor:lag.se,~round(.x,digits=3)))%>%
  rowwise()%>%
  mutate(median.cor = paste(c(median.cor,cor.se),collapse = " ± "),
         median.lag = paste(c(median.lag,lag.se),collapse = " ± ")) %>%
  select(-c(cor.se,lag.se))
write.csv(ccf.lagx.tab,file ="Results/ccf/ccf_tables/supplementary_table2.csv",row.names = F)

###########################################################################
## Convergent cross mapping supplementary violin plots ##
###########################################################################

pdf(file="Results/ccm/supplementary_figure6.pdf",
    width=10, height = 8)
ggplot(raw.ccm,aes(x = state.metric, y =  y_x.r0, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccm[summary.ccm$measure %in% "r0.skill",],
             aes(x = state.metric, y = y_x.obs_value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccm[summary.ccm$measure %in% "r0.skill",], 
            aes(x = state.metric, y = 0.98,label = y_x.sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  scale_y_continuous(breaks = c(0.0,0.5,1.0))+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross map skill") + xlab("System state proxy")+   ggtitle("Permuted cross skill of system state mapping FD (forward)")
dev.off()

pdf(file="Results/ccm/supplementary_figure7.pdf",
    width=10, height = 8)
ggplot(raw.ccm,aes(x = state.metric, y =  x_y.r0, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccm[summary.ccm$measure %in% "r0.skill",],
             aes(x = state.metric, y = x_y.obs_value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccm[summary.ccm$measure %in% "r0.skill",], 
            aes(x = state.metric, y = 0.98,label = x_y.sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  scale_y_continuous(breaks = c(0,0.5,1.0))+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross map skill") + xlab("System state proxy")+   ggtitle("Permuted cross skill of FD mapping system state (reverse)")
dev.off()

pdf(file="Results/ccm/supplementary_figure8.pdf",
    width=10, height = 8)
ggplot(raw.ccm,aes(x = state.metric, y =  y_x.skill, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccm[summary.ccm$measure %in% "max.skill",],
             aes(x = state.metric, y = y_x.obs_value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccm[summary.ccm$measure %in% "max.skill",], 
            aes(x = state.metric, y = 0.98,label = y_x.sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = summary.ccm[summary.ccm$measure %in% "t.max.skill",], 
            aes(x = state.metric, y = 0.9,label = y_x.obs_value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  scale_y_continuous(breaks = c(0.0,0.5,1.0))+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross map skill") + xlab("System state proxy")+   ggtitle("Permuted cross skill of system state mapping FD (forward)")
dev.off()

pdf(file="Results/ccm/supplementary_figure10.pdf",
    width=10, height = 8)
ggplot(raw.ccm,aes(x = state.metric, y =  x_y.skill, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccm[summary.ccm$measure %in% "max.skill",],
             aes(x = state.metric, y = x_y.obs_value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccm[summary.ccm$measure %in% "max.skill",], 
            aes(x = state.metric, y = 0.98,label = x_y.sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = summary.ccm[summary.ccm$measure %in% "t.max.skill",], 
            aes(x = state.metric, y = 0.9,label = x_y.obs_value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  scale_y_continuous(breaks = c(0,0.5,1.0))+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross map skill") + xlab("System state proxy")+   ggtitle("Permuted cross skill of FD mapping system state (reverse)")
dev.off()

###########################################################################
## Convergent cross mapping lag density plot ##
###########################################################################
ccm.plot.df <- summary.ccm %>%
  filter(measure != "r0.skill")%>% # keep absolute highest cross map skill
  group_by(FD.metric,system,state.metric) %>%
  select(-c(x_y.median_perm_value,y_x.median_perm_value))%>% #drop unnecessary variables
  nest(x_y = c(measure:x_y.obs_value,x_y.sig), y_x =  c(measure,y_x.quantile:y_x.obs_value,y_x.sig))%>% #nest into x_y and y_x for ease of manipulation
  mutate(x_y = map(x_y, ~.x %>% 
                     pivot_longer(cols = c(x_y.quantile:x_y.obs_value), #
                                  names_to = c(".value"), names_prefix = "x_y.",
                                  names_repair = "unique", values_to = "value") %>% #pivot just x_y
                     mutate(lag = obs_value[2])%>% #add lag from t.absmax.skill row
                     slice(-2) %>% #drop t.absmax.skill row
                     rename(sig = x_y.sig) %>% #remove excess prefix
                     setNames(paste0('x_y.', names(.))))) %>% #add prefix for downstream wrangling
  mutate(y_x = map(y_x, ~.x %>% 
                     pivot_longer(cols = c(y_x.quantile:y_x.obs_value), #repeat process for y_x
                                  names_to = c(".value"),names_prefix = "y_x.",
                                  names_repair = "unique", values_to = "value") %>%
                     mutate(lag = obs_value[2])%>%
                     slice(-2) %>%
                     rename(sig = y_x.sig) %>%
                     setNames(paste0('y_x.', names(.))))) %>%
  unnest(cols = c(x_y, y_x),names_repair = "unique") %>% #unnest
  pivot_longer(c(x_y.measure:y_x.lag), #pivot using prefix as reference
               names_to = c("direc",".value"),
               names_pattern = "(.*)\\.(.*)" ) %>%
  mutate(causality.direc = ifelse(grepl("^x_y",direc),"FD map State","State map FD"))%>% 
  #classify directions. Is counterintuitive as here x_y represents x map y, 
  #which if significant implies y causes x (x contains information on y and 
  #therefore is causative)
  mutate(across(quantile:lag, ~as.numeric(.))) %>% #ensure values are numeric
  ungroup()

pdf(file="Results/ccm/supplementary_figure9.pdf",
    width=7, height = 5)
ggplot(ccm.plot.df,aes(x=lag)) + 
  xlab("Lag") + ylab("Density")+
  geom_density(aes(linetype = "All associations"),col = "#90ADC6",fill="#90ADC6",alpha = 0.4,bw = 4,size=0.8)+
  geom_density(data = filter(ccm.plot.df,sig == "*"),aes(linetype = "Significant\nassociations"),col = "#FAD02C",fill="#FAD02C",alpha = 0.4,bw = 4,size=0.8) + 
  geom_vline(xintercept =quantile(ccm.plot.df$lag,probs = c(0.1,0.80)), linetype = "longdash",col="#90ADC6")+
  geom_vline(xintercept =quantile(filter(ccm.plot.df,sig == "*")$lag,probs = c(0.1,0.80)), linetype = "longdash",col="#FAD02C")+
  theme_bw() +  
  geom_vline(xintercept =0, linetype = "solid",col="black",alpha = 0.8)+
  facet_wrap(~causality.direc)+
  scale_linetype_manual(values = c(1,1),
                        guide = guide_legend(title = "Causality group",override.aes = list(size = 1,col=c("#90ADC6","#FAD02C"),fill = c("#90ADC6","#FAD02C"), alpha = 0.2)))

dev.off()

###########################################################################
## Convergent cross mapping statistic tables ##
###########################################################################
ccm.lag0.y_x.tab <- summary.ccm %>%
  dplyr::select(!starts_with("x_y"))%>%
  filter(measure == "r0.skill")%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(median.cor = median(y_x.obs_value),
            cor.se = sd(y_x.obs_value)/n(),
            nsig=sum(y_x.sig %in% "*"),prop.sig = sum(y_x.sig %in% "*")/length(y_x.sig)) %>%
  mutate(across(median.cor:cor.se,~round(.x,digits=3)))%>%
  rowwise()%>%
  mutate(median.cor = paste(c(median.cor,cor.se),collapse = " ± ")) %>%
  select(-c(cor.se))
write.csv(ccm.lag0.y_x.tab,file ="Results/ccm/ccm_tables/supplementary_table3.csv",row.names = F)

ccm.lag0.x_y.tab <- summary.ccm %>%
  dplyr::select(!starts_with("y_x"))%>%
  filter(measure == "r0.skill")%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(median.cor = median(x_y.obs_value),
            cor.se = sd(x_y.obs_value)/n(),
            nsig=sum(x_y.sig %in% "*"),prop.sig = sum(x_y.sig %in% "*")/length(x_y.sig)) %>%
  mutate(across(median.cor:cor.se,~round(.x,digits=3)))%>%
  rowwise()%>%
  mutate(median.cor = paste(c(median.cor,cor.se),collapse = " ± ")) %>%
  select(-c(cor.se))
write.csv(ccm.lag0.x_y.tab,file ="Results/ccm/ccm_tables/supplementary_table4.csv",row.names = F)

ccm.lag0.comp <- summary.ccm %>%
  filter(measure == "r0.skill")%>%  
  mutate(forward = ifelse(y_x.sig == "*" & x_y.sig != "*",TRUE,FALSE),
         reverse = ifelse(x_y.sig == "*" & y_x.sig != "*",TRUE,FALSE),
         bidirec =  ifelse(x_y.sig == "*" & y_x.sig == "*",TRUE,FALSE),
         none =  ifelse(x_y.sig != "*" & y_x.sig != "*",TRUE,FALSE))%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(prop.forward=sum(forward == TRUE)/length(forward),
            prop.reverse=sum(reverse == TRUE)/length(reverse),
            prop.bidirec=sum(bidirec == TRUE)/length(bidirec),
            prop.none=sum(none == TRUE)/length(none))
write.csv(ccm.lag0.comp,file ="Results/ccm/ccm_tables/supplementary_table5.csv",row.names = F)

ccm.lagx.y_x.tab <- summary.ccm %>%
  dplyr::select(!starts_with("x_y"))%>%
  filter(measure %in% c("max.skill","t.max.skill"))%>%
  pivot_wider(names_from = measure,values_from = y_x.obs_value)%>%
  ungroup()%>%
  mutate(t.max.skill=ifelse(is.numeric(max.skill) & !is.na(t.max.skill),t.max.skill,
                            dplyr::lead(t.max.skill)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
  na.omit() %>%#drop duplicate rows
  group_by(troph,FD.metric,state.metric) %>%
  summarise(median.cor = median(max.skill),
            cor.se = sd(max.skill)/n(),
            median.lag = median(t.max.skill),lag.se = sd(t.max.skill)/n(), 
            nsig=sum(y_x.sig %in% "*"),prop.sig = sum(y_x.sig %in% "*")/length(y_x.sig))%>%
  mutate(across(median.cor:lag.se,~round(.x,digits=3)))%>%
  rowwise()%>%
  mutate(median.cor = paste(c(median.cor,cor.se),collapse = " ± "),
         median.lag = paste(c(median.lag,lag.se),collapse = " ± ")) %>%
  select(-c(cor.se,lag.se))
write.csv(ccm.lagx.y_x.tab,file ="Results/ccm/ccm_tables/supplementary_table6.csv",row.names = F)

ccm.lagx.x_y.tab <- summary.ccm %>%
  dplyr::select(!starts_with("y_x"))%>%
  filter(measure %in% c("max.skill","t.max.skill"))%>%
  pivot_wider(names_from = measure,values_from = x_y.obs_value)%>%
  ungroup()%>%
  mutate(t.max.skill=ifelse(is.numeric(max.skill) & !is.na(t.max.skill),t.max.skill,
                            dplyr::lead(t.max.skill)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
  na.omit() %>%#drop duplicate rows
  group_by(troph,FD.metric,state.metric) %>%
  summarise(median.cor = median(max.skill),
            cor.se = sd(max.skill)/n(),
            median.lag = median(t.max.skill),lag.se = sd(t.max.skill)/n(), 
            nsig=sum(x_y.sig %in% "*"),prop.sig = sum(x_y.sig %in% "*")/length(x_y.sig))%>%
  mutate(across(median.cor:lag.se,~round(.x,digits=3)))%>%
  rowwise()%>%
  mutate(median.cor = paste(c(median.cor,cor.se),collapse = " ± "),
         median.lag = paste(c(median.lag,lag.se),collapse = " ± ")) %>%
  select(-c(cor.se,lag.se))
write.csv(ccm.lagx.x_y.tab,file ="Results/ccm/ccm_tables/supplementary_table7.csv",row.names = F)

ccm.lagx.comp <- summary.ccm %>%
  filter(measure == "max.skill")%>%  
  mutate(forward = ifelse(y_x.sig == "*" & x_y.sig != "*",TRUE,FALSE),
         reverse = ifelse(x_y.sig == "*" & y_x.sig != "*",TRUE,FALSE),
         bidirec =  ifelse(x_y.sig == "*" & y_x.sig == "*",TRUE,FALSE),
         none =  ifelse(x_y.sig != "*" & y_x.sig != "*",TRUE,FALSE),
         diff.lag = filter(summary.ccm,measure == "t.max.skill")$x_y.obs_value - filter(summary.ccm,measure == "t.max.skill")$y_x.obs_value)%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(prop.forward=sum(forward == TRUE)/length(forward),
            prop.reverse=sum(reverse == TRUE)/length(reverse),
            prop.bidirec=sum(bidirec == TRUE)/length(bidirec),
            prop.none=sum(none == TRUE)/length(none),
            mean.lag = mean(diff.lag))
write.csv(ccm.lagx.comp,file ="Results/ccm/ccm_tables/supplementary_table8.csv",row.names = F)
