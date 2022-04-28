### Figure 4 ###

## Preamble 
require(tidyverse) # dplyr, ggplot etc.
require(ggh4x) # facet_nested function
require(ggpubr) # multi-panel plots

summary.ccm <- read.csv(file ="Results/ccm/raw_data/ccm_summary.csv")

ccm.lag.plot.df <- summary.ccm %>%
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

count.lag.ccmdf <- ccm.lag.plot.df %>%
  filter(sig == "*")%>% #only keep significant relationships
  group_by(state.metric,causality.direc,FD.metric,troph)%>%
  summarise(N = n()) #significant count per group

###########################################################################
## Plot spread of strongest lagged cross map skills ##
###########################################################################
pdf(file="Results/ccm/Figure4.pdf",
    width=13, height = 6)
ggpubr::ggpaired(ccm.lag.plot.df  %>%
                   mutate(pairing = interaction(state.metric,system,FD.metric,troph)) ,
                 x = "causality.direc", y = "lag", fill = rep(c("#A1B4FE","#FFE7A1"),10*3), 
                 id = "pairing",alpha = 0.1,point.size =0,  line.size = 0.5,
                 ggtheme = theme_bw(),linetype = "dashed") +
  scale_x_discrete(labels =c("reverse","forward"))+
  geom_hline(yintercept = 0,alpha = 0.3)+
  geom_point(aes(group=causality.direc,shape = system, alpha = as.factor(sig)), 
             position = position_dodge(width=0.75),fill = "black",size = 2)+
  geom_point(aes(group=causality.direc,shape = system,fill=NULL),col="black",
             position = position_dodge(width=0.75),size = 2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FDis"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FEve"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FRic"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  xlab("Causality direction") + ylab("Strongest lag (months)") + 
  ggh4x::facet_nested(FD.metric~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  scale_y_continuous(breaks = seq(-40,40,by=40),limits = c(-60,75))+
  guides(shape = guide_legend(order = 1))+
  theme(panel.spacing.x = unit(0.25, "lines"),
        strip.background = element_rect(fill="white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) 
dev.off()
