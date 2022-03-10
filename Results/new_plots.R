test.dat <- summary.ccf.mth1 %>% filter(measure %in% "r0.ccf") %>%
  group_by(state.metric,FD.metric,troph) %>% mutate(upr = max(obs.value),lwr = min(obs.value),mid = median(obs.value))%>%
  ungroup()
  
str(summary.ccf.mth1)

pdf(file="C:/Users/ul20791/Downloads/sdob_boxplot.pdf",
    width=8, height = 5)  
pccf.lag0.fin <- ggplot(test.dat,
                        aes(x=state.metric,y=obs.value,col=FD.metric))+
  # geom_rect(data = layer_data(pccf.lag0.1, 1L),
  #              aes(x = x, y = y,xmin = xmin, xmax=xmax, ymin = ymin,ymax=ymax),
  #              color = "black", size = 0.5,position = position_dodge(width=0.75))
  #geom_pointrange(aes(y=mid,ymin =lwr,ymax=upr,group=FD.metric),position= position_dodge(width=0.75),shape = 22,col="black",size=0.5)+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3,position = position_dodge(width=0.75))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=FD.metric),alpha=1,col="black",size=0.3,outlier.shape = NA)+
  
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=FD.metric,shape = system),size=3.5,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=FD.metric,shape = system),size=3.5,,col="black")+
  facet_wrap(~troph)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,0.7),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  ylab("Cross correlation") + xlab("System state proxy")+
   guides(size = guide_legend(order = 1), 
          col = guide_legend(order = 2),
          fill = guide_legend(order = 2),
          shape = guide_legend(order = 3))+
  theme_bw()
pccf.lag0.fin
dev.off()




pccf.lag0.1 <- ggplot(filter(test.dat,measure %in% "r0.ccf"),aes(x=state.metric,y=obs.value,col=FD.metric))+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3)+
  #geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_line(aes(linetype = system),position = position_dodge(width = 0.75), 
            alpha = 0, show.legend = FALSE)+
  #facet_wrap(~troph)
facet_grid(system~troph)


pccf.lag0.fin <- pccf.lag0.1 + 
  #geom_linerange(aes(x=state.metric,y = 0,ymin = lwr,ymax=upr,group=FD.metric),position=position_dodge(width=0.75),col="black",linetype = "dashed")+
  geom_segment(data = layer_data(pccf.lag0.1, 1L),
                                            aes(x = xmin, xend=xmax, y = 0,yend=0, group = linetype),
                                            color = "black", size = 0.5,position = position_dodge(width=0.75))+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,fill=FD.metric,group=FD.metric),size=3.5,stroke = 1,shape=21)+
   geom_point(position=position_dodge(width=0.75),
              aes(fill=NULL,group=FD.metric),size=3.5,stroke = 1,shape=21)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_point(data = obs.cor.lag0.state.tab,aes(x=state.metric,y=mean.cor,group=FD.metric,size="mean"),position=position_dodge(width=0.75),col="black",alpha=1,shape = "-",fill="black")+
  geom_linerange(data = obs.cor.lag0.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric),position=position_dodge(width=0.75),col="black",linetype = "solid")+
  #geom_errorbar(data = obs.cor.lag0.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric), width=0.8,
  #             position=position_dodge(0.75),col="black",size=0.5) +
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  scale_size_manual(values = c(8),breaks = c("mean"),name = NULL, labels = c("Mean cross\ncorrelation"),
                    guide = guide_legend(override.aes = list(fill = c("black"),shape=c("*"),size=c(8),alpha = c(1))))+
  ylab("Cross correlation") + xlab("System state proxy")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme_bw()

pccf.lag0.fin



FDIS.ccf <- ggplot(filter(summary.ccf.mth1,FD.metric %in% "FDis" & measure %in% c("r0.ccf","absmax.ccf")),
                        aes(x=measure,y=obs.value,col=measure))+
  # geom_rect(data = layer_data(pccf.lag0.1, 1L),
  #              aes(x = x, y = y,xmin = xmin, xmax=xmax, ymin = ymin,ymax=ymax),
  #              color = "black", size = 0.5,position = position_dodge(width=0.75))
  #geom_pointrange(aes(y=mid,ymin =lwr,ymax=upr,group=FD.metric),position= position_dodge(width=0.75),shape = 22,col="black",size=0.5)+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3,position = position_dodge(width=0.75))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  # geom_point(position=position_dodge(width=0.75),
  #            aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  # geom_point(position=position_dodge(width=0.75),
  #            aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  #facet_grid(troph~state.metric)+
  
  
  
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
                      
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("Lagx","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("Lagx","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_x_discrete(labels=c("r0.ccf" = "Lag0", "absmax.ccf" = "Lagx"),limits=c("r0.ccf","absmax.ccf"))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  ylab("Cross correlation") + xlab("Time series comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FDis') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0,0), "pt"))


FEVE.ccf <- ggplot(filter(summary.ccf.mth1,FD.metric %in% "FEve" & measure %in% c("r0.ccf","absmax.ccf")),
                   aes(x=measure,y=obs.value,col=measure))+
  # geom_rect(data = layer_data(pccf.lag0.1, 1L),
  #              aes(x = x, y = y,xmin = xmin, xmax=xmax, ymin = ymin,ymax=ymax),
  #              color = "black", size = 0.5,position = position_dodge(width=0.75))
  #geom_pointrange(aes(y=mid,ymin =lwr,ymax=upr,group=FD.metric),position= position_dodge(width=0.75),shape = 22,col="black",size=0.5)+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3,position = position_dodge(width=0.75))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,,col="black")+
  #facet_grid(troph~state.metric)+
  
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("Lagx","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("Lagx","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_x_discrete(labels=c("r0.ccf" = "Lag0", "absmax.ccf" = "Lagx"),limits=c("r0.ccf","absmax.ccf"))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  ylab("Cross correlation") + xlab("Time series comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FEve') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0,0), "pt"))


FRIC.ccf <- ggplot(filter(summary.ccf.mth1,FD.metric %in% "FRic" & measure %in% c("r0.ccf","absmax.ccf")),
                   aes(x=measure,y=obs.value,col=measure))+
  # geom_rect(data = layer_data(pccf.lag0.1, 1L),
  #              aes(x = x, y = y,xmin = xmin, xmax=xmax, ymin = ymin,ymax=ymax),
  #              color = "black", size = 0.5,position = position_dodge(width=0.75))
  #geom_pointrange(aes(y=mid,ymin =lwr,ymax=upr,group=FD.metric),position= position_dodge(width=0.75),shape = 22,col="black",size=0.5)+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3,position = position_dodge(width=0.75))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,,col="black")+
  #facet_grid(troph~state.metric)+
  
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x",switch="x")+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("Lagx","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("Lagx","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_x_discrete(labels=c("r0.ccf" = "Lag0", "absmax.ccf" = "Lagx"),limits=c("r0.ccf","absmax.ccf"),position = "top")+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  ylab("Cross correlation") + xlab("Time series comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FRic') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0,0), "pt"))


pdf(file="C:/Users/ul20791/Downloads/sdob_boxplot.pdf",
    width=13, height = 13)
FDIS.ccf + FEVE.ccf + FRIC.ccf + plot_layout(nrow = 3,guides = "collect")+
                    plot_annotation(tag_levels = c('a')) & 
                    theme(plot.tag = element_text(face = "bold"))&
  theme(legend.position='right')
dev.off()  

