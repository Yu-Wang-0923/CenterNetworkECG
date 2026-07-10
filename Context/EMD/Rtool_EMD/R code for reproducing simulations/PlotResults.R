library(dplyr)
library(ggplot2)

###==================
### Plot EA for Unweighted network
## After saving the results: plot the results
## Plot Results
if(!dir.exists("Figure")){
  dir.create("Figure")
}
load("Res/SimUnweightedEA.RData")
Result.final <- Result.final %>%
  mutate(Method =variable,Vertex = N) %>%
  dplyr::select(-c("variable","N"))

result_summary <-
  summarySE(Result.final, measurevar="value",
            groupvars=c("Criterion","n","Vertex",
                        "out.prop","Scenario","Method"))

result_summary <- result_summary %>%
  mutate(out.prop=paste0("r=",out.prop))

result_summary.uw <- result_summary

pd <- position_dodge(3)
## Plot EA
for(Vertex.Num in c(30,200)){
  ggplot(result_summary %>% filter(Criterion=="EA"&Vertex==Vertex.Num ),
         aes(x=n, y=(mean), colour=Method, group=Method,
             shape=Method,linetype=Method)) +
    theme_bw()+
    facet_grid(out.prop~Scenario,scales = "free")+
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),  width=.1,position=pd) +
    geom_line() +
    geom_point()+
    labs(y="EA")+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    theme(legend.position="bottom")+
    theme(legend.text = element_text(size=15))+
    theme(strip.text.x = element_text(size=20),
          strip.text.y = element_text(size=20),
          axis.title=element_text(size=20),
          axis.title.x =element_text(size=20),
          axis.title.y=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20))+
    scale_color_manual(values=c("eMD"="#4DAF4A","Ham"="#984EA3","eMD3"="#B15928"))+
    scale_shape_manual(values=c("eMD"=19,"Ham"=17,"eMD3"=4))
  ggsave(paste0("Figure/EA_",Vertex.Num,".eps"),width = 10, height = 6)
}


###==================
### Plot EA for Unweighted network
## After saving the results: plot the results
## Plot Results
load("Res/SimWeightedEA.RData")
Result.final <- Result.final %>%
  mutate(Method =variable,Vertex = N) %>%
  dplyr::select(-c("variable","N"))

result_summary <-
  summarySE(Result.final, measurevar="value",
            groupvars=c("Criterion","n","Vertex",
                        "out.prop","Scenario","Method"))

result_summary <- result_summary %>%
  mutate(out.prop=paste0("r=",out.prop))

result_summary.w <- result_summary

## Plot EA
for(Vertex.Num in c(30,200)){
  ggplot(result_summary %>% filter(Criterion=="EA"&Vertex==Vertex.Num ),
         aes(x=n, y=(mean), colour=Method, group=Method,
             shape=Method,linetype=Method)) +
    theme_bw()+
    facet_grid(Scenario~out.prop,scales = "free")+
    #facet_wrap(vars(out.prop, Scenario), scales = "free_y")+
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),  width=.1) +
    geom_line() +
    geom_point()+
    labs(y="EA")+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    theme(legend.position="bottom")+
    theme(legend.text = element_text(size=15))+
    theme(strip.text.x = element_text(size=20),
          strip.text.y = element_text(size=20),
          axis.title=element_text(size=20),
          axis.title.x =element_text(size=20),
          axis.title.y=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20))+
    scale_color_manual(values=c("eMD"="#4DAF4A",Frob = "#FF7F00","Ham"="#984EA3",
                                "eMD3"="#B15928"))+
    scale_shape_manual(values=c("eMD"=19,Frob = 1,"Ham"=17,"eMD3"=4))
  ggsave(paste0("Figure/EA_W_",Vertex.Num,".eps"),width = 10, height = 6)
}

##================
# Time Comparison

result_summary.full <- rbind(data.frame(result_summary.uw,Type = "Unweight"),
                             data.frame(result_summary.w,Type = "Weight"))
result_summary.full$Vertex <- factor(paste0("N=",result_summary.full$Vertex),
                                     levels = c("N=30","N=100","N=200"))
ggplot(result_summary.full %>%
         filter(Criterion=="Time"&Scenario %in% c("Scenario I", "Setting I") &
                out.prop == "r=0.05"),
       aes(x=n, y=log(mean), colour=Method, group=Method,
           shape=Method,linetype=Method)) +
  theme_bw()+
  facet_grid(Type~Vertex,scales = "free")+
  geom_errorbar(aes(ymin=log(mean-ci), ymax=log(mean+ci)),  width=.1) +
  geom_line() +
  geom_point()+
  labs(y="Logarithm of Time")+
  theme(plot.title = element_text(hjust = 0.5,size=20))+
  theme(legend.position="bottom")+
  theme(legend.text = element_text(size=15))+
  theme(strip.text.x = element_text(size=20),
        strip.text.y = element_text(size=20),
        axis.title=element_text(size=20),
        axis.title.x =element_text(size=20),
        axis.title.y=element_text(size=20),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=20))+
  scale_color_manual(values=c("eMD"="#4DAF4A",Frob = "#FF7F00","Ham"="#984EA3",
                              "eMD3"="#B15928"))+
  scale_shape_manual(values=c("eMD"=19,Frob = 1,"Ham"=17,"eMD3"=4))
