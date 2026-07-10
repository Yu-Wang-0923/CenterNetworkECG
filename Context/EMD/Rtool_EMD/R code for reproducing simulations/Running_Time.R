## After saving the results: plot the results
## Plot Results
load("Res/SimUnweightedEA.RData")
library(dplyr)
library(ggplot2)
Result.final <- Result.final %>%
  mutate(Method =variable,Vertex = N) #%>%
#  select(-c("variable","N"))

result_summary <-
  summarySE(Result.final, measurevar="value",
            groupvars=c("Criterion","n","Vertex",
                        "out.prop","Scenario","Method"))

result_summary <- result_summary %>%
  mutate(out.prop=paste0("r=",out.prop))

if(!dir.exists("Figure")){
  dir.create("Figure")
}

# Plot Time Comparison
result_summary <- result_summary %>%
  mutate(Vertex = paste0("N=",Vertex))
ggplot(result_summary %>% filter(Criterion=="Time"&
                                   out.prop =="r=0.05"&
                                   Scenario =="Scenario I"),
       aes(x=n, y=log(mean), colour=Method, group=Method,
           shape=Method,linetype=Method)) +
  theme_bw()+
  facet_grid(.~Vertex,scales = "free")+
  geom_errorbar(aes(ymin=log(mean-ci), ymax=log(mean+ci)),  width=.1) +
  geom_line() +
  geom_point()+
  labs(y="Time for Ranking")+
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
  scale_color_manual(values=c("eMD"="#4DAF4A",
                              "Ham"="#984EA3","eMD3"="#B15928",
                              "eMDw" =  "#1F78B4" ,"Frob"="#FF7F00"))#+
#scale_shape_manual(values=c("eMD"=19,"Ham"=17,"eMD3"=4))
ggsave(paste0("Figure/Time_S1.eps"),width = 10, height = 6)



ggplot(result_summary %>% filter(Criterion=="EA"&
                                   out.prop =="r=0.05"&
                                   Scenario =="Scenario I"),
       aes(x=n, y=log(mean), colour=Method, group=Method,
           shape=Method,linetype=Method)) +
  theme_bw()+
  facet_grid(.~Vertex,scales = "free")+
  geom_errorbar(aes(ymin=log(mean-ci), ymax=log(mean+ci)),  width=.1) +
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
  scale_color_manual(values=c("eMD"="#4DAF4A",
                              "Ham"="#984EA3","eMD3"="#B15928",
                              "eMDw" =  "#1F78B4" ,"Frob"="#FF7F00"))
View(result_summary %>% filter(Criterion=="EA"&
                                 out.prop =="r=0.05"&
                                 Scenario =="Scenario I"))
#+
#scale_shape_manual(values=c("eMD"=19,"Ham"=17,"eMD3"=4))
ggsave(paste0("Figure/Time_EA.eps"),width = 10, height = 6)

