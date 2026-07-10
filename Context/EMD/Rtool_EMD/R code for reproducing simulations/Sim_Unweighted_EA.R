
#####Scenario I, r=0.05
library(netEMD)
source("main_function.R", echo=F)
source("Gen_unweighted_network.R", echo=F)

reptime = 200
Result.final = NULL
cl <- makePSOCKcluster(50)
registerDoParallel(cl)
# Scenario = "Scenario I";N=30;out.prop = 0.05;n=50
for(Scenario in c("Scenario I","Scenario II","Scenario III")){
  for(N in c(30,100,200)){
    for(out.prop in c(0.05,0.1)){
      for(n in seq(50,500,50)){
        if(N==100&out.prop!=0.05&Scenario!="Scenario I"){
          break
        }
        cat(Scenario," r=",out.prop," N=",N," n=",n,"\n")
        To.numeric.matrix <- function(X){return(matrix(as.numeric(X),dim(X)))}

        sbm.res.tmp <- foreach(i=1:reptime,.combine='rbind') %dopar% {
          library(MASS)
          #library(nevada)
          library(igraph)
          set.seed(i)
          if(Scenario == "Scenario I"){

            n1 <- ceiling(out.prop*n)
            x<-1:(n1-1)
            num<-sample(x,1)
            family_SBM<-Gen_Data.SBM(N, n, k=3, p.block=NULL,
                                     p.in=0.16, p.bt=0.075,
                                     p.in.per=0.05, p.bt.per=0.02)
            real_center<-family_SBM[1,,]
            real_center<-To.numeric.matrix(real_center)
            family_SBM[1:num,,]<-Gen_Data.SBM(N,num, k=3, p.block=NULL,
                                              p.in=0.5, p.bt=0.2,
                                              p.in.per=0.05, p.bt.per=0.02)
            family_SBM[(num+1):n1,,]<-Gen_Data.ER(N, n1-num, p=0.5, alpha= 0.1)


            for (i in 1:n) {
              diag(family_SBM[i,,])<-0
              family_SBM[i,,]<-To.numeric.matrix(family_SBM[i,,])
            }
            family_SBM<-family_SBM[2:n,,]
          }
          if(Scenario == "Scenario II"){
            n1 <- ceiling(out.prop*n)
            x<-1:(n1-1)
            num<-sample(x,1)
            family_SBM<-Gen_Data.SBM(N, n, k=3, p.block=NULL,
                                     p.in=0.16, p.bt=0.075,
                                     p.in.per=0.05, p.bt.per=0.02)
            real_center<-family_SBM[1,,]
            real_center<-To.numeric.matrix(real_center)
            family_SBM[1:n1,,]<-Gen_Data.SW(N, n1, d=3, p=0.5, alpha = 0.01)


            for (i in 1:n) {
              diag(family_SBM[i,,])<-0
              family_SBM[i,,]<-To.numeric.matrix(family_SBM[i,,])
            }
          }
          if(Scenario == "Scenario III"){

            n1 <- ceiling(out.prop*n)
            x<-1:(n1-1)
            num<-sample(x,1)
            family_SBM<-Gen_Data.SBM(N, n, k=3, p.block=NULL,
                                     p.in=0.16, p.bt=0.075,
                                     p.in.per=0.05, p.bt.per=0.02)
            real_center<-family_SBM[1,,]
            real_center<-To.numeric.matrix(real_center)
            family_SBM[1:num,,]<-Gen_Data.RGG(N, num, radius=0.175, sigma= 0.2)
            family_SBM[(num+1):n1,,]<-matrix(1,N,N)

            for (i in 1:n) {
              diag(family_SBM[i,,])<-0
              family_SBM[i,,]<-To.numeric.matrix(family_SBM[i,,])
            }

          }
          vec<-vec_matrix(n-1,N,family_SBM)

          t0 <- Sys.time()
          Fre<-min_hamming(n-1,family_SBM)
          t1 <- Sys.time()
          EMD<-netEMD::net_EMD(family_SBM,unweighted = T)
          t2 <- Sys.time()
          EMD3<-netEMD::net_EMD(family_SBM, K =3, unweighted = T)
          t3 <- Sys.time()

          index1<-Fre[[2]][1]
          index2<-which(EMD==max(EMD))[1]
          index3<-which(EMD3==max(EMD3))[1]

          fre_ham<-dis_hamming(family_SBM[index1,,],real_center)
          EMD<-dis_hamming(family_SBM[index2,,],real_center)
          EMD3<-dis_hamming(family_SBM[index3,,],real_center)

          sim<-cbind(fre_ham,EMD,EMD3,
                     difftime(t1,t0,units = c("secs")),
                     difftime(t2,t1,units = c("secs")),
                     difftime(t3,t2,units = c("secs")))

        }
        sbm_res.tmp2 <- data.frame(rbind(sbm.res.tmp[,1:3],sbm.res.tmp[,4:6]),
                                   Criterion = rep(c("EA","Time"),each = reptime))
        colnames(sbm_res.tmp2)[1:3] <- c("Ham","eMD","eMD3")

        result_out1.tmp <-melt(sbm_res.tmp2,id = "Criterion")
        result_out1.tmp <- data.frame(result_out1.tmp,n=n,N=N,out.prop=out.prop,
                                      Scenario = Scenario)
        Result.final <- rbind(Result.final,result_out1.tmp)
        # We can store the temporal results
        save(result_out1.tmp,file=paste0("Res/",Scenario,"N",N,"n",n,"r",out.prop*100,".RData"))
        save(Result.final,file=paste0("Res2/",Scenario,"N",N,"n",n,"r",out.prop*100,".RData"))
      }
    }
  }
}

if(!dir.exists("Res")){
  dir.create("Res")
}
save.image("Res/SimUnweightedEA.RData")
stopCluster(cl)


# Plot Time Comparison
# result_summary <- result_summary %>%
#   mutate(Vertex = paste0("N=",Vertex))
# ggplot(result_summary %>% filter(Criterion=="Time"&
#                                     out.prop =="r=0.05"&
#                                  Scenario =="Scenario I"),
#        aes(x=n, y=log(mean), colour=Method, group=Method,
#            shape=Method,linetype=Method)) +
#   theme_bw()+
#   facet_grid(.~Vertex,scales = "free")+
#   geom_errorbar(aes(ymin=log(mean-ci), ymax=log(mean+ci)),  width=.1) +
#   geom_line() +
#   geom_point()+
#   labs(y="Time for Ranking")+
#   theme(plot.title = element_text(hjust = 0.5,size=20))+
#   theme(legend.position="bottom")+
#   theme(legend.text = element_text(size=15))+
#   theme(strip.text.x = element_text(size=20),
#         strip.text.y = element_text(size=20),
#         axis.title=element_text(size=20),
#         axis.title.x =element_text(size=20),
#         axis.title.y=element_text(size=20),
#         axis.text.y=element_text(size=20),
#         axis.text.x=element_text(size=20))+
#   scale_color_manual(values=c("eMD"="#4DAF4A","Ham"="#984EA3","eMD3"="#B15928"))+
#   scale_shape_manual(values=c("eMD"=19,"Ham"=17,"eMD3"=4))
# ggsave(paste0("Figure/Time_S1.eps"),width = 10, height = 6)

