library(netEMD)
source("main_function.R", echo=F)
source("Gen_weighted_network.R", echo=F)

reptime = 200
Result.final = NULL
cl <- makePSOCKcluster(32)
registerDoParallel(cl)
#registerDoSEQ(cl)
# Scenario = "Scenario I";N=30;out.prop = 0.05;n=50
for(Scenario in c("Setting I","Setting II","Setting III")){
  for(N in c(30,100,200)){
    for(out.prop in c(0,0.05,0.1)){
      for(n in c(50,seq(100,500,100))){
        if(N==100&out.prop!=0.05&Scenario!="Setting I"){
          break
        }
        cat(Scenario," r=",out.prop," N=",N," n=",n,"\n")
        if(file.exists(paste0("Res/w_",Scenario,"N",N,"n",n,"r",out.prop*100,".RData"))){
          load(paste0("Res/w_",Scenario,"N",N,"n",n,"r",out.prop*100,".RData"))
          load(paste0("Res2/",Scenario,"N",N,"n",n,"r",out.prop*100,".RData"))
          next
        }
        To.numeric.matrix <- function(X){return(matrix(as.numeric(X),dim(X)))}

        res.tmp <- foreach(i=1:reptime,.combine='rbind') %dopar% {
          library(MASS)
          #library(nevada)
          library(igraph)
          set.seed(i)
          if(Scenario == "Setting I"){
            n1 <- ceiling(out.prop*n)
            if(out.prop==0){
              family<-WRG_family(n+1,N)
              real_center<-family[1,,]
              gen_family <- family[2:(n+1),,]
            }else{
              family<-WRG_family(n,N)
            family_out<-family
            x<-sample(2:(n1-2),1)
            family_out[1:x,,]<-weight_out_family(x+1,N)[2:(x+1),,]
            family_out[(x+1):n1,,]<-weight_out1_family(n1-x+1,N)[2:(n1-x+1),,]
            real_center<-family[1,,]

            gen_family <- family_out
            }
          }
          if(Scenario == "Setting II"){
            n1 <- ceiling(out.prop*n)

            family<-weight_out_family(n+1,N)
            if(out.prop==0){
              real_center<-family[1,,]
              gen_family <- family[2:(n+1),,]
            }else{
            family_out<-family[2:(n+1),,]
            x<- sample(2:(n1-2),1)
            family_out[1:x,6:10,]<-family_out[1:x,6:10,]*100
            family_out[1:x,,6:10]<-family_out[1:x,,6:10]*100
            family_out[(x+1):n1,,11:20]<-family_out[(x+1):n1,,11:20]*10
            family_out[(x+1):n1,11:20,]<-family_out[(x+1):n1,11:20,]*10
            #out2<-matrix(0,N,N)
            for (j in 1:n) {
              diag(family_out[j,,])<-0
            }
            real_center<-family[1,,]
            gen_family <- family_out
            }
          }
          if(Scenario == "Setting III"){

            real_center<-full_c(N)
            diag(real_center)<-0

            n1 <- ceiling(out.prop*n)
            weight_family<- weight_family_full(n,N)

            if(out.prop==0){
              gen_family <- weight_family
            }else{
            weight_family[1:n1,,]<-weight_family_out6(n1,N)

            gen_family <- weight_family
            }
          }


          t0 <- Sys.time()
          Fre<-min_hamming(n-1,gen_family)
          t1 <- Sys.time()
          EMD<-netEMD::net_EMD(gen_family,unweighted = F)
          t2 <- Sys.time()
          EMD3<-netEMD::net_EMD(gen_family, K =3, unweighted = F)
          t3 <- Sys.time()
          Frob<-min_froben(n-1,gen_family)
          t4 <- Sys.time()

          index1<-Fre[[2]][1]
          index2<-which(EMD==max(EMD))[1]
          index3<-which(EMD3==max(EMD3))[1]
          index4<-Frob[[2]][1]

          fre_ham<-dis_hamming(gen_family[index1,,],real_center)
          EMD<-dis_hamming(gen_family[index2,,],real_center)
          EMD3<-dis_hamming(gen_family[index3,,],real_center)
          Frob<-dis_hamming(gen_family[index4,,],real_center)

          sim<-cbind(fre_ham,EMD,EMD3,Frob,
                     difftime(t1,t0,units = c("secs")),
                     difftime(t2,t1,units = c("secs")),
                     difftime(t3,t2,units = c("secs")),
                     difftime(t4,t3,units = c("secs")))

        }
        res.tmp2 <- data.frame(rbind(res.tmp[,1:4],res.tmp[,5:8]),
                                   Criterion = rep(c("EA","Time"),each = reptime))
        colnames(res.tmp2)[1:4] <- c("Ham","eMD","eMD3","Frob")

        result_out1.tmp <-melt(res.tmp2,id = "Criterion")
        result_out1.tmp <- data.frame(result_out1.tmp,n=n,N=N,out.prop=out.prop,
                                      Scenario = Scenario)
        Result.final <- rbind(Result.final,result_out1.tmp)
        # We can store the temporal results
        save(result_out1.tmp,file=paste0("Res/w_",Scenario,"N",N,"n",n,"r",out.prop*100,".RData"))
        save(Result.final,file=paste0("Res2/",Scenario,"N",N,"n",n,"r",out.prop*100,".RData"))
      }
    }
  }
}

if(!dir.exists("Res")){
  dir.create("Res")
}
save.image("Res/SimWeightedEA.RData")
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

