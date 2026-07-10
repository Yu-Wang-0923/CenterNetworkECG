#######################
#This file is for reproducing the results from Table 2 and Figure 5 of the paper

source("main_function.R", echo=F)
source("Gen_unweighted_network.R", echo=F)
## convert the network array to a matrix
vec_matrix<-function(n,N,array){
  vec_mat<-matrix(0,N^2,n)
  for (i in 1:n) {
    vec_mat[,i]<-as.vector(array[i,,])
  }
  return(vec_mat)
}

##EMD with equal importance score (unweighted networks)
netEMD02<-function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  vecq<-rowSums(data)
  vecp<-n-vecq
  matp<-matrix(rep(vecp,n),p,n)
  matq<-matrix(rep(vecq,n),p,n)
  matcp2<-matp*(matp-1)/2
  matcq2<-matq*(matq-1)/2
  resmat<- (matrix(choose(n,2),p,n)-(data*matcp2-(data-1)*matcq2))/choose(n,2)
  #resdepth<-as.vector(w%*%resmat/sum(w))
  resdepth<-colSums(resmat)/p
  return(resdepth)
}

###EMD with importance scores
netEMD_score<-function(data,w){
  p=dim(data)[1]
  n=dim(data)[2]
  vecq<-rowSums(data)
  vecp<-n-vecq
  matp<-matrix(rep(vecp,n),p,n)
  matq<-matrix(rep(vecq,n),p,n)
  matcp2<-matp*(matp-1)/2
  matcq2<-matq*(matq-1)/2
  resmat<- (matrix(choose(n,2),p,n)-(data*matcp2-(data-1)*matcq2))/choose(n,2)
  resdepth<-as.vector(w%*%resmat/sum(w))
  #resdepth<-colSums(resmat)/p
  return(resdepth)
}


##node centrality
score<-function(sfamily){
  n=dim(sfamily)[1]
  N<-dim(sfamily)[2]
  w<-array(0,dim=c(n,N,N))
  w2<-matrix(0,N,N)
  aa<-list()
  for (k in 1:n) {
    aa[[k]]<-eigen_centrality(
      as.undirected(graph_from_adjacency_matrix(sfamily[k,,])),
      directed = FALSE,
      scale = TRUE,
      weights = NULL,
      options = arpack_defaults
    )
    for (i in 1:(N-1)) {
      for (j in (i+1):N ) {
        w[k,i,j]<-w[k,j,i]<-(aa[[k]]$vector[i]+aa[[k]]$vector[j])/2
      }
    }
  }
  for(h in 1:n){
    w2<-w2+w[h,,]
  }
  w2<-w2/n
  impor_vec<-as.vector(w2)/sum(w2)
  #impor_vec[which(impor_vec!=0)]<-1/impor_vec[which(impor_vec!=0)]
  return(impor_vec)
}

###edge centrality
score_edge<-function(sfamily){
  n=dim(sfamily)[1]
  N<-dim(sfamily)[2]
  w<-array(0,dim=c(n,N,N))
  w2<-matrix(0,N,N)
  gg<-list()
  for (k in 1:n) {
    gg[[k]]<- as.undirected(graph_from_adjacency_matrix(sfamily[k,,]))
    E(gg[[k]])$weight<-edge_betweenness(gg[[k]],directed = FALSE)
    w[k,,]<-as.matrix(as_adjacency_matrix(gg[[k]],attr="weight"))
  }
  for(i in 1:N){
    for (j in 1:N) {
      if (sum(w[,i,j])==0){
        w2[i,j]=0
      }else{
        w2[i,j]<-sum(w[,i,j])/length(which(w[,i,j]!=0))
      }
    }
  }
  impor_vec<-as.vector(w2)/sum(w2)
  #impor_vec[which(impor_vec!=0)]<-1/impor_vec[which(impor_vec!=0)]
  return(impor_vec)
}

##find the important module for PA model
find_same_stru<-function(pa){
  n=dim(pa)[1]
  ll<-rep(0,n)
  index_cen<-which(degree(as.undirected(graph_from_adjacency_matrix(pa[1,,])))>=5)
  for (i in 1:n) {
    if (identical(pa[i,index_cen,index_cen],pa[1,index_cen,index_cen])){
      ll[i]=i
    }else{
      ll[i]=0
    }
  }
  return(ll)
}


find_same_stru_dc<-function(pa,m){
  n=dim(pa)[1]
  ll<-rep(0,n)
  index_cen<-which(degree(as.undirected(graph_from_adjacency_matrix(pa[1,,])))>=m)
  for (i in 1:n) {
    if (identical(pa[i,index_cen,index_cen],pa[1,index_cen,index_cen])){
      ll[i]=i
    }else{
      ll[i]=0
    }
  }
  return(ll)
}



##generate the PA family in simulation
PA_family2<-function(N,n,r){
  g1 <- sample_pa_age(N, pa.exp=1, aging.exp=0, aging.bin=1000)
  cen=as.matrix(get.adjacency(g1))
  array<-array(0,dim=c(n,N,N))
  array[1,,]<-as.matrix(cen)
  for (j in 2:n) {
    for (k in 1:r) {
      i<-sample(3:N,1)
      index<-which(cen[i,]==1)
      cen[i,index]=0
      index1<-sample(c(1:(i-1))[-index],1)
      cen[i,index1]=1
    }
    array[j,,]<-as.matrix(cen)
    cen=as.matrix(get.adjacency(g1))
  }
  return(array)
}


##simulations for the average rank of PA model (Figure 5, n=50 and n=100)
registerDoParallel(cl)
s <- Sys.time()
times=200
n=50
N=100
pro=0.1
r=ceiling(pro*N)
PAmean_rank_50 <- foreach(i=1:times,.combine='rbind') %dopar% {
  library(MASS)
  #library(nevada)
  library(igraph)
  pa<-PA_family2(N,n,r)
  index_same<-find_same_stru(pa)
  ww<-score(pa)
  ww_mat<-matrix(ww,N,N)
  ww1<-score_edge(pa)
  ww1_mat<-matrix(ww1,N,N)
  n1 <- ceiling(0.1*n)
  x<-2:(n1-1)
  num<-sample(x,1)
  pa[1:n1,,]<-Gen_Data.ER(N, n1, p=0.5, alpha= 0.1)
  for (i in 1:n1) {
    pa[i,,][upper.tri(pa[i,,])]<-0
  }
  vec_family<-vec_matrix(n,N,pa)
  eMDr<-rank(netEMD02(vec_family))
  weMDr<-rank(netEMD_score(vec_family,ww))
  weMDr1<-rank(netEMD_score(vec_family,ww1))
  Hamr<-rank(-min_hamming(n,pa)[[1]])
  wHamr<-rank(-min_hamming_w(n,pa,ww_mat)[[1]])
  wHamr1<-rank(-min_hamming_w(n,pa,ww1_mat)[[1]])
  eMD<-c(mean(eMDr[which(index_same!=0)]),mean(eMDr[which(index_same==0)]))
  weMD<-c(mean(weMDr[which(index_same!=0)]),mean(weMDr[which(index_same==0)]))
  weMD1<-c(mean(weMDr1[which(index_same!=0)]),mean(weMDr1[which(index_same==0)]))
  Ham<-c(mean( Hamr[which(index_same!=0)]),mean( Hamr[which(index_same==0)]))
  wHam<-c(mean( wHamr[which(index_same!=0)]),mean( wHamr[which(index_same==0)]))
  wHam1<-c(mean( wHamr1[which(index_same!=0)]),mean( wHamr1[which(index_same==0)]))
  sim<-cbind(Ham,wHam,wHam1,eMD,weMD,weMD1)
}
e <- Sys.time()
print(e-s)
library(reshape2)
mean_rank_same1<-PAmean_rank_50
colnames(mean_rank_same1)<-c("Ham_equal","Ham_node","Ham_edge","eMD_equal","eMD_node","eMD_edge")
rank_res_same_50<-melt(mean_rank_same1)
rank_res_same_50$group<-rep("same IM",6*2*times)
rank_res_same_50$group[seq(2,6*2*times,2)]="diff IM"
ar1<-ggplot(rank_res_same_50,aes(x=Var2,y=value),color=Var2)+
  theme_bw()+geom_boxplot(aes(fill=group))+ylim(0,45)+
  labs(x=" Type of Importance score",y="Average Rank",title="",fill="")+
  theme(legend.position=c(0.8,0.8))+
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(legend.title =element_blank())+
  theme(legend.key.size = unit(20, "pt"))+
  theme(axis.title=element_text(size=15),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
  scale_fill_manual(values=c("#4DAF4A","#984EA3"))
 # scale_color_manual(values=c("#4DAF4A","#984EA3"),limits = c("diff IM", "same IM"))


registerDoParallel(cl)
s <- Sys.time()
times=200
n=100
N=100
pro=0.1
r=ceiling(pro*N)
PAmean_rank_100 <- foreach(i=1:times,.combine='rbind') %dopar% {
  library(MASS)
  #library(nevada)
  library(igraph)
  pa<-PA_family2(N,n,r)
  index_same<-find_same_stru(pa)
  ww<-score(pa)
  ww_mat<-matrix(ww,N,N)
  ww1<-score_edge(pa)
  ww1_mat<-matrix(ww1,N,N)
  n1 <- ceiling(0.1*n)
  x<-2:(n1-1)
  num<-sample(x,1)
  #family_SBM[1:num,,]<-SBM_array(n=num,N,p_in=0.5,p_bt=0.2,p_in_per=0.05,p_bt_per=0.02)
  #family_SBM[1:num,,]<-Gen_Data.RGG(N, num, radius=0.175, sigma= 0.2)
  pa[1:n1,,]<-Gen_Data.ER(N, n1, p=0.5, alpha= 0.1)
  for (i in 1:n1) {
    pa[i,,][upper.tri(pa[i,,])]<-0
  }
  vec_family<-vec_matrix(n,N,pa)
  eMDr<-rank(netEMD02(vec_family))
  weMDr<-rank(netEMD_score(vec_family,ww))
  weMDr1<-rank(netEMD_score(vec_family,ww1))
  Hamr<-rank(-min_hamming(n,pa)[[1]])
  wHamr<-rank(-min_hamming_w(n,pa,ww_mat)[[1]])
  wHamr1<-rank(-min_hamming_w(n,pa,ww1_mat)[[1]])
  eMD<-c(mean(eMDr[which(index_same!=0)]),mean(eMDr[which(index_same==0)]))
  weMD<-c(mean(weMDr[which(index_same!=0)]),mean(weMDr[which(index_same==0)]))
  weMD1<-c(mean(weMDr1[which(index_same!=0)]),mean(weMDr1[which(index_same==0)]))
  Ham<-c(mean( Hamr[which(index_same!=0)]),mean( Hamr[which(index_same==0)]))
  wHam<-c(mean( wHamr[which(index_same!=0)]),mean( wHamr[which(index_same==0)]))
  wHam1<-c(mean( wHamr1[which(index_same!=0)]),mean( wHamr1[which(index_same==0)]))
  sim<-cbind(Ham,wHam,wHam1,eMD,weMD,weMD1)
}
e <- Sys.time()
print(e-s)
library(reshape2)
mean_rank_same1<-PAmean_rank_100
colnames(mean_rank_same1)<-c("Ham_equal","Ham_node","Ham_edge","eMD_equal","eMD_node","eMD_edge")
rank_res_same_100<-melt(mean_rank_same1)
rank_res_same_100$group<-rep("same IM",6*2*times)
rank_res_same_100$group[seq(2,6*2*times,2)]="diff IM"
ar2<-ggplot(rank_res_same_100,aes(x=Var2,y=value),color=Var2)+
  theme_bw()+geom_boxplot(aes(fill=group))+ylim(0,85)+
  labs(x=" Type of Importance score",y="Average Rank",title="",fill="")+
  theme(legend.position="none")+
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(legend.title =element_blank())+
  theme(legend.key.size = unit(20, "pt"))+
  theme(axis.title=element_text(size=12),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
  scale_fill_manual(values=c("#4DAF4A","#984EA3"))




####generate blog family in simulations
library(sand)
plot(aidsblog)
blog_family<-function(N,n,pro){
  cen=as.vector(get.adjacency(aidsblog))
  index<-which(cen==1)
  array<-array(0,dim=c(n,N,N))
  array[1,,]<-as.matrix(get.adjacency(aidsblog))
  for (i in 2:n) {
    index1<-sample(index,ceiling(pro*N))
    cen[index1]=0
    array[i,,]<-matrix(cen,N,N)
    cen[index1]=1
  }
  return(array)
}


##simulations for the average rank of blog model (Figure 5, n=50 and n=100)
registerDoParallel(cl)
s <- Sys.time()
times=200
n=50
pro=0.1
#r=ceiling(pro*N)
blogmean_rank_50 <- foreach(i=1:times,.combine='rbind') %dopar% {
  library(MASS)
  #library(nevada)
  library(igraph)
  library(sand)
  N<-length(V(aidsblog))
  pa<-blog_family(N,n,pro)
  index_same<-find_same_stru_dc(pa,20)
  ww<-score(pa)
  ww_mat<-matrix(ww,N,N)
  ww1<-score_edge(pa)
  ww1_mat<-matrix(ww1,N,N)
  n1 <- ceiling(0.1*n)
  x<-2:(n1-1)
  num<-sample(x,1)
  #family_SBM[1:num,,]<-SBM_array(n=num,N,p_in=0.5,p_bt=0.2,p_in_per=0.05,p_bt_per=0.02)
  #family_SBM[1:num,,]<-Gen_Data.RGG(N, num, radius=0.175, sigma= 0.2)
  pa[1:n1,,]<-Gen_Data.ER(N, n1, p=0.5, alpha= 0.1)
  for (i in 1:n1) {
     pa[i,,][upper.tri(pa[i,,])]<-0
    }
  vec_family<-vec_matrix(n,N,pa)
  eMDr<-rank(netEMD02(vec_family))
  weMDr<-rank(netEMD_score(vec_family,ww))
  weMDr1<-rank(netEMD_score(vec_family,ww1))
  Hamr<-rank(-min_hamming(n,pa)[[1]])
  wHamr<-rank(-min_hamming_w(n,pa,ww_mat)[[1]])
  wHamr1<-rank(-min_hamming_w(n,pa,ww1_mat)[[1]])
  eMD<-c(mean(eMDr[which(index_same!=0)]),mean(eMDr[which(index_same==0)]))
  weMD<-c(mean(weMDr[which(index_same!=0)]),mean(weMDr[which(index_same==0)]))
  weMD1<-c(mean(weMDr1[which(index_same!=0)]),mean(weMDr1[which(index_same==0)]))
  Ham<-c(mean( Hamr[which(index_same!=0)]),mean( Hamr[which(index_same==0)]))
  wHam<-c(mean( wHamr[which(index_same!=0)]),mean( wHamr[which(index_same==0)]))
  wHam1<-c(mean( wHamr1[which(index_same!=0)]),mean( wHamr1[which(index_same==0)]))
  sim<-cbind(Ham,wHam,wHam1,eMD,weMD,weMD1)
}
e <- Sys.time()
print(e-s)
library(reshape2)
mean_rank_same1<-blogmean_rank_50
colnames(mean_rank_same1)<-c("Ham_equal","Ham_node","Ham_edge","eMD_equal","eMD_node","eMD_edge")
rank_blog50<-melt(mean_rank_same1)
rank_blog50$group<-rep("same IM",6*2*times)
rank_blog50$group[seq(2,6*2*times,2)]="diff IM"
arblog1<-ggplot(rank_blog50,aes(x=Var2,y=value),color=Var2)+
  theme_bw()+geom_boxplot(aes(fill=group))+ylim(0,45)+
  labs(x=" Type of Importance score",y="Average Rank",title="",fill="")+
  theme(legend.position="none")+
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(legend.title =element_blank())+
  theme(legend.key.size = unit(20, "pt"))+
  theme(axis.title=element_text(size=15),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
  scale_fill_manual(values=c("#4DAF4A","#984EA3"))
# scale_color_manual(values=c("#4DAF4A","#984EA3"),limits = c("diff IM", "same IM"))



registerDoParallel(cl)
stopCluster(cl)
cl<-makeCluster(50)
s <- Sys.time()
times=200
n=100
pro=0.1
#r=ceiling(pro*N)
blogmean_rank_100 <- foreach(i=1:times,.combine='rbind') %dopar% {
  library(MASS)
  #library(nevada)
  library(igraph)
  library(sand)
  N<-length(V(aidsblog))
  pa<-blog_family(N,n,pro)
  index_same<-find_same_stru_dc(pa,20)
  ww<-score(pa)
  ww_mat<-matrix(ww,N,N)
  ww1<-score_edge(pa)
  ww1_mat<-matrix(ww1,N,N)
  n1 <- ceiling(0.1*n)
  x<-2:(n1-1)
  num<-sample(x,1)
  #family_SBM[1:num,,]<-SBM_array(n=num,N,p_in=0.5,p_bt=0.2,p_in_per=0.05,p_bt_per=0.02)
  #family_SBM[1:num,,]<-Gen_Data.RGG(N, num, radius=0.175, sigma= 0.2)
  pa[1:n1,,]<-Gen_Data.ER(N, n1, p=0.5, alpha= 0.1)
  for (i in 1:n1) {
    pa[i,,][upper.tri(pa[i,,])]<-0
  }
  vec_family<-vec_matrix(n,N,pa)
  eMDr<-rank(netEMD02(vec_family))
  weMDr<-rank(netEMD_score(vec_family,ww))
  weMDr1<-rank(netEMD_score(vec_family,ww1))
  Hamr<-rank(-min_hamming(n,pa)[[1]])
  wHamr<-rank(-min_hamming_w(n,pa,ww_mat)[[1]])
  wHamr1<-rank(-min_hamming_w(n,pa,ww1_mat)[[1]])
  eMD<-c(mean(eMDr[which(index_same!=0)]),mean(eMDr[which(index_same==0)]))
  weMD<-c(mean(weMDr[which(index_same!=0)]),mean(weMDr[which(index_same==0)]))
  weMD1<-c(mean(weMDr1[which(index_same!=0)]),mean(weMDr1[which(index_same==0)]))
  Ham<-c(mean( Hamr[which(index_same!=0)]),mean( Hamr[which(index_same==0)]))
  wHam<-c(mean( wHamr[which(index_same!=0)]),mean( wHamr[which(index_same==0)]))
  wHam1<-c(mean( wHamr1[which(index_same!=0)]),mean( wHamr1[which(index_same==0)]))
  sim<-cbind(Ham,wHam,wHam1,eMD,weMD,weMD1)
}

e <- Sys.time()
print(e-s)
library(reshape2)
mean_rank_same1<-blogmean_rank_100
colnames(mean_rank_same1)<-c("Ham_equal","Ham_node","Ham_edge","eMD_equal","eMD_node","eMD_edge")
rank_blog100<-melt(mean_rank_same1)
rank_blog100$group<-rep("same IM",6*2*times)
rank_blog100$group[seq(2,6*2*times,2)]="diff IM"
arblog2<-ggplot(rank_blog100,aes(x=Var2,y=value),color=Var2)+
  theme_bw()+geom_boxplot(aes(fill=group))+ylim(0,85)+
  labs(x=" Type of Importance score",y="Average Rank",title="",fill="")+
  theme(legend.position="none")+
  ggtitle("")+theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(legend.title =element_blank())+
  theme(legend.key.size = unit(20, "pt"))+
  theme(axis.title=element_text(size=15),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
  scale_fill_manual(values=c("#4DAF4A","#984EA3"))
# scale_color_manual(values=c("#4DAF4A","#984EA3"),limits = c("diff IM", "same IM"))





##simulations for estimating centers for PA and blog model ( results of Table 2)
times=200
n=100
N=100
pro=0.1
r=ceiling(pro*N)
PA_cen<- foreach(i=1:times,.combine='rbind') %dopar% {
  library(randnet)
  library(igraph)
  library(sand)
  blog<-PA_family2(N,n,r)
  real_center<-blog[1,,]
  #ww<-score_edge(family_SBM)
  #ww_mat<-matrix(ww,N,N)
  #ww1<-score(family_SBM)
  #ww1_mat<-matrix(ww1,N,N)
  n1 <- ceiling(0.1*n)
  x<-2:(n1-1)
  num<-sample(x,1)
  blog[1:n1,,]<-Gen_Data.ER(N, n1, p=0.5, alpha= 0.1)
  for (i in 1:n1) {
    blog[i,,][upper.tri(blog[i,,])]<-0
  }
  ww<-score_edge(blog)
  ww_mat<-matrix(ww,N,N)
  ww1<-score(blog)
  ww1_mat<-matrix(ww1,N,N)
  Fre<-min_hamming(n,blog)
  Fre_w<-min_hamming_w(n,blog,ww_mat)
  Fre_w1<-min_hamming_w(n,blog,ww1_mat)
  vec<-vec_matrix(n,N,blog)
  EMD<-netEMD02(vec)
  EMD_w<-netEMD_score(vec,ww)
  EMD_w1<-netEMD_score(vec,ww1)
  index1<-Fre[[2]][1]
  index5<-Fre_w[[2]][1]
  index6<-Fre_w1[[2]][1]
  index2<-which(EMD==max(EMD))[1]
  index3<-which(EMD_w==max(EMD_w))[1]
  index4<-which(EMD_w1==max(EMD_w1))[1]
  ham<-c(dis_hamming_w(blog[index1,,],real_center,ww_mat),dis_hamming_w(blog[index5,,],real_center,ww_mat))
  EMD<-c(dis_hamming_w(blog[index2,,],real_center,ww_mat),dis_hamming_w(blog[index3,,],real_center,ww_mat))
  ham2<-c(dis_hamming_w(blog[index1,,],real_center,ww1_mat),dis_hamming_w(blog[index6,,],real_center,ww1_mat))
  EMD2<-c(dis_hamming_w(blog[index2,,],real_center,ww1_mat),dis_hamming_w(blog[index4,,],real_center,ww1_mat))
  sim<-cbind(ham,EMD,ham2,EMD2)
}

PA_cenres<-PA_cen
colnames(PA_cenres)<-c("Ham","eMD","Ham","eMD")
plot_para_edge<-melt(PA_cenres[,c(1:2)])
plot_para_edge$group<-rep("1equal",2*2*times)
plot_para_edge$group[seq(2,2*2*times,2)]="2edge"

pa_edge<-ggplot(plot_para_edge,aes(x=Var2,y=value),color=Var2)+theme_bw()+
  geom_boxplot(aes(fill=group))+labs(x="",y="EA(wEA)",title="",fill="")+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(legend.title =element_blank())+
  theme(legend.key.size = unit(20, "pt"))+
  theme(axis.title=element_text(size=15),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
  scale_fill_manual(values=c("#4876FF","#EEB422"),limits = c("1equal","2edge"))
#scale_color_manual(values=c("#4DAF4A","#984EA3"),limits = c("equal","edge"))

###

mean(plot_para_edge[which(plot_para_edge$group[1:400]=="1equal"),]$value)
sd(plot_para_edge[which(plot_para_edge$group[1:400]=="1equal"),]$value)

mean(plot_para_edge[which(plot_para_edge$group[1:400]=="2edge"),]$value)
sd(plot_para_edge[which(plot_para_edge$group[1:400]=="2edge"),]$value)

mean(plot_para_edge[which(plot_para_edge$group[401:800]=="1equal")+400,]$value)
sd(plot_para_edge[which(plot_para_edge$group[401:800]=="1equal")+400,]$value)

mean(plot_para_edge[which(plot_para_edge$group[401:800]=="2edge")+400,]$value)
sd(plot_para_edge[which(plot_para_edge$group[401:800]=="2edge")+400,]$value)



plot_para_node<-melt(PA_cenres[,c(3:4)])
plot_para_node$group<-rep("equal IS",2*2*times)
plot_para_node$group[seq(2,2*2*times,2)]="IS node"

pa_node<-ggplot(plot_para_node,aes(x=Var2,y=value),color=Var2)+theme_bw()+
  geom_boxplot(aes(fill=group))+labs(x="",y="EA(wEA)",title="",fill="")+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size=12))+
  theme(legend.title =element_blank())+
  theme(legend.key.size = unit(20, "pt"))+
  theme(axis.title=element_text(size=15),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
  scale_fill_manual(values=c("#4876FF","#EEB422"))

##
mean(plot_para_node[which(plot_para_node$group[1:400]=="equal IS"),]$value)
sd(plot_para_node[which(plot_para_node$group[1:400]=="equal IS"),]$value)

mean(plot_para_node[which(plot_para_node$group[1:400]=="IS node"),]$value)
sd(plot_para_node[which(plot_para_node$group[1:400]=="IS node"),]$value)

mean(plot_para_node[which(plot_para_node$group[401:800]=="equal IS")+400,]$value)
sd(plot_para_node[which(plot_para_node$group[401:800]=="equal IS")+400,]$value)


mean(plot_para_node[which(plot_para_node$group[401:800]=="IS node")+400,]$value)
sd(plot_para_node[which(plot_para_node$group[401:800]=="IS node")+400,]$value)


times=200
n=100
pro=0.1
r=ceiling(pro*N)
 blog_cen<- foreach(i=1:times,.combine='rbind') %dopar% {
  library(randnet)
  library(igraph)
   library(sand)
   N<-length(V(aidsblog))
  blog<-blog_family(N,n,pro)
  real_center<-blog[1,,]
  n1 <- ceiling(0.1*n)
  x<-2:(n1-1)
  num<-sample(x,1)
  blog[1:n1,,]<-Gen_Data.ER(N, n1, p=0.5, alpha= 0.1)
  for (i in 1:n1) {
    blog[i,,][upper.tri(blog[i,,])]<-0
  }
  ww<-score_edge(blog)
  ww_mat<-matrix(ww,N,N)
  ww1<-score(blog)
  ww1_mat<-matrix(ww1,N,N)
  Fre<-min_hamming(n,blog)
  Fre_w<-min_hamming_w(n,blog,ww_mat)
  Fre_w1<-min_hamming_w(n,blog,ww1_mat)
  vec<-vec_matrix(n,N,blog)
  EMD<-netEMD02(vec)
  EMD_w<-netEMD_score(vec,ww)
  EMD_w1<-netEMD_score(vec,ww1)
  #sd<-SD(family_SBM)
  index1<-Fre[[2]][1]
  index5<-Fre_w[[2]][1]
  index6<-Fre_w1[[2]][1]
  index2<-which(EMD==max(EMD))[1]
  index3<-which(EMD_w==max(EMD_w))[1]
  index4<-which(EMD_w1==max(EMD_w1))[1]
  ham<-c(dis_hamming_w(blog[index1,,],real_center,ww_mat),dis_hamming_w(blog[index5,,],real_center,ww_mat))
  EMD<-c(dis_hamming_w(blog[index2,,],real_center,ww_mat),dis_hamming_w(blog[index3,,],real_center,ww_mat))
  #ham_edge<-dis_hamming_w(blog[index5,,],real_center,ww_mat)
  #EMD_edge<-dis_hamming_w(blog[index3,,],real_center,ww_mat)
  ham2<-c(dis_hamming_w(blog[index1,,],real_center,ww1_mat),dis_hamming_w(blog[index6,,],real_center,ww1_mat))
  #ham_node<-dis_hamming_w(blog[index6,,],real_center,ww1_mat)
  EMD2<-c(dis_hamming_w(blog[index2,,],real_center,ww1_mat),dis_hamming_w(blog[index4,,],real_center,ww1_mat))
  #EMD_node<-dis_hamming_w(blog[index4,,],real_center,ww1_mat)
  sim<-cbind(ham,EMD,ham2,EMD2)
}

 blog_cenres<-blog_cen
 colnames(blog_cenres)<-c("Ham","eMD","Ham","eMD")
blog_edge_res<-melt(blog_cenres[,c(1:2)])
 blog_edge_res$group<-rep("equal IS",2*2*times)
blog_edge_res$group[seq(2,2*2*times,2)]="IS edge"

 blog_edge<-ggplot(blog_edge_res,aes(x=Var2,y=value),color=Var2)+theme_bw()+
   geom_boxplot(aes(fill=group))+labs(x="",y="EA(wEA)",title="",fill="")+
   theme(legend.position="none")+
   theme(plot.title = element_text(hjust = 0.5,size=12))+
   theme(legend.title =element_blank())+
   theme(legend.key.size = unit(20, "pt"))+
   theme(axis.title=element_text(size=15),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
   scale_fill_manual(values=c("#4876FF","#EEB422"))


 ##
 mean(blog_edge_res[which(blog_edge_res$group[1:400]=="equal IS"),]$value)
 sd(blog_edge_res[which(blog_edge_res$group[1:400]=="equal IS"),]$value)

 mean(blog_edge_res[which(blog_edge_res$group[1:400]=="IS edge"),]$value)
 sd(blog_edge_res[which(blog_edge_res$group[1:400]=="IS edge"),]$value)

 mean(blog_edge_res[which(blog_edge_res$group[401:800]=="equal IS")+400,]$value)
 sd(blog_edge_res[which(blog_edge_res$group[401:800]=="equal IS")+400,]$value)


 mean(blog_edge_res[which(blog_edge_res$group[401:800]=="IS edge")+400,]$value)
 sd(blog_edge_res[which(blog_edge_res$group[401:800]=="IS edge")+400,]$value)




 blog_node_res<-melt(blog_cenres[,c(3:4)])
 blog_node_res$group<-rep("equal IS",2*2*times)
 blog_node_res$group[seq(2,2*2*times,2)]="IS node"

 blog_node<-ggplot(blog_node_res,aes(x=Var2,y=value),color=Var2)+theme_bw()+
   geom_boxplot(aes(fill=group))+labs(x="",y="EA(wEA)",title="",fill="")+
   theme(legend.position="none")+
   theme(plot.title = element_text(hjust = 0.5,size=12))+
   theme(legend.title =element_blank())+
   theme(legend.key.size = unit(20, "pt"))+
   theme(axis.title=element_text(size=15),axis.title.x =element_text(size=15), axis.title.y=element_text(size=15),axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
   scale_fill_manual(values=c("#4876FF","#EEB422"))



 ##
 mean(blog_node_res[which(blog_node_res$group[1:400]=="equal IS"),]$value)
 sd(blog_node_res[which(blog_node_res$group[1:400]=="equal IS"),]$value)

 mean(blog_node_res[which(blog_node_res$group[1:400]=="IS node"),]$value)
 sd(blog_node_res[which(blog_node_res$group[1:400]=="IS node"),]$value)

 mean(blog_node_res[which(blog_node_res$group[401:800]=="equal IS")+400,]$value)
 sd(blog_node_res[which(blog_node_res$group[401:800]=="equal IS")+400,]$value)


 mean(blog_node_res[which(blog_node_res$group[401:800]=="IS node")+400,]$value)
 sd(blog_node_res[which(blog_node_res$group[401:800]=="IS node")+400,]$value)

