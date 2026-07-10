n=200
N=30
#cl <- makeCluster(50)
source("main_function.R", echo=F)
source("Gen_weighted_network.R", echo=F)
registerDoParallel(cl)
s <- Sys.time()
averange_rank_weighted <- foreach(i=1:200,.combine='rbind') %dopar% {
  library(MASS)
  library(nevada)
  family<-WRG_family(n,N)
  family_out<-family
  x<-sample(2:18,1)
  family_out[1:x,,]<-weight_out_family(x,N)
  family_out[(x+1):20,,]<-weight_out1_family(20-x,N)
  Fre<-min_hamming(n,family_out)
  froben<-min_froben(n,family_out)
  vec<-vec_matrix(n,N,family_out)
  mbd<-netEMD(vec)
  Ham<-mean(rank(-Fre[[1]])[1:20])
  Frob<-mean(rank(-froben[[1]])[1:20])
  eMD<-mean(rank(mbd)[1:20])
  sim<-cbind(Ham,Frob,eMD)
}
e <- Sys.time()
print(e-s)

colMeans(averange_rank_weighted)
sd(averange_rank_weighted[,1])
sd(averange_rank_weighted[,2])
sd(averange_rank_weighted[,3])


n=200
N=50
#cl <- makeCluster(50)
registerDoParallel(cl)
s <- Sys.time()
averange_rank_weighted50 <- foreach(i=1:200,.combine='rbind') %dopar% {
  library(MASS)
  library(nevada)
  family<-WRG_family(n,N)
  family_out<-family
  x<-sample(2:18,1)
  family_out[1:x,,]<-weight_out_family(x,N)
  family_out[(x+1):20,,]<-weight_out1_family(20-x,N)
  Fre<-min_hamming(n,family_out)
  froben<-min_froben(n,family_out)
  vec<-vec_matrix(n,N,family_out)
  mbd<-netEMD(vec)
  Ham<-mean(rank(-Fre[[1]])[1:20])
  Frob<-mean(rank(-froben[[1]])[1:20])
  eMD<-mean(rank(mbd)[1:20])
  sim<-cbind(Ham,Frob,eMD)
}
e <- Sys.time()
print(e-s)
colMeans(averange_rank_weighted50)
sd(averange_rank_weighted50[,1])
sd(averange_rank_weighted50[,2])
sd(averange_rank_weighted50[,3])
