

cd(./Rtool_EMD/R_package)  ########change this to the path of the downloaded file
devtools::install("netEMD")
library(netEMD)
help(net_EMD)
help(net_EMD.IS)

n=100
N=50
###generate the network ensember of n networks.
data<-array(rbinom(n*N*N,1,0.1),c(n,N,N))

###generate the importance score matrix
weight<-runif(N*N,min=0,max=1)/(N*N)
#####calculate the EMD with equal importance scores
EMD<-net_EMD(data,K=2,unweighted=T)
EMD_3<-net_EMD(data,K=3,unweighted=T)
#####calculate the EMD with importance score matrix "weight"
EMD_IS<-net_EMD.IS(data,weight,K=2,unweighted=T)
EMD_IS_3<-net_EMD.IS(data,weight,K=3,unweighted=T)



