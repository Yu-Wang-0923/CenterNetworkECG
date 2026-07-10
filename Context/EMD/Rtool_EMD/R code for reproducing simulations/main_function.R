# necessary packages
library(doParallel) 
library(parallel)
library(reshape2)
library(ggplot2)
####Function for EMD 
netEMD=function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=rmat-1
  up=n-rmat
  (rowSums(up*down)/p+n-1)/choose(n,2)
}


## convert the network array to a matrix
vec_matrix<-function(n,N,array){
  vec_mat<-matrix(0,N^2,n)
  for (i in 1:n) {
    vec_mat[,i]<-as.vector(array[i,,])
  }
  return(vec_mat)
}


SD=function(data)
{
  DIM=dim(data)
  sd=seq(DIM[1])
  for (i in 1:DIM[1])
  { temp=NULL
  for (k in 1:(DIM[1]-1))
  {
    for (j in (k+1):DIM[1])
    {
      for (l in 1:DIM[2])
      {
        #temp=apply(data[i,,],1,function(x) (min(x-data[k,i,]*data[j,i,])>-1)*(min(data[k,i,]+data[j,i,]-x)>-1))
        temp=c(temp,(min(data[i,l,]-data[k,l,]*data[j,l,])>-1)*(min(data[k,l,]+data[j,l,]-data[i,l,])>-1))
      }
    }
  }
  sd[i]=mean(temp)
  }
  return(sd)
}



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  #datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##Hamming distance
dis_hamming<-function(x,y){
  dis<-sum(abs(x-y))/(nrow(x)*(nrow(x)-1))
  return(dis)
}


##frobenius distance
dis_frobenius<-function(x,y){
  dis<-sqrt(sum((x-y)^2))
  return(dis)
}

##weighted Hamming distance
dis_hamming_w<-function(x,y,w){
  dis<-sum(abs(x-y)*w)
  return(dis)
}


##Fréchet mean with Hamming distance
min_hamming <- function(n,weight_family){
  dist<-matrix(0,n,n)
  dist_m<-rep(0,n)
  for (i in 1:n) {
    for (j in 1:n) {
      dist[i,j]<-dis_hamming(weight_family[i,,],weight_family[j,,])^2
    }
    dist_m[i] <- sum(dist[i,])/(n-1)
  }
  min_in <- which(dist_m==min(dist_m))
  max_in <- which(dist_m==max(dist_m))
  return(list(dist_m,min_in,max_in))
}

##Fréchet mean with weighted Hamming distance
min_hamming_w <- function(n,weight_family,w){
  dist<-matrix(0,n,n)
  dist_m<-rep(0,n)
  for (i in 1:n) {
    for (j in 1:n) {
      dist[i,j]<-dis_hamming_w(weight_family[i,,],weight_family[j,,],w)^2
    }
    dist_m[i] <- sum(dist[i,])/(n-1)
  }
  min_in <- which(dist_m==min(dist_m))
  max_in <- which(dist_m==max(dist_m))
  return(list(dist_m,min_in,max_in))
}



##Fréchet mean with frobenius distance
min_froben <- function(n,weight_family){
  dist<-matrix(0,n,n)
  dist_m<-rep(0,n)
  for (i in 1:n) {
    for (j in 1:n) {
      dist[i,j]<-dis_frobenius(weight_family[i,,],weight_family[j,,])^2
    }
    dist_m[i] <- sum(dist[i,])/(n-1)
  }
  min_in <- which(dist_m==min(dist_m))
  max_in <- which(dist_m==max(dist_m))
  return(list(dist_m,min_in,max_in))
}



