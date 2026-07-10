#'
#' This function calculates the Edgewise Multiplex Depth for network data.
#'
#'
#'
#' @param data an n \* N \* N array, the ensember of network, where n is the number of networks, and N is the number of vertex of each network;
#' @param K an integer, the parameter for depth.
#' @param unweighted a bool, whether the networks are unweighted or weighted.
#' @return A vector of the estimated Edgewise Multiplex Depth for the n network samples.
#' @export
net_EMD=function(data,K=2,unweighted=F){
  n=dim(data)[1]
  N=dim(data)[2]
  vec_mat<-matrix(0,N^2,n)
  for (i in 1:n) {
    vec_mat[,i]<-as.vector(data[i,,])
  }
  p=N^2

  # If we use unweighted matrix, we can use our fast algorithm
  if(unweighted){
    # calculate how many samples are zero/one for each edge
    num.one <- rowSums(vec_mat)
    num.zero <- n - num.one
    # Calculate depth according to our fast algorithm
    depth_val = rep(0,n)
    for(k in 2:K){
      # The number of combination that choosed edge are all one
      choose.allone = choose(num.one,k)
      # The number of combination that choosed edge are all zero
      choose.allzero = choose(num.zero,k)
      # (1-vec_mat)*choose.allone represents:
      depth_val = depth_val+(colSums(choose(n,k)-(1-vec_mat)*choose.allone-vec_mat*choose.allzero))/choose(n,k)/p
    }

  }else{
    # If we use weighted matrix, we can use rank to accelarate our algorithm
    rmat=apply(vec_mat,1,rank)
    down=rmat-1
    up=n-rmat
    depth_val = rep(0,n)
    for(k in 2:K){
      depth_val = depth_val+(rowSums(choose(n,k)-choose(up,k)-choose(down,k)))/choose(n,k)/p
    }
  }
  return(depth_val)
}



#' This function calculates the Edgewise Multiplex Depth for network data with importance score matrix.
#'
#'
#'
#' @param data an n \* N \* N array, the ensember of network, where n is the number of networks, and N is the number of vertex of each network;
#' @param w an N\*N importance score matrix, corresponding to the weight of each edge in the network samples.
#' @param K an integer, the parameter for depth.
#' @param unweighted a bool, whether the networks are unweighted or weighted.
#' @return A vector of the estimated Edgewise Multiplex Depth for the n network samples.
#' @export
net_EMD.IS=function(data,w,K=2,unweighted=F){
  n=dim(data)[1]
  N=dim(data)[2]
  vec_mat<-matrix(0,N^2,n)

  w = w/sum(w) #normalize
  w<-as.vector(w)
  for (i in 1:n) {
    vec_mat[,i]<-as.vector(data[i,,])
  }
  p=N^2
  # If we use unweighted matrix, we can use our fast algorithm
  if(unweighted){
    # calculate how many samples are zero/one for each edge
    num.one <- rowSums(vec_mat)
    num.zero <- n - num.one
    # Calculate depth according to our fast algorithm
    depth_val = rep(0,n)
    for(k in 2:K){
      # The number of combination that choosed edge are all one
      choose.allone = choose(num.one,k)
      # The number of combination that choosed edge are all zero
      choose.allzero = choose(num.zero,k)
      # (1-vec_mat)*choose.allone represents:
      depth_val = depth_val+(t(choose(n,k)-(1-vec_mat)*choose.allone-vec_mat*choose.allzero)%*%w)/choose(n,k)
    }

  }else{
    # If we use weighted matrix, we can use rank to accelarate our algorithm
  rmat=apply(vec_mat,1,rank)
  down=rmat-1
  up=n-rmat
  depth_val = rep(0,n)
  for(k in 2:K){
    depth_val = depth_val+((choose(n,k)-choose(up,k)-choose(down,k))%*%w)/choose(n,k)
  }
  }
  depth_val = as.vector(depth_val)
  return(depth_val)
}

