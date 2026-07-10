#==== Data Generation
library(igraph)

#=== ER Model
#' @param N The number of vertex;
#' @param n The number of network;
#' @param p The probability of generating the central network of each edge;
#' @param alpha The probability of perturbation
#' @param verbose a logical value indicating whether to print out progress messages during computation.
#'
#' @return A n \* N \* N array, the ensember of network. The first network is the true center.
#' This function returns a directed network ensemble; self loop is allowed here
#' Test: ER.ensemble <- Gen_Data.ER(N=3, n=4, p=0.1, alpha= 0.01)

Gen_Data.ER <- function(N=10, n=10, p=0.5, alpha= 0.1){
  adj_center.ER <- matrix(rbinom(N*N,1,p),N,N)
  #--- Or simply use
  #G.ensemble.Ham_Pur <- Gen_Data.Ham_Per(adj_center.ER, n=n, alpha= alpha)
  #return(Network.Ensemble=G.ensemble.Ham_Pur)

  Mm_family.ER <- array(rbinom(n*N*N,1,alpha), dim = c(n, N, N))
  Mm_family.ER[1,,] <- matrix(0,N,N)
  Mm_family.ER <- t(apply(Mm_family.ER, 1, function(X){return(xor(X,adj_center.ER))}))

  dim(Mm_family.ER) <- c(n,N,N)

  return(Network.Ensemble=Mm_family.ER)
}

#=== RGG Model
#' @param N The number of vertex;
#' @param n The number of network;
#' @param radius The distance of connection;
#' @param sigma The noise level;
#' @param verbose a logical value indicating whether to print out progress messages during computation.
#'
#' @return A n \* N \* N array, the ensember of network. The first network is the true center.
#' This function returns an undirected network ensemble
#' Test: RGG.ensemble <- Gen_Data.RGG(N=50, n=3, radius=0.175, sigma= 0.1)

Gen_Data.RGG <- function(N=10, n=10, radius=0.175, sigma= 0.1){
  p.loc <- matrix(runif(N*2),N,2)

  p_family.loc <- array(rnorm(n*N*2,0,sigma), dim = c(n, N, 2))

  p_family.loc[1,,] <- matrix(0,N,2)

  #We can use trunc to contrain our point into [0,1]X[0,1]
  #p_family.loc <- t(apply(p_family.loc, 1, function(X){return(trunc(X+p.loc,0,1))}))
  p_family.loc <- t(apply(p_family.loc, 1, function(X){return(X+p.loc)}))

  # Once touch the boundary, mirrors back
  while(min(c(p_family.loc))<0 | max(c(p_family.loc))>1){
    less_0 <- which(p_family.loc<0)
    plus_1 <- which(p_family.loc>1)
    p_family.loc[less_0] <- -p_family.loc[less_0]
    p_family.loc[plus_1] <- 1-(p_family.loc[plus_1]-1)
  }

  dim(p_family.loc) <- c(n,N,2)

  Connect <- function(p.loc){
    Dis_1 <- (t(matrix(rep(p.loc[,1],N),N,N))-matrix(rep(p.loc[,1],N),N,N))^2
    Dis_2 <- (t(matrix(rep(p.loc[,2],N),N,N))-matrix(rep(p.loc[,2],N),N,N))^2

    G <- sqrt(Dis_1+Dis_2)<radius
    #Remove self-loop
    diag(G) <- F
    return(G)
  }

  Mm_family.RGG <- t(apply(p_family.loc, 1, Connect))
  dim(Mm_family.RGG) <- c(n,N,N)

  return(Network.Ensemble=Mm_family.RGG)
}


#=== SBM Model
#' @param N The number of vertex;
#' @param n The number of network;
#' @param k The number of blocks;
#' @param p.block The probability of distributing vertexes into each block.
#' @param p.in The probability of linking an edge within the same group;
#' @param p.bt The probability of linking an edge between group;
#' @param p.in.per The probability of perturbing edges within the same group;
#' @param p.bt.per The probability of perturbing edges between group;
#' @param verbose a logical value indicating whether to print out progress messages during computation.
#'
#' @return A n \* N \* N array, the ensember of network. The first network is the true center.
#' There are several possible types of perturbations in SBM model, we only consider the reconnected edge.
#' We can also consider the block varies.
#' This function returns an undirected network ensemble
#' When: p.in=1, p.bt=1, p.in.per=1, p.bt.per=0. The network generated is totally blocked.
#' Test: SBM.ensemble <- Gen_Data.SBM(N=50, n=3,k=3, p.block=NULL, p.in=0.16, p.bt=0.075, p.in.per=0.05, p.bt.per=0.02)

Gen_Data.SBM <- function(N=10, n=10, k=3, p.block=NULL,
                         p.in=0.16, p.bt=0.075,
                         p.in.per=0.05, p.bt.per=0.02){
  # Calculate the size of each block
  if(is.null(p.block)) p.block <- rep(1/k,k)
  size.block <- rmultinom(1,N,p.block)

  # Generate pref.matrix for center
  pm <- matrix(p.bt,k,k)
  diag(pm) <- p.in

  adj_center.SBM <- as.matrix(as_adjacency_matrix(sample_sbm(N, pref.matrix=pm, block.sizes=size.block)))

  # Generate pref.matrix for perturbation
  pm.per <- matrix(p.bt.per,k,k)
  diag(pm.per) <- p.in.per
  Perturb <- function(G){
    G <- sample_sbm(N, pref.matrix=pm.per, block.sizes=size.block)
    return(as.matrix(as_adjacency_matrix(G)))
  }
  Mm_family.SBM <- t(sapply(c(1:n), Perturb))
  dim(Mm_family.SBM) <- c(n,N,N)
  Mm_family.SBM[1,,] <- matrix(0,N,N)

  # Add perturbation
  Mm_family.SBM <- t(apply(Mm_family.SBM, 1, function(X){return(xor(X,adj_center.SBM))}))
  dim(Mm_family.SBM) <- c(n,N,N)

  return(Network.Ensemble=Mm_family.SBM)
}


#=== SW Model
#' @param N The number of vertex;
#' @param n The number of network;
#' @param d The dimension of lattice
#' @param p The probability of rewiring(In center network);
#' @param alpha The probability of rewiring from the center network;
#' @param verbose a logical value indicating whether to print out progress messages during computation.
#'
#' @return A n \* N \* N array, the ensember of network. The first network is the true center.
#' This function returns an undirected network ensemble;without loops and multi-edges.
#' Test: SW.ensemble <- Gen_Data.SW(N=50, n=3, d=2, p=0.2, alpha = 0.1)

Gen_Data.SW <- function(N=10, n=10, d=2, p=0.2, alpha = 0.1){
  g.center <- sample_smallworld(1, N, d, p)
  adj_center.SW <- as.matrix(as_adjacency_matrix(g.center))

  Rewire <- function(i){
    G <- g.center %>% rewire(each_edge(p = alpha, loops = FALSE,multiple = FALSE))  #Remove self-loop
    return(as.matrix(as_adjacency_matrix(G)))
  }

  Mm_family.SW <- t(sapply(c(1:n), Rewire))
  dim(Mm_family.SW) <- c(n,N,N)

  Mm_family.SW[1,,] <- adj_center.SW
  # Visualization
  # G <- as.undirected(graph.adjacency(Mm_family.SW[1,,]))
  # plot(G)
  return(Network.Ensemble=Mm_family.SW)
}

#=== Hamming Distance Mask
#' This function is going to add noise via the hamming distance between the target graph and the center graph.
#' @param G The center graph(N * N matrix);
#' @param n The number of network;
#' @param alpha The probability of perturbation;
#' @param verbose a logical value indicating whether to print out progress messages during computation.
#'
#' @return A n \* N \* N array, the ensember of network. The first network is the true center.
#' Test:
#' N <- 50; p <- 0.5; G <- matrix(rbinom(N*N,1,p),N,N)
#' G.ensemble.Ham_Pur <- Gen_Data.Ham_Per(G, n=10, alpha= 0.1)

Gen_Data.Ham_Per <- function(G=NULL, n=10, alpha= 0.1){
  if(is.null(G)){
    cat("There is no center network.\n")
  }

  if(is.matrix(G)){
    if(dim(G)[1]!=dim(G)[2]){
      cat("The adjacent matric of the center network should be a square matrix.\n")
    }
    else{
      N <- dim(G)[1]
    }
  }else{
    cat("The adjacent matric of the center network should be a matrix.\n")
  }

  Mm_family <- array(rbinom(n*N*N,1,alpha), dim = c(n, N, N))
  Mm_family[1,,] <- matrix(0,N,N)
  Mm_family <- t(apply(Mm_family, 1, function(X){return(xor(X,G))}))

  dim(Mm_family) <- c(n,N,N)

  return(Network.Ensemble=Mm_family)
}

#==== The unified approach of generate data
Gen_Data <- function(N,n,type=c('ER','RGG','SBM','SW'),
                     p=0.1, alpha= 0.01, #For ER
                     radius=0.175, sigma= 0.1, #For RGG
                     k=3, p.block=NULL, p.in=0.16, p.bt=0.075, p.in.per=0.05, p.bt.per=0.02, #For SBM
                     d=2 #For SW p=0.2, alpha = 0.1
                     ){
  type <- match.arg(type)
  #==== Generate Data according to our type
  if (type == 'ER') {
    Network.ensemble <- Gen_Data.ER(N, n, p, alpha)
  }

  if (type == 'RGG') {
    Network.ensemble <- Gen_Data.RGG(N, n, radius, sigma)
  }

  if (type == 'SBM') {
    Network.ensemble <- Gen_Data.SBM(N, n, k, p.block,
                                     p.in, p.bt, p.in.per, p.bt.per)
  }

  if (type == 'SW') {
    Network.ensemble <- Gen_Data.SW(N, n, d, p, alpha)
  }

  return(Network.ensemble)
}


#' Generate Mixture Data
#' @param mix.prop The proportion of mixture

Gen_Mix_Data <- function(N,n,type=c('ER','RGG','SBM','SW'),
                         p=0.1, alpha= 0.01, #For ER
                         radius=0.175, sigma= 0.1, #For RGG
                         k=3, p.block=NULL, p.in=0.16, p.bt=0.075, p.in.per=0.05, p.bt.per=0.02, #For SBM
                         d=2, #For SW p=0.2, alpha = 0.1
                         mix.prop = c(0.8,0.2)){
  #--- Initialize
  mix.prop <- sort(mix.prop, decreasing = TRUE)
  subsize <- round(n * mix.prop)
  if(sum(subsize)!=n) return("The proportion is not reasonable.")

  Mm_family_Mix <- array(dim=c(n,N,N))
  center_Mix <- array(dim=c(length(subsize),N,N))
  subindex <- cumsum(c(0,subsize))

  #--- Generate them separately
  for(i in 1:length(subsize)){
    Mm_family <- Gen_Data(N=N, n=subsize[i], type=type,
                          p=p, alpha=alpha, #For ER
                          radius=radius, sigma=sigma, #For RGG
                          k=k, p.block=p.block, p.in=p.in, p.bt=p.bt, p.in.per=p.in.per, p.bt.per=p.bt.per, #For SBM
                          d=d)
    Mm_family_Mix[(subindex[i]+1):subindex[i+1],,] <- Mm_family
    center_Mix[i,,] <- Mm_family[1,,]
  }

  return(list(Network.Ensemble = Mm_family_Mix,
              Network.center.ind = subindex[1:(length(subsize))]+1))
}

