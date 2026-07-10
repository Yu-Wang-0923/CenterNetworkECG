
###generate the matrix from Model 1

WRG_family<- function(n, N) {
  w <- runif(N * N, -0.5, 0.5)
  center <- matrix(w, N, N)
  center <- (center + t(center)) / 2
  diag(center) <- 0

  my_array <- array(0, dim = c(n, N, N))
  my_array[1, , ] <- center

  for (i in 1:N) {
    for (j in 1:i) {
      vn1 <- rnorm(n, center[i, j], 0.1)
      my_array[-1, i, j] <- vn1[-1]
      #my_array[-1, j, i] <- vn1[-1]
    }
  }

  for (k in 2:n) {
    my_array[k, , ] <- my_array[k, , ] + t(my_array[k, , ])
    diag(my_array[k, , ]) <- 0
  }
  return(my_array)
}




###generate the matrix from Model 2


weight_out1_family<-function(n,N){
  w<-runif(N*N,0.5,1.5)
  center<-matrix(w,N,N)
  center <- (center + t(center)) / 2
  diag(center) <- 0
  my_array <- array(0, dim = c(n, N, N))
  my_array[1,,]<-center
  for (i in 1:N) {
    for (j in 1:i) {
      vn1 <- rnorm(n, center[i, j], 0.1)
      my_array[-1, i, j] <- vn1[-1]
    }
  }
  for (k in 2:n) {
    my_array[k, , ] <- my_array[k, , ] + t(my_array[k, , ])
    diag(my_array[k, , ]) <- 0
  }
  return(my_array)
}



###generate the matrix from Model 3
weight_out_family<-function(n,N){
  w<-runif(N*N,15,20)
  center<-matrix(w,N,N)
  center <- (center + t(center)) / 2
  diag(center) <- 0
  my_array <- array(0, dim = c(n, N, N))
  my_array[1,,]<-center
  for (i in 1:N) {
    for (j in 1:i) {
      vn1 <- rnorm(n, center[i, j], 0.5)
      my_array[-1, i, j] <- vn1[-1]
    }
  }
  for (k in 2:n) {
    my_array[k, , ] <- my_array[k, , ] + t(my_array[k, , ])
    diag(my_array[k, , ]) <- 0
  }
  return(my_array)
}



###generate the covariance matrix from Model 6
full_c<-function(N){
  c <- matrix(1,N,N)
  diag(c)<-2
  return(c)
}

###generate the covariance matrix from Model 7
Circle_C<- function(N){
  c <- diag(1,N)
  for (i in 1:(N-1)) {
    c[i+1,i]<- 0.5
    c[i,i+1]<- 0.5
  }
  c[1,N] <- c[N,1] <- 0.4
  return(c)
}


###generate the covariance matrix from Model 7
AR_C<- function(N){
  c <- diag(1,N)
  for (i in 1:(N-1)) {
    c[i+1,i]<- 0.5
    c[i,i+1]<- 0.5
  }
  for (i in 1:(N-2)) {
    c[i+2,i]<- 0.25
    c[i,i+2]<- 0.25
  }
  return(c)
}

###generate networks from Model 6
weight_family_full<-function(n,N){
  full_c_inv <- solve(full_c(N))
  my_array <- array(0, dim = c(n, N, N))
  for (i in 1:n) {
    data<-mvrnorm(3000, rep(0, N), full_c_inv)
    my_array[i,,]<-solve(var(data))
  }
  for (j in 1:n) {
    #my_array[j,,]<-my_array[j,,]+1
    diag(my_array[j,,])<-0
    my_array[j,,]<-round(my_array[j,,],3)
  }
  return(my_array)
}


###generate networks from Model 7
weight_family_out6<-function(n,N){
  full_c_inv <- solve(full_c(N))
  Circel_c_inv <- solve(Circle_C(N))
  AR_c_inv <- solve(AR_C(N))
  my_array <-array(0,dim = c(n, N, N))
  for (i in 1:n) {
    data<-mvrnorm(3000, rep(0, N), full_c_inv)
    #data[4:20,]<-data[4:20,]*10
    data[1:500,]<-mvrnorm(500, rep(0, N), Circel_c_inv)
    data[1500:2000,]<-mvrnorm(501, rep(0, N), AR_c_inv)
    my_array[i,,]<-solve(var(data))
    #my_array[i,,]<-my_array[i,,]
    diag(my_array[i,,])<-0
    my_array[i,,]<-round(my_array[i,,],3)
  }
  return(my_array)
}
