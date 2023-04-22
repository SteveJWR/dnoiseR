
#' Simulate test outcomes
#' @export
#'
sample_bivariate <- function(n,p.mat){
  # column stacking
  p.vec <- as.vector(p.mat)
  N <- dim(p.mat)[1] - 1
  k <- sample(seq((N + 1)^2), n,p.vec, replace = T)
  y1 <- (k - 1) %% (N + 1)
  y2 <- floor((k - 1) /(N + 1))
  Y <- matrix(c(y1,y2), ncol = 2, byrow = F)
  return(Y)
}

#' Simulate test outcomes
#' @export
#'
simulate_test <- function(n2,n1, cond, mixture){
  R.bins <- length(mixture)
  A.tens <- compute_A_tensor_2(R.bins, cond)
  A.mat <- compute_A_matrix_2(R.bins, cond)

  p.ma <- A.mat %*% mixture
  p.ma2 <- tensor_prod(A.tens, mixture)

  N = length(cond(0.5)) - 1
  y1 <- sample(seq(0,N),n1,prob = p.ma, replace = T)
  Y2 <- sample_bivariate(n2,p.ma2)
  Y1 <- matrix(c(y1,rep(NA,length(y1))), ncol = 2, byrow = F)
  Y <- rbind(Y2,Y1)

  return(Y)
}

#' Simulate test outcomes
#' @export
#'
simulate_test_cond <- function(obs.set,cond, gamma.set){
  Y <- matrix(NA, nrow = 0, ncol = 2)
  for(i in seq(length(gamma.set))){
    gamma = gamma.set[i]
    obs = obs.set[i]

    p <- cond(gamma)
    p2 <- outer(p,p,"*")
    N <- length(cond(0.5)) - 1
    if(obs == 1){
      y <- sample(seq(0,N),1,prob = p, replace = T)
      y <- c(y,NA)
    } else if(obs == 2){
      y <- sample_bivariate(1,p2)
    }
    Y <- rbind(Y,y)
  }
  return(Y)
}



