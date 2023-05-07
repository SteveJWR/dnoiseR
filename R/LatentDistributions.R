#TODO: Change the name in the simulations to the fit npem algorithm

# Function: applies the nonparametric EM algorithm to compute the possibly regularized npmle
# Input: p.hat: empirical distribution of the observed scores
#        cond: function defining the conditional distribution of Y|gamma
#        R_bins:  number of bins (quantiles) to sample from
#        mu:  regularization parameter
#        n.mc.samp: number of monte-carlo samples per iteration of the nonparametric em algorithm
#        verbose: indicator for whether rolling number of iterations will be displayed
# Output: sequence of quantiles for a uniform partition of the latent distribution


fit_nonpar_em <- function(p.hat, cond,  R_bins = 300, mu = 0, n.mc.samp = 1000, verbose = T){
  N <- length(cond(0.5)) - 1
  # looking at a uniform spacing of the inverse quantiles.
  tau <- seq(0,1,length.out = R_bins)
  # initialize with a uniform distribution
  latent.trait.quantiles.init <- tau

  # likelihood change threshold
  threshold <- 10^(-6)


  # initialize with a uniform distribution
  current.quantiles <- latent.trait.quantiles.init

  n.quantiles <- length(current.quantiles)

  # computing the implied distribution on the observed data
  p.ma <- compute_p_ma(tau = tau,
                       latent.trait.quantiles = current.quantiles,
                       cond = cond,
                       numeric.points = 100)


  prev.likelihood <- compute_loglikelihood_from_latent(p.hat = p.hat,
                                                       p.ma = p.ma,
                                                       tau = tau,
                                                       latent.trait.quantiles = current.quantiles,
                                                       mu = mu)

  diff.likelihood <- Inf
  em.steps <- 1

  while(diff.likelihood >= threshold){
    latent.sample <- sample_mc_em(p.hat = p.hat,
                                  mu = mu,
                                  tau = tau,
                                  latent.trait.quantiles = current.quantiles,
                                  cond = cond,
                                  n.mc.samp = n.mc.samp)

    # computing updated set of quantiles
    current.quantiles <- quantile(latent.sample, tau, names = FALSE)
    current.quantiles[1] <- 0
    current.quantiles[n.quantiles] <- 1

    p.ma <- compute_p_ma(tau = tau,
                         latent.trait.quantiles = current.quantiles,
                         cond = cond,
                         numeric.points = 100)

    new.likelihood <- compute_loglikelihood_from_latent(p.hat = p.hat,
                                                        p.ma = p.ma,
                                                        tau = tau,
                                                        latent.trait.quantiles = current.quantiles,
                                                        mu = mu)

    diff.likelihood <- new.likelihood - prev.likelihood
    prev.likelihood <- new.likelihood
    if(verbose){
      print(paste0("EM Steps: ", em.steps, " || Likelihood Change: ", diff.likelihood))
    }
    em.steps <- em.steps + 1
  }
  if(verbose){
    print(paste0("Stochastic EM converged"))
  }
  return(current.quantiles)
}



#' Select Number of Bins
#'
#' @param p.hat Distribution to approximate
#' @param cond Conditional distribution (Empirical Bayes f model)
#' @param mu Regularization Parameter
#' @param R.min Minimum number of bins
#' @param R.max Maximum number of bins
#' @param threshold Threshold for log-likelihood difference
#' @param num.steps Number of steps to take for the R (number of bins) selection
#' @param max.iter Maximum iterations of the NPEM algorithm.
#'
#' @return Smallest viable value of R
#'
select_R <- function(p.hat, cond, mu = 0.001, R.min = 100, R.max = 5000, threshold = 10**(-4), num.steps = 20, max.iter = 50){

  R.start <- R.min
  R.seq = round(seq(R.min,R.max, length.out = num.steps))
  diff <- Inf
  A.matrix <- compute_A_matrix_2(R.min, cond)
  lat.list <- estimate_mixing_npem(p.hat,A.matrix, mu = mu, threshold = threshold, max.iter = max.iter)
  mixture <- lat.list$latent
  p.ma = A.matrix %*% mixture
  uniform.latent = rep(1/(R.start), R.start)
  if (mu != 0) {
    like.old <- -kl_divergence(p.hat, p.ma) - mu *
      kl_divergence(uniform.latent, mixture)
  }
  else {
    like.old <- -kl_divergence(p.hat, p.ma)
  }
  for(r in seq(2,num.steps)){
    R.new = R.seq[r]
    A.matrix <- compute_A_matrix_2(R.new, cond)
    lat.list <- estimate_mixing_npem(p.hat,A.matrix, mu = mu, threshold = threshold, max.iter = max.iter)
    mixture <- lat.list$latent
    p.ma = A.matrix %*% mixture
    uniform.latent = rep(1/(R.new), R.new)

    if (mu != 0) {
      like.new <- -kl_divergence(p.hat, p.ma) - mu *
        kl_divergence(uniform.latent, mixture)
    }
    else {
      like.new <- -kl_divergence(p.hat, p.ma)
    }
    diff <- like.new - like.old
    like.old <- like.new
    cat(paste("Loglikelihood error differenc:", diff, "Num Bins: ", R.new), end = "\r")
    if(diff < threshold){
      break
    }
  }

  paste("Selected number of bins:",R.new)
  return(R.new)
}




#' CVX implementation
#'
#' @param p.hat Discrete empirical distribution function
#' @param A.matrix  Conditional distribution Matrix
#' @param mu Regularization Parameter
#' @param cvx.solver Which solver to use
#'
#' @return latent trait and observed distribution
#' @export
#'
estimate_mixing_numeric <- function(p.hat, A.matrix, mu = 0, cvx.solver = "SCS"){

  R.bins <- ncol(A.matrix) # number of bins in the latent distribution
  N <- nrow(A.matrix) - 1 # number of scores
  theta <- CVXR::Variable(R.bins, name = "latent discretized distribution") # values of the weight vector
  obs.dist <- CVXR::Variable(N + 1, name = "model observed distribution") # values of the observed distribution

  data.obj <-  t(p.hat) %*% log(obs.dist)
  pen.obj <- (mu/R.bins)*t(rep(1, R.bins)) %*% log(theta)

  constraints <- list(
    obs.dist == A.matrix %*% theta,
    sum(theta) <= 1,
    mu/(R.bins*(1 + mu)) <= theta
  )

  obj.arg <- data.obj + pen.obj
  obj <- CVXR::Maximize(obj.arg)
  p <- CVXR::Problem(obj, constraints)


  #value(theta) <- rep(1/R.bins, R.bins) # initial guess of a uniform distribution
  result <- CVXR::solve(p, solver = cvx.solver)
  #result <- solve(p, verbose = TRUE)

  p.m <- result$getValue(theta)
  p.ma <- result$getValue(obs.dist)
  out.list <- list("latent" = p.m, "observed" = p.ma)
  return(out.list)
}




#' Nonparametric EM Algorithm for latent trait
#'
#' @param p.hat Discrete empirical distribution function
#' @param A.matrix  Conditional distribution Matrix
#' @param mu Regularization Parameter
#' @param threshold Threshold for the change in the likelihood
#' @param verbose Whether to display progress
#' @param max.iter Maximum number of iterations of the EM algorithm
#' @param init.latent Initialize latent distribution.  This must be
#'
#' @return latent trait and observed distribution
#' @export
#'
estimate_mixing_npem <- function (p.hat, A.matrix, mu = 0, threshold = 5 * 10^(-5), verbose = F,
                                  max.iter = 50, init.latent)
{
  R_bins = ncol(A.matrix)
  uniform.latent <- rep(1, R_bins)
  uniform.latent <- uniform.latent/sum(uniform.latent)
  if(missing(init.latent)){
    latent.trait.init <- uniform.latent
  } else {
    latent.trait.init <- init.latent
  }

  latent.trait <- latent.trait.init
  p.ma.init <- A.matrix %*% latent.trait.init
  like.old <- -Inf
  diff = Inf
  iter = 1
  while (abs(diff) > threshold & iter <= max.iter) {
    p.y.gamma.joint <- t(t(A.matrix) * latent.trait)
    p.gamma.given.y <- p.y.gamma.joint/rowSums(p.y.gamma.joint)
    latent.trait.new <- as.numeric(t(p.gamma.given.y) %*%
                                     p.hat * (1/(1 + mu))) + (mu/(1 + mu)) * uniform.latent
    latent.trait.new = latent.trait.new/sum(latent.trait.new)
    p.ma.new = A.matrix %*% latent.trait.new
    if (mu != 0) {
      like.new <- -kl_divergence(p.hat, p.ma.new) - mu *
        kl_divergence(uniform.latent, latent.trait.new)
    }
    else {
      like.new <- -kl_divergence(p.hat, p.ma.new)
    }
    diff = like.new - like.old
    if (verbose) {
      if (iter%%10 == 1) {
        cat(paste0("Likelihood Change: ", diff), end = "\r")
      }
    }
    delta.gamma.vec = latent.trait.new - latent.trait.new
    latent.trait = latent.trait.new
    like.old = like.new
    iter = iter + 1
  }
  out.list <- list(latent = latent.trait, observed = as.numeric(p.ma.new))
  return(out.list)
}


# npem algorithm for the bivariate likelihood part 2
#' @export
#'
estimate_mixing_npem_2 <- function(Y, A.matrix, A.tensor, weights, mu = 0, threshold = 5*10**(-5),  verbose = F, max.iter = 50) {

  if(ncol(Y) != 2){
    stop("Y must have 2 columns" )
  }
  if(any(is.na(Y[,1]))){
    stop("Y cannot have missingness in first column")
  }
  if(missing(weights)){
    weights <- rep(1,nrow(Y))
  }
  if(missing(mu)){
    mu = 0
  }
  n.obs <- nrow(Y)
  N <- nrow(A.matrix) - 1 # number of scores
  R.bins <- ncol(A.matrix)
  count.vector <- rep(0,N + 1)
  count.matrix <- matrix(0,N+1, N+1)


  #
  uniform.latent <- rep(1,R_bins)
  uniform.latent <- uniform.latent/sum(uniform.latent)

  latent.trait.init <- uniform.latent



  A.ten.mat <- A.tensor
  # converting the dimension of the tensor to a matrix

  dim(A.ten.mat) <- c((N + 1)^2,R.bins)

  #NA's can be allowed in second column only
  for(i in seq(nrow(Y))){
    if(!any(is.na(Y[i,]))){
      count.matrix[Y[i,1] + 1,Y[i,2] + 1] <- count.matrix[Y[i,1] + 1,Y[i,2] + 1] + weights[i]
    } else {
      count.vector[Y[i,1] + 1] <- count.vector[Y[i,1] + 1] + weights[i]
    }

  }

  # need to vectorize the matrix counts
  count.mat.vec <- count.matrix
  dim(count.mat.vec) <- (N + 1)^2
  n.w = sum(weights)

  n.w1 = sum(count.vector)
  n.w2 = sum(count.matrix)

  p.hat.1 = as.numeric(count.vector/n.w1)
  p.hat.2 = as.numeric(count.mat.vec/n.w2)

  like.old <- -Inf
  latent.trait = latent.trait.init
  delta.gamma.vec = rep(0,R.bins)

  diff = Inf
  iter = 1
  while(abs(diff) > threshold & iter <= max.iter){
    p.y.gamma.joint1 <- t(t(A.matrix) * latent.trait)
    p.y.gamma.joint2 <- t(t(A.ten.mat) * latent.trait)

    p.gamma.given.y1 <- p.y.gamma.joint1 / rowSums(p.y.gamma.joint1)
    p.gamma.given.y2 <- p.y.gamma.joint2 / rowSums(p.y.gamma.joint2)

    latent.trait.new <- as.numeric((t(p.gamma.given.y1) %*% count.vector + t(p.gamma.given.y2) %*% count.mat.vec)) * (1/(n.w*(1 + mu))) + (mu/(1 + mu))*uniform.latent #+ momentum*delta.gamma.vec
#    latent.trait.new[latent.trait.new < (mu/((1 + mu)*R.bins))] = (mu/((1 + mu)*R.bins))
    latent.trait.new = latent.trait.new/sum(latent.trait.new)


    p.ma.new1 = A.matrix %*% latent.trait.new
    p.ma.new2 = A.ten.mat %*% latent.trait.new

    like.new <- -kl_divergence(p.hat.1,p.ma.new1)  -kl_divergence(p.hat.2,p.ma.new2) - mu* kl_divergence(uniform.latent, latent.trait.new)

    diff = like.new - like.old
    if(verbose){
      if(iter %% 10 == 1){
        cat(paste0("Likelihood Change: ", diff), end = "\r")
      }

      # if(iter %% 100 == 1){
      #   plot(latent.trait)
      #   plot(latent.trait.new)
      # }
    }
    delta.gamma.vec = latent.trait.new - latent.trait
    latent.trait <- latent.trait.new
    like.old <- like.new
    iter = iter + 1
  }

  out.list <- list("latent" = latent.trait, "observed" = as.numeric(p.ma.new1))
  return(out.list)
}


#' @export
#'
estimate_mixing_numeric_2 <- function(Y, A.matrix, A.tensor, mu, weights,cvx.solver = "SCS"){


  if(ncol(Y) != 2){
    stop("Y must have 2 columns" )
  }
  if(any(is.na(Y[,1]))){
    stop("Y cannot have missingness in first column")
  }
  if(missing(weights)){
    weights <- rep(1,nrow(Y))
  }
  if(missing(mu)){
    mu = 0
  }
  n.obs <- nrow(Y)
  N <- nrow(A.matrix) - 1 # number of scores
  R.bins <- ncol(A.matrix)
  count.vector <- rep(0,N + 1)
  count.matrix <- matrix(0,N+1, N+1)


  A.ten.mat <- A.tensor
  # converting the dimension of the tensor to a matrix

  dim(A.ten.mat) <- c((N + 1)^2,R.bins)

  #NA's can be allowed in second column only
  for(i in seq(nrow(Y))){
    if(!any(is.na(Y[i,]))){
      count.matrix[Y[i,1] + 1,Y[i,2] + 1] <- count.matrix[Y[i,1] + 1,Y[i,2] + 1] + weights[i]
    } else {
      count.vector[Y[i,1] + 1] <- count.vector[Y[i,1] + 1] + weights[i]
    }

  }
  # need to vectorize the matrix counts
  count.mat.vec <- count.matrix
  dim(count.mat.vec) <- (N + 1)^2

  theta <- CVXR::Variable(R.bins, name = "latent discretized distribution") # values of the weight vector

  data.obj1 <-  t(count.vector) %*% log(A.matrix %*% theta)
  data.obj2 <-  t(count.mat.vec) %*% log(A.ten.mat %*% theta)
  pen.obj <- n.obs*(mu/R.bins)*t(rep(1, R.bins)) %*% log(theta)

  constraints <- list(
    #obs.dist == A.matrix %*% theta,
    sum(theta) <= 1#,
    #mu/(R.bins*(1 + mu)) <= theta
  )

  obj.arg <- data.obj1 + data.obj2 + pen.obj
  obj <- CVXR::Maximize(obj.arg)
  prob <- CVXR::Problem(obj, constraints)


  #value(theta) <- rep(1/R.bins, R.bins) # initial guess of a uniform distribution
  result <- CVXR::solve(prob, solver = cvx.solver)
  #result <- solve(p, verbose = TRUE)

  p.m <- result$getValue(theta)
  p.ma <- A.matrix %*% p.m
  out.list <- list("latent" = p.m, "observed" = p.ma)
  return(out.list)
}


#' @export
population_bivariate_likelihood <- function(Y,A.matrix,A.tensor,mixture, weights, mu = 0){
  if (ncol(Y) != 2) {
    stop("Y must have 2 columns")
  }
  if (any(is.na(Y[, 1]))) {
    stop("Y cannot have missingness in first column")
  }
  if (missing(weights)) {
    weights <- rep(1, nrow(Y))
  }
  if (missing(mu)) {
    mu = 0
  }
  n.obs <- nrow(Y)
  N <- nrow(A.matrix) - 1
  R.bins <- ncol(A.matrix)
  count.vector <- rep(0, N + 1)
  count.matrix <- matrix(0, N + 1, N + 1)
  uniform.latent <- rep(1, R.bins)
  uniform.latent <- uniform.latent/sum(uniform.latent)
  latent.trait.init <- uniform.latent
  A.ten.mat <- A.tensor
  dim(A.ten.mat) <- c((N + 1)^2, R.bins)
  for (i in seq(nrow(Y))) {
    if (!any(is.na(Y[i, ]))) {
      count.matrix[Y[i, 1] + 1, Y[i, 2] + 1] <- count.matrix[Y[i,
                                                               1] + 1, Y[i, 2] + 1] + weights[i]
    }
    else {
      count.vector[Y[i, 1] + 1] <- count.vector[Y[i, 1] +
                                                  1] + weights[i]
    }
  }
  count.mat.vec <- count.matrix
  dim(count.mat.vec) <- (N + 1)^2
  n.w = sum(weights)
  n.w1 = sum(count.vector)
  n.w2 = sum(count.matrix)
  p.ma1 = as.numeric(A.matrix %*% mixture)
  p.ma2 = as.numeric(A.ten.mat %*% mixture)
  uniform = rep(1/R.bins, R.bins)
  if(mu != 0 ){
    like <- sum(count.vector * log(p.ma1)) + sum(count.mat.vec * log(p.ma2)) + mu*kl_divergence(uniform,mixture)
  } else {
    like <- sum(count.vector * log(p.ma1)) + sum(count.mat.vec * log(p.ma2))
  }
  return(like)
}



