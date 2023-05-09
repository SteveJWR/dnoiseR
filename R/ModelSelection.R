
#' Generate List of Measurement Kernel Models
#'
#' @param N Number of questions
#' @param ker Kernel function
#' @param h.set Vector of bandwidths
#'
#' @return List of kernel functions
#' @export
#'
generate_mkm_list <- function(N,ker,h.set){
  out.list <- list()
  for(i in 1:length(h.set)){
    out.list[[i]] <- conditional_mkm(N,ker,h.set[i])
  }
  return(out.list)
}


#' @export
#'
error_model_selection_bivariate <- function(cond.set, Y, R.bins, cond.names, verbose = F){
  if(missing(cond.names)){
    cond.names = NULL
  }
  res <- rep(NA, length(cond.set))
  for(i in seq(length(cond.set))){
    if(verbose){
      cat(paste("Model",i,"/",length(cond.set)), end = "\r")
    }
    cond <- cond.set[[i]]
    A.matrix <- compute_A_matrix_2(R.bins, cond)
    A.tensor <- compute_A_tensor_2(R.bins, cond)
    lat.list <- estimate_mixing_npem_2(Y,A.matrix,A.tensor, mu = 0)
    mixture <- lat.list$latent
    res[i] <- population_bivariate_likelihood(Y,A.matrix,A.tensor,mixture, mu = 0)
  }
  i.max <- which.max(res)
  out.list <- list("opt_model" = cond.set[[i.max]], "lik_vals" = res, "opt_model_name" = cond.names[i.max])
  return(out.list)
}

#' @export
#'
mu_selection <- function(mu.set, cond, Y, R.bins, folds = 5, verbose = T){
  res <- rep(NA, length(mu.set))
  A.matrix <- compute_A_matrix_2(R.bins, cond)
  N = length(cond(0.5)) - 1
  n = nrow(Y)
  idx.folds <- caret::createFolds(seq(n),k = folds)
  idx = seq(n)

  for(i in seq(length(mu.set))){
    if(verbose){
      cat(paste("Model",i,"/",length(mu.set)), end = "\r")
    }

    mu <- mu.set[[i]]
    obj <- 0
    for(k in seq(folds)){
      val.idx = idx.folds[[k]]
      Y.val <- Y[val.idx,]
      train.idx <- setdiff(idx,val.idx)
      Y.train <- Y[train.idx,]
      p.hat.train <- compute_edf(Y.train[,1],N)
      p.hat.val <- compute_edf(Y.val[,1],N)
      model <- estimate_mixing_npem(p.hat.train,A.matrix, mu = mu)
      mixture <- as.vector(model$latent)
      p.ma <- as.numeric(A.matrix %*% mixture)
      obj <- obj - (1/folds)*kl_divergence(p.hat.val,p.ma)
    }
    res[i] <- obj
  }
  i.min <- which.min(res)
  out.list <- list("opt.mu" = mu.set[[i.min]], "cv.lik" = res)
  if(verbose){
    plot(mu.set,res)
  }
  return(out.list)
}

#' @export
#'
mu_selection_2 <- function(mu.set, cond, Y, R.bins, folds = 5, verbose = T){
  res <- rep(NA, length(mu.set))
  A.matrix <- compute_A_matrix_2(R.bins, cond)
  A.tensor <- compute_A_tensor_2(R.bins,cond)
  N = length(cond(0.5)) - 1
  n = nrow(Y)
  idx.folds <- caret::createFolds(seq(n),k = folds)
  idx = seq(n)

  for(i in seq(length(mu.set))){
    if(verbose){
      cat(paste("Model",i,"/",length(mu.set)), end = "\r")
    }

    mu <- mu.set[[i]]
    obj <- 0
    for(k in seq(folds)){
      val.idx = idx.folds[[k]]
      Y.val <- Y[val.idx,]
      train.idx <- setdiff(idx,val.idx)
      Y.train <- Y[train.idx,]

      model <- estimate_mixing_npem_2(Y.train,A.matrix,A.tensor, mu = mu)
      mixture <- as.vector(model$latent)
      obj <- obj + population_bivariate_likelihood(Y.val, A.matrix,A.tensor,mixture, mu = 0)
    }
    res[i] <- obj
  }
  i.min <- which.min(res)
  out.list <- list("opt.mu" = mu.set[[i.min]], "cv.lik" = res)
  if(verbose){
    plot(mu.set,res)
  }
  return(out.list)
}


##   Feasibility tests ---------

# Function: computes the first order feasibility test

#' First order feasibility test
#'
#' @param p.ma Nonparametric Maximum Likelihood Marginal Distribution
#' @param p.hat Empirical Distribution
#' @param sample.size Sample size
#'
#' @return p value for first order feasibility test
#' @export
#'
first_order_feasibility_test <- function (p.ma, p.hat, sample.size)
{
  lr <- 2 * sample.size * kl_divergence(p.hat, p.ma)
  k <- length(p.hat)
  lr.thresh <- multChernoff::criticalValue(k, sample.size,
                                           p = 10^(-8))
  if (lr >= lr.thresh) {
    p.feasibility <- 10^(-8)
  }
  tryCatch(p.feasibility <- multChernoff::tailProbBound(x = lr,
                                                        k = k, n = sample.size))
  if(is.nan(p.feasibility)){
    lr <- 2 * sample.size * kl_divergence(p.hat, p.ma)
    k <- nrow(p.hat) * ncol(p.hat)
    lr.thresh <- multChernoff::criticalValue(k, sample.size,
                                             p = 10**(-8), verbose = T)
    if (lr >= lr.thresh) {
      p.feasibility <- 10^(-8)
    }
    else {
      p.feasibility <- multChernoff::tailProbBound(x = lr,
                                                   k = k, n = sample.size)
      if (p.feasibility > 1) {
        p.feasibility <- 1
      }
      if (p.feasibility < 0) {
        p.feasibility <- 0
      }
    }
  }
  return(p.feasibility)
}






#' @export
#'
compute_edf_2 <- function (X, N, weights)
{
  n = nrow(X)
  if (missing(weights)) {
    weights <- rep(1, n)
  }
  p.hat <- matrix(data = 0, nrow = N + 1, ncol = N + 1)
  for (i in seq(n)) {
    p.hat[X[i, 1] + 1, X[i, 2] + 1] <- p.hat[X[i, 1] +
                                               1, X[i, 2] + 1] + weights[i]
  }
  p.hat = p.hat/sum(p.hat)
  return(p.hat)
}

function(x,N, weights){

  return(p.hat)
}

#' Second Order Feasibility Test
#'
#' @param latent.mixture Latent Distribution
#' @param A.tensor Conditional distribution tensor
#' @param p.hat Empirical Distribution
#' @param sample.size Sample Size
#'
#' @return p value for second order feasibility test
#' @export
#'
second_order_feasibility_test <- function(latent.mixture, A.tensor, p.hat, sample.size){
  # model implied distribution on Y
  R.bins <- length(latent.mixture)
  p.ma.list <- lapply(1:R.bins, function(z){
    out <- A.tensor[,,z]*latent.mixture[z]
    return(out)
  })
  p.ma <- matrix(data = 0, nrow = nrow(A.tensor[,,1]),
                 ncol = ncol(A.tensor[,,1]))

  for(i in 1:R.bins){
    p.ma <- p.ma + p.ma.list[[i]]
  }

  # likelihood ratio statistic
  lr <- 2*sample.size*kl_divergence(p.hat, p.ma)

  k <- nrow(p.hat)*ncol(p.hat)
  # threshold for a very low p.value.
  # numerical errors can occur if the likelihood ratio statistic is too large
  # for values which would correspond to p-values below 10^(-8) we simply use 10^(-8)
  tryCatch(p.feasibility <- multChernoff::tailProbBound(x = lr,
                                                        k = k, n = sample.size))
  if(is.nan(p.feasibility)){
    lr <- 2 * sample.size * kl_divergence(p.hat, p.ma)
    k <- nrow(p.hat) * ncol(p.hat)
    lr.thresh <- multChernoff::criticalValue(k, sample.size,
                                             p = 10**(-8), verbose = T)
    if (lr >= lr.thresh) {
      p.feasibility <- 10^(-8)
    }
    else {
      p.feasibility <- multChernoff::tailProbBound(x = lr,
                                                   k = k, n = sample.size)
      if (p.feasibility > 1) {
        p.feasibility <- 1
      }
      if (p.feasibility < 0) {
        p.feasibility <- 0
      }
    }
  }
  return(p.feasibility)
}




##   Model Selection Functions ----------------------------------------

# Function: Numerically approximate the distribution of second test score, given the first
# Input: y.obs, the observed score;  must be an integer
#        latent.mixture, weights assigning to the discretized latent distribution; vector summing to 1
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> P{0, ..., N}
# Output: probability : weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value



second_score_conditional <- function(y.obs, latent.mixture, cond){

  tau <- inv_quantiles_from_weights(latent.mixture)

  # corresponding quantiles to the latent mixture
  latent.quantiles <- seq(0,1, length.out = length(tau))

  R.bins <- length(latent.mixture)
  latent.points <- sapply(1:R.bins, function(z){
    gam <- (latent.quantiles[z] + latent.quantiles[z+1])/2
    return(gam)
  })


  joint.dist.grid  <- sapply(1:R.bins, function(gam.idx){
    gam <- latent.points[gam.idx]
    weight <- latent.mixture[gam.idx]

    out.y1y2.joint <- cond(gam) %o% cond(gam)
    prob.row <- out.y1y2.joint[,y.obs + 1]
    out <- prob.row*weight
    return(out)
  })

  cond.prob <- rowSums(joint.dist.grid)
  cond.prob <- cond.prob/sum(cond.prob)
  return(cond.prob)
}


# Function: Compute the likelihood of a mixing distributions based on two sequential observations for a constant gamma value
# Input: y.paired, a data.frame of the pairs of observed scores
#        A.tensor, a tensor which is used to map the latent gamma to the joint bivariate distribution of scores
#        latent.mixture.list,  a list of latent mixtures to use for each pair of observations.
# Output: log-likelihood value


compute_two_obs_loglikelihood <- function(y.paired, A.tensor, latent.mixture.list){
  M <- nrow(y.paired)
  log.likelihood <- 0
  for(m in 1:M){
    pair.prob <- t(A.tensor[y.paired[m,1] +1,y.paired[m,2] + 1,]) %*% latent.mixture.list[[m]]
    log.likelihood <- log.likelihood + as.numeric(log(pair.prob))
  }
  return(log.likelihood)
}


# Function: Compute the joint distribution of the observed scores for a harmonizable set of beta distribuions
# Input: beta*.model.w, beta * (1 or 2) parameters for the w (y or z branch ) of the distribution
#        cond.w, w (y or z branch ) function defining the conditional distribution of Y|gamma; function [0,1] -> P{0, ..., N}
#        grid.size,  a grid size for the numerical approximation of the distribution
# Output: joint distribution matrix
compute_joint_dist_beta <- function(beta1.model.y, beta2.model.y,
                                    beta1.model.z, beta2.model.z,
                                    cond.y, cond.z, grid.size = 10000){

  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1

  latent.grid <- seq(0,1, length.out = grid.size)

  joint.dist <- lapply(latent.grid, function(x){
    latent.y <- qbeta(x, shape1 = beta1.model.y, shape2 = beta2.model.y)
    latent.z <- qbeta(x, shape1 = beta1.model.z, shape2 = beta2.model.z)
    iy = 0:Ny
    iz = 0:Nz
    assumption.weights.y <- cond.y(latent.y)
    assumption.weights.z <- cond.z(latent.z)

    prob.y <- assumption.weights.y/sum(assumption.weights.y)
    prob.z <- assumption.weights.z/sum(assumption.weights.z)

    p.yz <- prob.y %*% t(prob.z)
    out <- p.yz

  })

  p.yz <- matrix(data = 0, nrow = Ny + 1, ncol = Nz + 1)
  for(i in 1:length(joint.dist)){
    p.yz <- p.yz + joint.dist[[i]]
  }
  p.yz <- p.yz/grid.size
  return(p.yz)
}







