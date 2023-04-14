
##   Feasibility tests ---------

# Function: computes the first order feasibility test
# Input: latent.mixture, a vector of weights for the latent distribution; numeric
#        A.matrix, matrix which maps a vector of uniform bins to the observed data; matrix
#        p.hat, a discrete empirical distribution function;  must be an vector of numeric values
# Output: p.feasibility: p value for the feasibility test, is at minimum e-8 due to numerical rounding issues

test_feasibility_first_order <- function(latent.mixture, A.matrix, p.hat, sample.size){
  # model implied distribution on Y
  p.ma <- A.matrix %*% latent.mixture
  p.ma <- as.numeric(p.ma)
  # likelihood ratio statistic
  lr <- 2*sample.size*kl_divergence(p.hat, p.ma)
  k <- length(p.hat)

  # threshold for a very low p.value.
  # numerical errors can occur if the likelihood ratio statistic is too large
  # for values which would correspond to p-values below 10^(-8) we simply use 10^(-8)
  lr.thresh <- criticalValue(k,sample.size, p = 10^(-8))
  if(lr >= lr.thresh){
    p.feasibility <- 10^(-8)
  } else {
    p.feasibility <- tailProbBound(x = lr, k = k, n = sample.size)
    if(p.feasibility > 1){
      p.feasibility <- 1
    }
    # handling rounding errors
    if(p.feasibility < 0){
      p.feasibility <- 0
    }
  }
  return(p.feasibility)
}


# Function: computes the second order feasibility test
# Input: latent.mixture, a vector of weights for the latent distribution; numeric
#        A.two.sample.tensor, tensor which maps a vector of uniform bins to the bivariate observed data distribution; 3D tensor
#        p.hat, a discrete empirical distribution function;  must be an vector of numeric values
# Output: p.feasibility: p value for the feasibility test, is at minimum e-8 due to numerical rounding issues

test_feasibility_second_order <- function(latent.mixture, A.two.sample.tensor, p.hat, sample.size){
  # model implied distribution on Y
  R_bins <- length(latent.mixture)
  p.ma.list <- lapply(1:R_bins, function(z){
    out <- A.two.sample.tensor[,,z]*latent.mixture[z]
    return(out)
  })
  p.ma <- matrix(data = 0, nrow = nrow(A.two.sample.tensor[,,1]),
                 ncol = ncol(A.two.sample.tensor[,,1]))

  for(i in 1:R_bins){
    p.ma <- p.ma + p.ma.list[[i]]
  }

  # likelihood ratio statistic
  lr <- 2*sample.size*kl_divergence(p.hat, p.ma)

  k <- nrow(p.hat)*ncol(p.hat)
  # threshold for a very low p.value.
  # numerical errors can occur if the likelihood ratio statistic is too large
  # for values which would correspond to p-values below 10^(-8) we simply use 10^(-8)
  lr.thresh <- criticalValue(k,sample.size, p = 10^(-8))
  if(lr >= lr.thresh){
    p.feasibility <- 10^(-8)
  } else {
    p.feasibility <- tailProbBound(x = lr, k = k, n = sample.size)
    if(p.feasibility > 1){
      p.feasibility <- 1
    }
    # handling rounding errors
    if(p.feasibility < 0){
      p.feasibility <- 0
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

  R_bins <- length(latent.mixture)
  latent.points <- sapply(1:R_bins, function(z){
    gam <- (latent.quantiles[z] + latent.quantiles[z+1])/2
    return(gam)
  })


  joint.dist.grid  <- sapply(1:R_bins, function(gam.idx){
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
#        A.two.sample.tensor, a tensor which is used to map the latent gamma to the joint bivariate distribution of scores
#        latent.mixture.list,  a list of latent mixtures to use for each pair of observations.
# Output: log-likelihood value


compute_two_obs_loglikelihood <- function(y.paired, A.two.sample.tensor, latent.mixture.list){
  M <- nrow(y.paired)
  log.likelihood <- 0
  for(m in 1:M){
    pair.prob <- t(A.two.sample.tensor[y.paired[m,1] +1,y.paired[m,2] + 1,]) %*% latent.mixture.list[[m]]
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







