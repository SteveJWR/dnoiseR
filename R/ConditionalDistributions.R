# Function: conditional distribution generated according to the measurement kernel model
# Input: N;  must be an integer
#        ker, number of repeated observations per individual; integer
#        h, bandwidth; > 0
# Output: Conditional Distribution Function for the Measurement Kernel Model

#' Conditional distribution function for the measurement kernel model
#'
#' @param N Number of Questions
#' @param ker Measurement Kernel
#' @param h Bandwidth
#'
#' @return Conditional distribution function for the measurement kernel model
#' @export
#'
conditional_mkm <- function(N, ker, h){
  #ensuring these don't change when changing definitions of N,ker and h
  force(N)
  force(ker)
  force(h)
  i <- 0:N
  out.function <- function(gam){
    out <- ker((i - N*gam)/h)/sum(ker((i - N*gam)/h))
    return(out)
  }
  return(out.function)
}


#' Conditional distribution function for the binomial model
#'
#' @param N Number of Questions
#'
#' @return Conditional distribution function for the binomial model
#' @export
#'
generate_cond_binomial <- function(N){
  out.function <- function(x){
    dbinom(0:N,size = N, prob = x)
  }
}


##   Latent model fitting functions  ----------------------------------------


# Function: sample from the latent distribution given an observed score.
# Input: y.obs, the observed score value
#        n.mc.samp,  number of samples from gamma|Y
#        tau: a set of inverse quantiles corresponding to the edges of the bins of the latent distribution
#        latent.mixture, input is a list of weights assigned to a uniformly binned latent distribution
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
# Output: a vector of samples from gamma|Y

# sample_latent_conditional <- function(y.obs, n.mc.samp, tau, latent.trait.quantiles, cond){
#
#   N <- length(cond(0.5)) - 1
#
#   # approximate the observed frequencies
#   p.ma <- compute_p_ma(tau,
#                        latent.trait.quantiles,
#                        cond)
#
#   p.approx <- p.ma[y.obs + 1] # probability of conditional distribution
#
#   # approximate number of times required to sample
#   n.parallel <- round(2*n.mc.samp/p.approx)
#
#   # number of bins
#   n.bins <- length(tau) - 1
#   weights <- c()
#   for(i in 1:n.bins){
#     weights[i] <- tau[i + 1] - tau[i]
#   }
#   if(n.parallel >= 500000){
#     n.parallel <- 500000
#   }
#   # outcome
#   out <- c()
#   while (length(out) < n.mc.samp){
#
#     # sampling from the binned latent distribution
#     latent.idx <- sample(1:n.bins, size = n.parallel, replace = TRUE, prob = weights)
#     low.bounds <- latent.trait.quantiles[latent.idx]
#     high.bounds <- latent.trait.quantiles[latent.idx + 1]
#     latent.sample <- runif(n.parallel, min = low.bounds, max = high.bounds)
#
#     y.model.samp <- sapply(latent.sample, function(x){
#       i = 0:N
#       assumption.weights <- cond(x)
#       out <- sample(0:N, size = 1, replace = TRUE, prob = assumption.weights)
#       return(out)
#     })
#
#     # keeping the latent variables which generated the observed score
#     keep.latent.idx <- (y.model.samp == y.obs)
#
#     out <- c(out, latent.sample[keep.latent.idx])
#   }
#   # returning only the first exact number of samples
#   out <- out[1:n.mc.samp]
#   return(out)
#
# }


# Function: sample from the next EM step for updating the latent distribution.
# Input: y.obs: the observed score value
#        mu:  regularization parameter
#        tau: a set of inverse quantiles corresponding to the edges of the bins of the latent distribution
#        latent.mixture: input is a list of weights assigned to a uniformly binned latent distribution
#        cond: function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        n.mc.samp:  number of monte carlo samples from the em update step
# Output: a vector of samples from gamma|Y

# sample_mc_em <- function(p.hat, mu, tau, latent.trait.quantiles, cond, n.mc.samp){
#   # tau and quantiles must be of the same length
#   # must be ordered
#
#   N <- length(cond(0.5)) - 1
#
#   n.tau <- length(tau) - 1
#   n.parallel <- n.mc.samp # blocking sampling
#
#   weights <- c()
#   for(i in 1:n.tau){
#     weights[i] <- tau[i + 1] - tau[i]
#   }
#   # weights based on differences of quantiles
#
#   gamma.em.step <- c()
#
#   while (length(gamma.em.step) < n.mc.samp){
#     y.hat.samp <- sample(0:N, size = n.parallel, replace = TRUE, prob = p.hat)
#
#     frequencies <- n.parallel*compute_edf(y.hat.samp,N)
#     # monte carlo sampled from the nonparametric em algorithm
#     gamma.em.step.samp <- lapply(0:N, function(x){
#       num.samples <- frequencies[x + 1]
#
#       if(num.samples > 0){
#         out <- sample_latent_conditional(y.obs = x,
#                                          n.mc.samp = num.samples,
#                                          tau = tau,
#                                          latent.trait.quantiles = latent.trait.quantiles,
#                                          cond = cond)
#       } else {
#         out <- c()
#       }
#
#       return(out)
#     })
#     gamma.em.step.samp <- unlist(gamma.em.step.samp)
#     gamma.em.step <- c(gamma.em.step, gamma.em.step.samp)
#
#   }
#   gamma.em.step <- gamma.em.step[1:n.mc.samp]
#
#   if(mu > 0){
#     n.mixture <- round(n.mc.samp*mu) # round to nearest integer
#     uniform.mixture <- runif(n.mixture)
#     out <- c(gamma.em.step,uniform.mixture)
#   } else {
#     out <- gamma.em.step
#   }
#   return(out)
# }




# Oldname: compute_loglikelihood_from_latent
# Function: computes the loglikelihood value from
# Input: y.obs: the observed score value
#        mu:  regularization parameter
#        tau: a set of inverse quantiles corresponding to the edges of the bins of the latent distribution
#        latent.mixture: input is a list of weights assigned to a uniformly binned latent distribution
#        cond: function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        n.mc.samp:  number of monte carlo samples from the em update step
# Output: a vector of samples from gamma|Y



# compute_loglikelihood_from_latent <- function(p.hat, p.ma, tau,
#                                               latent.trait.quantiles,
#                                               mu = 0){
#
#   R_bins <- length(tau) - 1
#
#   heights <- c()
#   for(i in 1:R_bins){
#     heights[i] <- (tau[i + 1] - tau[i])/(latent.trait.quantiles[i + 1] - latent.trait.quantiles[i])
#   }
#
#   widths <- c()
#   for(i in 1:R_bins){
#     widths[i] <- latent.trait.quantiles[i + 1] - latent.trait.quantiles[i]
#   }
#
#   if(mu > 0){
#     reg.term <- mu*sum(log(heights)*widths)
#   } else {
#     reg.term <- 0
#   }
#
#
#   data.term.vec <- p.hat*log(p.ma)
#   # defining nan (0*log(0)) values as 0
#   data.term.vec[is.nan(data.term.vec)] <- 0
#
#   out <- sum(data.term.vec) + reg.term
#
#   return(out)
#
# }




# Function: compute the A matrix for a binned latent approximation
# Input: R_bins, the number of uniform latent bins  must an integer
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        numeric.points,  number of points used in the numeric approximation of A
# Output: latent: weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value



#' Compute conditional distribution matrix
#'
#' @param R_bins The number of latent bins to choose
#' @param cond Conditional distribution
#' @param verbose Verbose parameter
#'
#' @return Conditional Distribution Matric Form
#' @export
#'
compute_A_matrix <- function(R_bins, cond, numeric.points = 100, verbose = F){
  N <- length(cond(0.5)) - 1
  A_matrix <- matrix(data = NA, nrow = N + 1, ncol = R_bins)
  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
    A.row <- rep(0, N+1)
    for(j in 1:numeric.points){
      y = 0:N
      weights <- cond(design.points[j])
      A.row <- A.row + weights/(sum(weights))
    }
    A.row <- A.row/sum(A.row)
    A_matrix[,i] = A.row
    if(verbose){
      cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
    }
  }
  return(A_matrix)
}




# Function: compute the A tensor using a binned latent approximation
# Input: R_bins, the number of uniform latent bins  must an integer
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        numeric.points,  number of points used in the numeric approximation of A
# Output: latent: weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value

#' @export
#'
compute_A_tensor <- function(R_bins, cond, numeric.points = 100, verbose = F){
  # 3d Array for faster computation.
  N <- length(cond(0.5)) - 1
  A_3D <- array(NA, dim = c(N+1,N+1,R_bins))

  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)

    A.block <- array(0, dim = c(N+1,N+1))
    y1 = 0:N
    y2 = 0:N

    for(j in 1:numeric.points){

      weights <- outer(cond(design.points[j]), cond(design.points[j]), "*")
      A.block <- A.block + weights/(sum(weights))
    }
    A_3D[,,i] <- A.block/(sum(A.block))
    if(verbose){
      cat(paste0("A Tensor Computed Row: ", i,"/",R_bins), end="\r")
    }
  }

  return(A_3D)
}



# faster than the previous version
#' @export
#'
compute_A_matrix_2 <- function(R.bins, cond){
  grid <- seq(0,R.bins -1)/(R.bins - 1)
  A.out <- sapply(grid,function(z){
    out <- cond(z)
    out <- out/sum(out)
    return(out)
  })
  return(A.out)
}

#' @export
#'
compute_A_tensor_2 <- function(R.bins, cond){
  grid <- seq(0,R.bins -1)/(R.bins - 1)
  #grid <- array(grid,c(length(grid),1,1))
  N <- length(cond(0.5)) - 1
  A.out <- vapply(grid,function(z){
    out <- cond(z)
    out.mat <- outer(out,out,"*")
    out.mat <- out.mat/sum(out.mat)
    return(out.mat)
  }, FUN.VALUE = array(0,c(N + 1,N + 1)))
  return(A.out)
}


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
