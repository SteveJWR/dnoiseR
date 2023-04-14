


# Function: applies the nonparametric EM algorithm to compute the possibly regularized npmle
# Input: p.hat: empirical distribution of the observed scores
#        cond: function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
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








# Function: Estimate the latent distribution using a convex solver
# Input: p.hat, a discrete empirical distribution function;  numeric
#        A.matrix, matrix which maps a vector of uniform bins to the observed data; matrix
#        mu,  regularization parameter; >= 0
# Output: latent: weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value

estimate_mixing_numeric <- function(p.hat, A.matrix, mu){

  R.bins <- ncol(A.matrix) # number of bins in the latent distribution
  N <- nrow(A.matrix) - 1 # number of scores
  theta <- Variable(R.bins, name = "latent discretized distribution") # values of the weight vector
  obs.dist <- Variable(N + 1, name = "model observed distribution") # values of the observed distribution

  data.obj <-  t(p.hat) %*% log(obs.dist)
  pen.obj <- (mu/R.bins)*t(rep(1, R.bins)) %*% log(theta)

  constraints <- list(
    obs.dist == A.matrix %*% theta,
    sum(theta) <= 1,
    mu/(R.bins*(1 + mu)) <= theta
  )

  obj.arg <- data.obj + pen.obj
  obj <- Maximize(obj.arg)
  p <- Problem(obj, constraints)


  #value(theta) <- rep(1/R.bins, R.bins) # initial guess of a uniform distribution
  result <- solve(p, solver = "MOSEK")
  #result <- solve(p, verbose = TRUE)

  p.m <- result$getValue(theta)
  p.ma <- result$getValue(obs.dist)
  out.list <- list("latent" = p.m, "observed" = p.ma)
  return(out.list)

}









