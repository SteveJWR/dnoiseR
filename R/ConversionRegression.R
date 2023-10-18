

# Computing Quantile Function for joint distribution
#' @export
#'
QuantileFunction <- function(q,cdf, x = seq(length(cdf))){
  i = min(which(cdf >= q))
  return(x[i])
}

# function for creating latent joint distribution.
#' @export
#'
JointDistribution <- function (gamma, zeta, n.quantiles = 10 * max(c(length(gamma),
                                                                     length(zeta))))
{
  q.grid = seq(0, n.quantiles)/n.quantiles
  i.set = rep(NA, length(q.grid))
  j.set = rep(NA, length(q.grid))
  x.set = rep(1/(n.quantiles + 1), length(q.grid))
  cdf.g = cumsum(gamma)
  cdf.g = cdf.g/max(cdf.g)
  cdf.z = cumsum(zeta)
  cdf.z = cdf.z/max(cdf.z)
  R.bins.g = length(gamma)
  R.bins.z = length(zeta)
  for (k in seq(length(q.grid))) {
    q = q.grid[k]
    i.set[k] = QuantileFunction(q, cdf.g)
    j.set[k] = QuantileFunction(q, cdf.z)
    x.set[k] = 1/(n.quantiles + 1)
  }
  p.gz <- Matrix::sparseMatrix(i = i.set, j = j.set, x = x.set,
                               dims = c(R.bins.g, R.bins.z))
  return(p.gz)
}

#' @export
#'
weight_vec <- function(x,X, ker.set){
  if(length(ker.set) != length(x)){
    stop("ker.set and x must have the same length")
  }
  if(ncol(X) != length(x)){
    stop("x and ncol(X) must have the same length")
  }
  w.mat <- matrix(NA,nrow = nrow(X), ncol = ncol(X))
  for(i in seq(length(ker.set))){
    ker.tmp <- ker.set[[i]]
    w.tmp <- ker.tmp(x[i] - X[,i])
    w.mat[,i] <- w.tmp
  }
  w.vec <- exp(rowSums(log(w.mat)))
  #w.vec <- nrow(X)*w.vec/sum(w.vec)
  return(w.vec)
}

#' @export
#'
unique_X <- function(X){
  X.unique <- unique(X, MARGIN = 1)
  return(X.unique)
}


#' Impute Outcomes
#'
#' @param X.ref Covariates
#' @param y.ref Outcome 1
#' @param z.ref Outcome 2 (imputed outcome)
#' @param n.impute Number of times to impute
#' @param Y.train Training data
#' @param Z.train Testing data
#' @param cond.y Conditional distribution 1
#' @param cond.z Conditional distribution 2
#' @param mu.y Tuning parameter 1
#' @param mu.z Tuning parameter 2
#' @param ref.cols reference columns
#' @param ker.set Kernels for smoothing over X parameters
#' @param R.bins Number of bins for model parameter
#' @param B.boot Number of Bootstrap iterations
#'
#' @export
#'
ImputeOutcomes <- function (X.ref, y.ref, z.ref, n.impute, Y.train, Z.train, cond.y,
                            cond.z, mu.y, mu.z, ref.cols, ker.set, R.bins = 1000, threshold = 5 *
                              10^(-5), max.iter = 50, verbose = F, init.latents = F,
                            latent.set.covariates = NA, init.latent.set.y = NA, init.latent.set.z = NA)
{
  if (missing(ref.cols) | missing(ker.set)) {
    if (missing(n.impute)) {
      stop("must indicate number of imputations")
    }
    if (missing(cond.y) | missing(cond.y)) {
      stop("must specify conditional distributions")
    }
    if (missing(mu.y) | missing(mu.z)) {
      stop("must specify regularization parameters")
    }
    if (missing(Y.train) | missing(Z.train)) {
      stop("must specify training data")
    }
    Ny <- length(cond.y(0.5)) - 1
    Nz <- length(cond.z(0.5)) - 1
    A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
    A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
    outcome.cols.y = stringr::str_detect(colnames(Y.train),
                                         "y\\d+") | stringr::str_detect(colnames(Y.train),
                                                                        "y")
    outcome.cols.z = stringr::str_detect(colnames(Z.train),
                                         "z\\d+") | stringr::str_detect(colnames(Z.train),
                                                                        "z")
    y.tr <- Y.train[, outcome.cols.y]
    z.tr <- Z.train[, outcome.cols.z]
    cond.block <- array(NA, c(nrow(X.ref), Ny + 1, Nz + 1))
    cond.block.Boot <- array(NA, c(nrow(X.ref), Ny + 1, Nz +
                                     1))
    p.hat.y <- compute_edf(y.tr, Ny)
    p.hat.z <- compute_edf(z.tr, Nz)
    mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y, mu.y,
                                  threshold = threshold, max.iter = max.iter)
    mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z, mu.z,
                                  threshold = threshold, max.iter = max.iter)
    mixture.y <- mix.y$latent
    mixture.z <- mix.z$latent
    p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,
                                                       cond.z, mixture.y, mixture.z)
    for (id in 1:nrow(X.ref)) {
      cond.block[id, , ] <- p.z.cond.y.slice
    }
    cond.mat <- matrix(NA, nrow(X.ref), Nz + 1)
    for (j in seq(length(y.ref))) {
      if (!is.na(y.ref[j])) {
        cond.mat[j, ] <- cond.block[j, y.ref[j] + 1,
        ]
      }
      else {
        cond.mat[j, ] <- rep(1/(Nz + 1), Nz + 1)
      }
    }
    Z.imp <- matrix(NA, length(y.ref), n.impute)
    for (k in 1:length(y.ref)) {
      if (!is.na(y.ref[k])) {
        p.vec <- as.numeric(cond.mat[k, ])
        Z.imp.row <- sample(seq(0, Nz), n.impute, prob = p.vec,
                            replace = T)
        Z.imp[k, ] <- Z.imp.row
      }
      else {
        Z.imp[k, ] <- rep(NA, n.impute)
      }
    }
    na.idx <- which(is.na(Z.imp[, 1]))
    for (i in na.idx) {
      Z.imp[i, ] <- z.ref[i]
    }
  }
  else {
    if (length(ref.cols) == 1) {
      X.un <- data.frame(tmp = unique((X.ref[, ref.cols])))
    }
    else {
      X.un <- unique_X(X.ref[, ref.cols])
    }
    colnames(X.un) = ref.cols
    if (length(ref.cols) != length(ker.set)) {
      stop("ref.cols and ker.set must be the same length")
    }
    if (missing(n.impute)) {
      stop("must indicate number of imputations")
    }
    if (missing(cond.y) | missing(cond.y)) {
      stop("must specify conditional distributions")
    }
    if (missing(mu.y) | missing(mu.z)) {
      stop("must specify regularization parameters")
    }
    if (missing(Y.train) | missing(Z.train)) {
      stop("must specify training data")
    }
    Ny <- length(cond.y(0.5)) - 1
    Nz <- length(cond.z(0.5)) - 1
    A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
    A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
    outcome.cols.y = stringr::str_detect(colnames(Y.train),
                                         "y\\d+") | stringr::str_detect(colnames(Y.train),
                                                                        "y")
    outcome.cols.z = stringr::str_detect(colnames(Z.train),
                                         "z\\d+") | stringr::str_detect(colnames(Z.train),
                                                                        "z")
    y.tr <- as.matrix(Y.train[, outcome.cols.y])
    z.tr <- as.matrix(Z.train[, outcome.cols.z])
    X.train.Y <- as.data.frame(Y.train[, ref.cols])
    X.train.Z <- as.data.frame(Z.train[, ref.cols])
    colnames(X.train.Y) = ref.cols
    colnames(X.train.Z) = ref.cols
    cond.block <- array(NA, c(nrow(X.ref), Ny + 1, Nz + 1))
    for (i in seq(nrow(X.un))) {
      if (verbose) {
        cat(paste("Unique mixture estimate:", i, "/",
                  nrow(X.un)), end = "\r")
      }
      x <- as.numeric(X.un[i, ])
      weights.y <- weight_vec(x, X.train.Y, ker.set)
      weights.z <- weight_vec(x, X.train.Z, ker.set)
      p.hat.y <- compute_edf(y.tr[, 1], Ny, weights.y)
      p.hat.z <- compute_edf(z.tr[, 1], Nz, weights.z)
      if (init.latents) {
        j = which(!colSums(t(latent.set.covariates) !=
                             x))
        init.latent.y = init.latent.set.y[j, ]
        init.latent.z = init.latent.set.z[j, ]
        mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y,
                                      mu.y, threshold = threshold, max.iter = max.iter,
                                      init.latent = init.latent.y)
        mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z,
                                      mu.z, threshold = threshold, max.iter = max.iter,
                                      init.latent = init.latent.z)
      }
      else {
        mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y,
                                      mu.y, threshold = threshold, max.iter = max.iter)
        mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z,
                                      mu.z, threshold = threshold, max.iter = max.iter)
      }
      mixture.y <- mix.y$latent
      mixture.z <- mix.z$latent
      p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,
                                                         cond.z, mixture.y, mixture.z, R.bins)
      if (length(ref.cols) > 1) {
        match.idx <- which(apply(X.ref[, ref.cols], 1,
                                 function(z) return(all(z == x))))
      }
      else {
        match.idx <- which(X.ref[, ref.cols] == x)
      }
      for (id in match.idx) {
        cond.block[id, , ] <- p.z.cond.y.slice
      }
    }
    cond.mat <- matrix(NA, nrow(X.ref), Nz + 1)
    for (j in seq(length(y.ref))) {
      if (!is.na(y.ref[j])) {
        cond.mat[j, ] <- cond.block[j, y.ref[j] + 1,
        ]
      }
      else {
        cond.mat[j, ] <- rep(1/(Nz + 1), Nz + 1)
      }
    }
    Z.imp <- matrix(NA, length(y.ref), n.impute)
    for (k in 1:length(y.ref)) {
      if (!is.na(y.ref[k])) {
        p.vec <- as.numeric(cond.mat[k, ])
        Z.imp.row <- sample(seq(0, Nz), n.impute, prob = p.vec,
                            replace = T)
        Z.imp[k, ] <- Z.imp.row
      }
      else {
        Z.imp[k, ] <- rep(NA, n.impute)
      }
    }
    na.idx <- which(is.na(Z.imp[, 1]))
    for (i in na.idx) {
      Z.imp[i, ] <- z.ref[i]
    }
  }
  return(Z.imp)
}


#' @export
#'
#TODO: Add to R package
ImputeOutcomesTrueLatent <- function (X.ref, y.ref, z.ref, n.impute, cond.y,cond.z, ref.cols,
                                      mixture.y.set, mixture.z.set)
{
  R.bins.y = ncol(mixture.y.set)
  R.bins.z = ncol(mixture.z.set)
  if (length(ref.cols) == 1) {
    X.un <- data.frame(tmp = unique((X.ref[, ref.cols])))
  }
  else {
    X.un <- unique_X(X.ref[, ref.cols])
  }
  colnames(X.un) = ref.cols
  if (length(ref.cols) != length(ker.set)) {
    stop("ref.cols and ker.set must be the same length")
  }
  if (missing(n.impute)) {
    stop("must indicate number of imputations")
  }
  if (missing(cond.y) | missing(cond.y)) {
    stop("must specify conditional distributions")
  }
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  A.matrix.y <- compute_A_matrix_2(R.bins.y, cond.y)
  A.matrix.z <- compute_A_matrix_2(R.bins.z, cond.z)

  cond.block <- array(NA, c(nrow(X.ref), Ny + 1, Nz + 1))
  for (i in seq(nrow(X.un))) {
    cat(paste("Unique mixture estimate:", i, "/", nrow(X.un)),
        end = "\r")
    x <- as.numeric(X.un[i, ])
    mixture.y <- mixture.y.set[i,]
    mixture.z <- mixture.z.set[i,]
    p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,
                                                       cond.z, mixture.y, mixture.z)
    if (length(ref.cols) > 1) {
      match.idx <- which(apply(X.ref[, ref.cols], 1,
                               function(z) return(all(z == x))))
    }
    else {
      match.idx <- which(X.ref[, ref.cols] == x)
    }
    for (id in match.idx) {
      cond.block[id, , ] <- p.z.cond.y.slice
    }
  }
  cond.mat <- matrix(NA, nrow(X.ref), Nz + 1)
  for (j in seq(length(y.ref))) {
    if (!is.na(y.ref[j])) {
      cond.mat[j, ] <- cond.block[j, y.ref[j] + 1,
      ]
    }
    else {
      cond.mat[j, ] <- rep(1/(Nz + 1), Nz + 1)
    }
  }
  Z.imp <- matrix(NA, length(y.ref), n.impute)
  for (k in 1:length(y.ref)) {
    if (!is.na(y.ref[k])) {
      p.vec <- as.numeric(cond.mat[k, ])
      Z.imp.row <- sample(seq(0, Nz), n.impute, prob = p.vec,
                          replace = T)
      Z.imp[k, ] <- Z.imp.row
    }
    else {
      Z.imp[k, ] <- rep(NA, n.impute)
    }
  }
  na.idx <- which(is.na(Z.imp[, 1]))
  for (i in na.idx) {
    Z.imp[i, ] <- z.ref[i]
  }
  return(Z.imp)
}


#TODO: Can we port this in to the Impute function for consistency of versions
#' Predict Distributions
#'
#' @param X.ref Covariates Reference
#' @param y.ref Observed Data Reference (to convert)
#' @param Y.train Training data
#' @param Z.train Testing data
#' @param cond.y Conditional distribution 1
#' @param cond.z Conditional distribution 2
#' @param mu.y Tuning parameter 1
#' @param mu.z Tuning parameter 2
#' @param ref.cols reference columns
#' @param ker.set Kernels for smoothing over X parameters
#' @param R.bins Number of bins for model parameter
#' @param threshold Threshold for NPEM Iterations
#' @param max.iter Maximum number of iterations in the NPEM algorithm
#' @param verbose Print fitting details
#'
#' @export
#'
predictedDistributions <- function (X.ref, y.ref,  Y.train, Z.train, cond.y,
                                    cond.z, mu.y, mu.z, ref.cols, ker.set, R.bins = 1000, threshold = 5 *
                                      10^(-5), max.iter = 50, verbose = F)
{
  if (missing(ref.cols) | missing(ker.set)) {
    if (missing(cond.y) | missing(cond.y)) {
      stop("must specify conditional distributions")
    }
    if (missing(mu.y) | missing(mu.z)) {
      stop("must specify regularization parameters")
    }
    if (missing(Y.train) | missing(Z.train)) {
      stop("must specify training data")
    }
    Ny <- length(cond.y(0.5)) - 1
    Nz <- length(cond.z(0.5)) - 1
    A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
    A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
    outcome.cols.y = stringr::str_detect(colnames(Y.train),
                                         "y\\d+") | stringr::str_detect(colnames(Y.train),
                                                                        "y")
    outcome.cols.z = stringr::str_detect(colnames(Z.train),
                                         "z\\d+") | stringr::str_detect(colnames(Z.train),
                                                                        "z")
    y.tr <- Y.train[, outcome.cols.y]
    z.tr <- Z.train[, outcome.cols.z]
    cond.block <- array(NA, c(nrow(X.ref), Ny + 1, Nz + 1))
    cond.block.Boot <- array(NA, c(nrow(X.ref), Ny + 1, Nz +
                                     1))
    p.hat.y <- compute_edf(y.tr, Ny)
    p.hat.z <- compute_edf(z.tr, Nz)
    mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y, mu.y,
                                  threshold = threshold, max.iter = max.iter)
    mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z, mu.z,
                                  threshold = threshold, max.iter = max.iter)
    mixture.y <- mix.y$latent
    mixture.z <- mix.z$latent
    p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,
                                                       cond.z, mixture.y, mixture.z)
    for (id in 1:nrow(X.ref)) {
      cond.block[id, , ] <- p.z.cond.y.slice
    }
    cond.mat <- matrix(NA, nrow(X.ref), Nz + 1)
    for (j in seq(length(y.ref))) {
      if (!is.na(y.ref[j])) {
        cond.mat[j, ] <- cond.block[j, y.ref[j] + 1,
        ]
      }
      else {
        cond.mat[j, ] <- rep(1/(Nz + 1), Nz + 1)
      }
    }
    return(cond.mat)
  }
  else {
    if (length(ref.cols) == 1) {
      X.un <- data.frame(tmp = unique((X.ref[, ref.cols])))
    }
    else {
      X.un <- unique_X(X.ref[, ref.cols])
    }
    colnames(X.un) = ref.cols
    if (length(ref.cols) != length(ker.set)) {
      stop("ref.cols and ker.set must be the same length")
    }
    if (missing(cond.y) | missing(cond.y)) {
      stop("must specify conditional distributions")
    }
    if (missing(mu.y) | missing(mu.z)) {
      stop("must specify regularization parameters")
    }
    if (missing(Y.train) | missing(Z.train)) {
      stop("must specify training data")
    }
    Ny <- length(cond.y(0.5)) - 1
    Nz <- length(cond.z(0.5)) - 1
    A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
    A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
    outcome.cols.y = stringr::str_detect(colnames(Y.train),
                                         "y\\d+") | stringr::str_detect(colnames(Y.train),
                                                                        "y")
    outcome.cols.z = stringr::str_detect(colnames(Z.train),
                                         "z\\d+") | stringr::str_detect(colnames(Z.train),
                                                                        "z")
    y.tr <- as.matrix(Y.train[, outcome.cols.y])
    z.tr <- as.matrix(Z.train[, outcome.cols.z])
    X.train.Y <- as.data.frame(Y.train[, ref.cols])
    X.train.Z <- as.data.frame(Z.train[, ref.cols])
    colnames(X.train.Y) = ref.cols
    colnames(X.train.Z) = ref.cols
    cond.block <- array(NA, c(nrow(X.ref), Ny + 1, Nz + 1))
    for (i in seq(nrow(X.un))) {
      if (verbose) {
        cat(paste("Unique mixture estimate:", i, "/",
                  nrow(X.un)), end = "\r")
      }
      x <- as.numeric(X.un[i, ])
      weights.y <- weight_vec(x, X.train.Y, ker.set)
      weights.z <- weight_vec(x, X.train.Z, ker.set)
      p.hat.y <- compute_edf(y.tr[, 1], Ny, weights.y)
      p.hat.z <- compute_edf(z.tr[, 1], Nz, weights.z)

      mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y,
                                    mu.y, threshold = threshold, max.iter = max.iter)
      mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z,
                                    mu.z, threshold = threshold, max.iter = max.iter)

      mixture.y <- mix.y$latent
      mixture.z <- mix.z$latent
      p.z.cond.y.slice <- conditional_imputation_model_2(cond.y,
                                                         cond.z, mixture.y, mixture.z, R.bins)
      if (length(ref.cols) > 1) {
        match.idx <- which(apply(X.ref[, ref.cols], 1,
                                 function(z) return(all(z == x))))
      }
      else {
        match.idx <- which(X.ref[, ref.cols] == x)
      }
      for (id in match.idx) {
        cond.block[id, , ] <- p.z.cond.y.slice
      }
    }
    cond.mat <- matrix(NA, nrow(X.ref), Nz + 1)
    for (j in seq(length(y.ref))) {
      if (!is.na(y.ref[j])) {
        cond.mat[j, ] <- cond.block[j, y.ref[j] + 1,
        ]
      }
      else {
        cond.mat[j, ] <- rep(1/(Nz + 1), Nz + 1)
      }
    }
  }
  return(cond.mat)
}



#' Evaluation of Average Cross Entropy On A Test Set
#'
#' @param p.mat Matrix of predicted probabilities
#' @param z.vec Vector of outcomes (Train or test)
#'
#' @return Normalized Cross Entropy
#' @export
#'
empirical_cross_entropy <- function(p.mat, z.vec){
  out <- 0
  n = nrow(p.mat)
  for(j in seq(n)){
    out <- out + log(p.mat[j,z.vec[j] + 1])
  }
  return(-out/n)
}




#' Naive Conversion Method
#'
#' @export
#'
normalScoreConversionProb <- function(y.vec,muy,sdy,muz,sdz, Nz, n.samp = 20000){

  z.prob <- matrix(NA, nrow = length(y.vec), ncol = Nz + 1)

  for(i in seq(length(y.vec))){
    y <- y.vec[i]
    z.pred <- muz + (sdz/sdy)*(y - muy)
    eps <- rnorm(n.samp, mean = 0, sd = sdz)

    # replace the sampled version with the cdf
    reference.points <- c(seq(1,Nz) - 0.5, Inf)
    cdf.vals <- c(0,pnorm(reference.points, mean = z.pred, sd = sdz))
    out <- rep(0,Nz + 1)
    for(j in seq(length(cdf.vals) - 1)){
      out[j] = cdf.vals[j + 1] - cdf.vals[j]
    }
    z.prob[i,] = out
  }
  return(z.prob)
}






#' @export
#'
conditional_imputation_model_2 <- function (cond.y, cond.z, mixture.y, mixture.z,
                                            n.quantiles = 10 * max(c(length(mixture.y),length(mixture.z))))
{
  R.bins <- length(mixture.y)
  A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
  A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  p.gamma.zeta <- JointDistribution(mixture.y, mixture.z, n.quantiles)
  p.yz <- A.matrix.y %*% p.gamma.zeta %*% t(A.matrix.z)
  p.yz <- p.yz/sum(p.yz)
  p.marg.y <- Matrix::rowSums(p.yz)
  p.z.cond.y <- p.yz/p.marg.y
  p.z.cond.y <- matrix(0, nrow = Ny + 1, ncol = Nz + 1)
  for (k in seq(Ny + 1)) {
    p.z.cond.y[k, ] <- p.yz[k, ]/p.marg.y[k]
  }
  return(p.z.cond.y)
}



# requires an indicator of sequential observations.
# requires #id and #visit variables
# only for a specific application
# X.ref <- clean.tests.lag.3y
# y.ref <- y.ref.3y
# z.ref <- z.ref.3y

#' Impute Outcomes of Differences over time
#'
#' @param X.ref Covariates
#' @param y.ref Outcome 1
#' @param z.ref Outcome 2
#' @param n.impute Number of times to impute
#' @param Y.train Training data
#' @param Z.train Testing data
#' @param cond.y Conditional distribution 1
#' @param cond.z Conditional distribution 2
#' @param mu.y Tuning parameter 1
#' @param mu.z Tuning parameter 2
#' @param ref.cols reference columns
#' @param ker.set Kernels for smoothing over X parameters
#' @param R.bins Number of bins for model parameter
#'
#' @export
#'
ImputeOutcomeDifferences <- function(X.ref,y.ref,z.ref,n.impute,
                                     Y.train,Z.train,cond.y,cond.z,
                                     mu.y,mu.z,ref.cols,ker.set,R.bins = 1000){

  Z.imp <- ImputeOutcomes(X.ref,y.ref,z.ref,n.impute,
                          Y.train,Z.train,cond.y,cond.z,
                          mu.y,mu.z,ref.cols,ker.set,R.bins)

  id.set <- unique(X.ref$id)
  X.out <- data.frame(matrix(NA, nrow = 0, ncol = ncol(X.ref)))
  Z.diff <- matrix(NA, nrow = 0, ncol = ncol(Z.imp))
  colnames(X.out) <- colnames(X.ref)
  k = 1
  cat(end = "\n")
  complete.vec <- c()
  for(id in id.set){

    if(k %% 100 == 0 ){
      cat(paste0("Computing Score Differences: ",k,"/",length(id.set)), end = "\r")
    }

    k = k + 1
    idx <- which(X.ref$id == id & X.ref$visit == 1)
    idx2 <- which(X.ref$id == id & X.ref$visit == 2)
    if(length(idx) != 1){
      break
    }
    if(length(idx2) != 1){
      break
    }
    # fixed difference vector
    Z.diff <- rbind(Z.diff,  Z.imp[idx2,] - Z.imp[idx,])
    X.out <- rbind(X.out, X.ref[idx,])
    z.obs1 <- !is.na(X.ref$z[idx])
    z.obs2 <- !is.na(X.ref$z[idx2])
    if(z.obs1 & z.obs1){
      complete.vec <- c(complete.vec, 1)
    }else {
      complete.vec <- c(complete.vec, 0)
    }

  }
  X.out$complete <- complete.vec
  return(list("X" = X.out, "Z" = Z.diff))
}


#' @export
ZScoreMatchDifferences <- function(X.ref,y.train,z.train, Ny = 30,Nz = 30) {
  mu.y <- mean(y.train$y1, na.rm = T)
  sd.y <- sd(y.train$y1, na.rm = T)

  mu.z <- mean(z.train$z1, na.rm = T)
  sd.z <- sd(z.train$z1, na.rm = T)


  #quantile.vec <- rep(NA, Ny + 1)
  conversion.vec <- rep(NA, Ny + 1)

  for(y in seq(0,Ny)){
    q <- round((sd.z/sd.y)*(y - mu.y) + mu.z)
    conversion.vec[y + 1] = max(min(q,Nz),0)
  }

  X.ref.new <- X.ref
  X.ref.new <- X.ref.new %>% arrange(id,visit)
  idx.set.1 <- 2*(1:(nrow(X.ref.new)/2) ) - 1
  idx.set.2 <- 2*(1:(nrow(X.ref.new)/2) )

  y.out <- X.ref.new$y
  z.out <- X.ref.new$z
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing]+ 1]

  z1.vec <- z.out[idx.set.1]
  z2.vec <- z.out[idx.set.2]
  z.diff <- z2.vec - z1.vec

  X.out <- X.ref.new[idx.set.1,]
  return(list("X" = X.out, "Z" = z.diff))
}

#' @export
ZScoreConversion <- function (X.ref, y.train, z.train, Ny = 30, Nz = 30)
{
  outcome.cols.y = stringr::str_detect(colnames(y.train),
                                       "y\\d+") | stringr::str_detect(colnames(y.train),
                                                                      "y")
  outcome.cols.z = stringr::str_detect(colnames(z.train),
                                       "z\\d+") | stringr::str_detect(colnames(z.train),
                                                                      "z")
  y.tr <- Y.train[, outcome.cols.y]
  z.tr <- Z.train[, outcome.cols.z]

  mu.y <- mean(y.tr, na.rm = T)
  sd.y <- sd(y.tr, na.rm = T)
  mu.z <- mean(z.tr, na.rm = T)
  sd.z <- sd(z.tr, na.rm = T)
  conversion.vec <- rep(NA, Ny + 1)
  for (y in seq(0, Ny)) {
    q <- round((sd.z/sd.y) * (y - mu.y) + mu.z)
    conversion.vec[y + 1] = max(min(q, Nz), 0)
  }
  y.out <- X.ref$Y
  z.out <- X.ref$Z
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing] +
                                         1]
  X.out <- X.ref
  return(list(X = X.out, Z = z.out))
}

#TODO: Add the column name agnostic version to the Match Differences Functions

#' @export
QuantileMatchDifferences <- function(X.ref,y.train,z.train, Ny = 30,Nz = 30) {

  #quantile.vec <- rep(NA, Ny + 1)
  conversion.vec <- rep(NA, Ny + 1)

  for(y in seq(0,Ny)){
    p <- mean(y.train$y1 <= y)
    q <- round(quantile(z.train$z1, p))
    conversion.vec[y + 1] = q
  }

  X.ref.new <- X.ref
  X.ref.new <- X.ref.new %>% arrange(id,visit)
  idx.set.1 <- 2*(1:(nrow(X.ref.new)/2) ) - 1
  idx.set.2 <- 2*(1:(nrow(X.ref.new)/2) )

  y.out <- X.ref.new$y
  z.out <- X.ref.new$z
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing]+ 1]

  z1.vec <- z.out[idx.set.1]
  z2.vec <- z.out[idx.set.2]
  z.diff <- z2.vec - z1.vec

  X.out <- X.ref.new[idx.set.1,]
  return(list("X" = X.out, "Z" = z.diff))
}


#' @export
QuantileConversion <- function (X.ref, y.train, z.train, Ny = 30, Nz = 30)
{
  outcome.cols.y = stringr::str_detect(colnames(y.train),
                                       "y\\d+") | stringr::str_detect(colnames(y.train),
                                                                      "y")
  outcome.cols.z = stringr::str_detect(colnames(z.train),
                                       "z\\d+") | stringr::str_detect(colnames(z.train),
                                                                      "z")
  y.tr <- y.train[, outcome.cols.y]
  z.tr <- z.train[, outcome.cols.z]
  conversion.vec <- rep(NA, Ny + 1)
  y.out <- X.ref$Y
  z.out <- X.ref$Z
  for (y in seq(0, Ny)) {
    p <- mean(y.tr <= y, na.rm = T)
    q <- ceiling(quantile(z.tr, p, na.rm = T))
    conversion.vec[y + 1] = q
  }
  idx.missing <- which(is.na(z.out))
  z.out[idx.missing] <- conversion.vec[y.out[idx.missing] +
                                         1]
  X.out <- X.ref
  return(list(X = X.out, Z = z.out))
}


# should be able to pass in all the info of GEE
# Rely on GEE package
# pass on all information to the GLM function
#' @export
#'
ImputationRegressionGLM <- function (formula, X, Z.impute,fit.cc = T, verbose = F, ...)
{
  outcome = strsplit(format(formula), split = " ")[[1]][1]
  n.impute <- ncol(Z.impute)
  if (!("complete" %in% colnames(X))) {
    stop("Must indicate which rows have complete data")
  }
  n.cc <- nrow(X[X$complete == 1, ])
  if (fit.cc) {
    data.full = cbind(as.numeric(Z.impute[, 1]), X)
    data.full <- data.full[data.full$complete == 1, ]
    n.cc <- nrow(data.full)
    colnames(data.full)[1] <- outcome
    naive.fit <- tryCatch(expr = {
      glm(formula, data = data.full, ...)
    }, error = function(e) {
      message("Error in Complete Cases Fit")
    }, warning = function(w) {
      message("Warning in Complete Cases Fit")
    }, finally = {
    })
  }
  else {
    data.full = cbind(as.numeric(Z.impute[, 1]), X)
    if (!("complete" %in% colnames(data.full))) {
      stop("Must indicate which rows have complete data")
    }
    n.cc <- nrow(data.full)
    colnames(data.full)[1] <- outcome
    naive.fit <- tryCatch(expr = {
      glm(formula, data = data.full, ...)
    }, error = function(e) {
      message("Error in Complete Cases Fit")
    }, warning = function(w) {
      message("Warning in Complete Cases Fit")
    }, finally = {
    })
  }
  n.vars <- length(naive.fit$coefficients)
  Vg <- matrix(0, n.vars, n.vars)
  impute.betas <- matrix(0, n.vars, n.impute)
  for (j in 1:n.impute) {
    data.imp = cbind(as.numeric(Z.impute[, j]), X)
    colnames(data.imp)[1] <- outcome
    imp.fit <- glm(formula, data = data.imp, ...)
    Vg <- Vg + sandwich::sandwich(imp.fit)
    impute.betas[, j] <- imp.fit$coefficients
  }
  Vg <- Vg/n.impute
  beta.impute.mean <- rowMeans(impute.betas)
  names(beta.impute.mean) <- names(imp.fit$coefficients)
  Vb <- matrix(0, n.vars, n.vars)
  mean.centered.betas <- impute.betas - beta.impute.mean
  rownames(mean.centered.betas) <- names(imp.fit$coefficients)
  colnames(Vb) <- names(imp.fit$coefficients)
  rownames(Vb) <- names(imp.fit$coefficients)
  colnames(Vg) <- names(imp.fit$coefficients)
  rownames(Vg) <- names(imp.fit$coefficients)
  for (j in 1:n.impute) {
    Vb <- Vb + outer(mean.centered.betas[, j], mean.centered.betas[,
                                                                   j], "*")
  }
  Vb <- (1/(n.impute - 1)) * Vb
  rubin.var <- Vg + (1 + 1/n.impute) * Vb
  impute.varfrac <- max(diag((((1/n.impute) * Vb)/rubin.var)))
  if(verbose){
    print(paste0("Residual Fraction of Variance due to imputation: ",
                 round(impute.varfrac, 3)))
  }

  z.scores <- beta.impute.mean/sqrt(diag(rubin.var))
  p.vals.imp <- 2 * pnorm(-abs(z.scores))
  if (fit.cc) {
    cc.var <- sandwich::sandwich(naive.fit)
    cc.beta <- naive.fit$coefficients
    cc.z.scores <- cc.beta/sqrt(diag(cc.var))
    cc.p.vals <- 2 * pnorm(-abs(cc.z.scores))
  }
  else {
    cc.var <- NA
    cc.beta <- NA
    cc.z.scores <- NA
    cc.p.vals <- NA
    n.cc <- NA
  }
  out.list <- list("coefficients" = beta.impute.mean,
                   "variance" = rubin.var,
                   "z-scores" = z.scores,
                   "p-values" = p.vals.imp,
                   "cc-coefficients" = cc.beta,
                   "cc-variance" = cc.var,
                   "cc-z-scores" = cc.z.scores,
                   "cc-p-values" = cc.p.vals,
                   "cc-n" = n.cc,
                   "impute.variance.fraction" = impute.varfrac)

  return(out.list)
}

#' @export
#'
ImputationRegressionGLMBootstrap <- function (formula, X.ref, y.ref, z.ref, n.impute, Y.train, Z.train,
                                              cond.y, cond.z, mu.y, mu.z, ref.cols, ker.set, R.bins = 1000,
                                              threshold = 5 * 10^(-5), max.iter = 3, B.boot = 200, verbose = F)
{
  if (length(ref.cols) == 1) {
    X.un <- data.frame(tmp = unique((X.ref[, ref.cols])))
  }
  else {
    X.un <- unique_X(X.ref[, ref.cols])
  }
  colnames(X.un) = ref.cols
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  A.matrix.y <- compute_A_matrix_2(R.bins, cond.y)
  A.matrix.z <- compute_A_matrix_2(R.bins, cond.z)
  outcome.cols.y = stringr::str_detect(colnames(Y.train), "y\\d+") |
    stringr::str_detect(colnames(Y.train), "y")
  outcome.cols.z = stringr::str_detect(colnames(Z.train), "z\\d+") |
    stringr::str_detect(colnames(Z.train), "z")
  y.tr <- as.matrix(Y.train[, outcome.cols.y])
  z.tr <- as.matrix(Z.train[, outcome.cols.z])
  X.train.Y <- as.data.frame(Y.train[, ref.cols])
  X.train.Z <- as.data.frame(Z.train[, ref.cols])
  colnames(X.train.Y) = ref.cols
  colnames(X.train.Z) = ref.cols
  latent.set.covariates = X.un
  init.latent.set.y = matrix(NA, nrow = nrow(X.un), ncol = R.bins)
  init.latent.set.z = matrix(NA, nrow = nrow(X.un), ncol = R.bins)
  for (i in seq(nrow(X.un))) {
    if (verbose) {
      cat(paste("Unique mixture estimate:", i, "/", nrow(X.un)),
          end = "\r")
    }
    x <- as.numeric(X.un[i, ])
    weights.y <- weight_vec(x, X.train.Y, ker.set)
    weights.z <- weight_vec(x, X.train.Z, ker.set)
    p.hat.y <- compute_edf(y.tr[, 1], Ny, weights.y)
    p.hat.z <- compute_edf(z.tr[, 1], Nz, weights.z)
    mix.y <- estimate_mixing_npem(p.hat.y, A.matrix.y, mu.y,
                                  threshold = threshold)
    mix.z <- estimate_mixing_npem(p.hat.z, A.matrix.z, mu.z,
                                  threshold = threshold)
    mixture.y <- mix.y$latent
    mixture.z <- mix.z$latent
    init.latent.set.y[i, ] = mixture.y
    init.latent.set.z[i, ] = mixture.z
  }
  idx.y = seq(nrow(Y.train))
  idx.z = seq(nrow(Z.train))
  idx.ref = seq(length(y.ref))

  X.frame = X.ref
  X.frame$complete = 1 * !is.na(X.frame$Z)
  coef.boot = NULL
  var.boot = NULL
  for (b in seq(B.boot)) {
    boot.idx.y = sample(idx.y, replace = T)
    boot.idx.z = sample(idx.z, replace = T)
    Y.train.boot = Y.train[boot.idx.y, ]
    Z.train.boot = Z.train[boot.idx.z, ]

    # Added the simple bootstrap conversion here
    boot.idx.ref = sample(idx.ref, replace = T)
    X.ref.boot = X.ref[boot.idx.ref,]
    y.ref.boot = y.ref[boot.idx.ref]
    z.ref.boot = z.ref[boot.idx.ref]

    Z.impute <- ImputeOutcomes(X.ref.boot, y.ref.boot, z.ref.boot, n.impute,
                               Y.train.boot, Z.train.boot, cond.y, cond.z, mu.y, mu.z, ref.cols,
                               ker.set, R.bins, verbose = F, max.iter = max.iter,
                               init.latents = TRUE, latent.set.covariates = latent.set.covariates,
                               init.latent.set.y = init.latent.set.y, init.latent.set.z = init.latent.set.z)

    res.tmp <- ImputationRegressionGLM(formula, X.frame,
                                       Z.impute, fit.cc = F)
    if (is.null(coef.boot)) {
      coef.boot <- matrix(NA, ncol = length(res.tmp$coefficients),
                          nrow = B.boot)
    }
    coef.boot[b, ] <- res.tmp$coefficients

    if (is.null(var.boot)) {
      var.boot <- matrix(NA, ncol = length(res.tmp$variance),
                         nrow = B.boot)
    }
    var.boot[b, ] <- res.tmp$variance

    if (verbose) {
      m1 = (round(20 * b/B.boot))
      m2 = 20 - m1
      progress.bar = paste0("|", strrep("=", m1), strrep("-",
                                                         m2), "|")
      cat(paste0("Bootstrap:", b, "/", B.boot, "  ", progress.bar),
          end = "\r")
    }
  }
  beta.est <- colMeans(coef.boot, na.rm = T)
  # total.var.block <- matrix(data = colMeans(var.boot, na.rm = T),
  #                           nrow = length(beta.est), ncol = length(beta.est)) + var(coef.boot, na.rm = T)
  total.var.block <- var(coef.boot,na.rm = T)
  z.scores <- beta.est/sqrt(diag(total.var.block))
  p.vals.imp <- 2 * pnorm(-abs(z.scores))
  out.list <- list(coefficients = beta.est, variance = total.var.block,
                   `z-scores` = z.scores, `p-values` = p.vals.imp)
  return(out.list)
}




#' @export
fit_to_table <- function(fit){
  n = length(fit$coefficients)
  out.table <- matrix(NA,nrow = n, ncol = 8)
  rownames(out.table) = names(fit$coefficients)
  colnames(out.table) = c("C.C. coef",
                          "C.C. sd",
                          "C.C. Zscore",
                          "C.C. p-values",
                          "Imputation coef",
                          "Imputation sd",
                          "Imputation Zscore",
                          "Imputation p-values",)

  out.table[,1] <-  fit$`cc-coefficients`
  out.table[,2] <-  sqrt(diag(fit$`cc-variance`))
  out.table[,3] <-  fit$`cc-z-scores`
  out.table[,4] <-  fit$`cc-p-values`


  out.table[,5] <-  fit$`coefficients`
  out.table[,6] <-  sqrt(diag(fit$`variance`))
  out.table[,7] <-  fit$`z-scores`
  out.table[,8] <-  fit$`p-values`
  return(out.table)
}



#' GEE Regression from imputed data
#'
#' @param formula regression formula
#' @param X covariates
#' @param Z.impute imputed outcomes
#'
# #' @return
#' @export
#'
# #' @examples
ImputationRegressionGEE <- function(formula, X, Z.impute, ...){
  n.impute <- ncol(Z.impute)
  data.full = cbind(as.numeric(Z.impute[,1]),X)
  if(!("complete" %in% colnames(data.full))){
    stop("Must indicate which rows have complete data")
  }
  data.full <- data.full[data.full$complete == 1,]
  colnames(data.full)[1] <- "z"
  naive.fit <- tryCatch(                       # Applying tryCatch
    expr = {                      # Specifying expression
      gee(formula, data = data.full, ...)
      #message("Everything was fine.")
    },

    error = function(e){          # Specifying error message
      message("Error in Complete Cases Fit")
    },

    warning = function(w){        # Specifying warning message
      message("Warning in Complete Cases Fit")
    },

    finally = {                   # Specifying final message
      #message("tryCatch is finished.")
    }
  )


  n.vars <- length(naive.fit$coefficients)



  Vg <- matrix(0,n.vars,n.vars)
  impute.betas <- matrix(0,n.vars,n.impute)
  for(j in 1:n.impute){
    data.imp = cbind(as.numeric(Z.impute[,j]),X)
    colnames(data.imp)[1] <- "z"
    imp.fit <- gee(formula, data = data.imp, ...)
    Vg <- Vg + imp.fit$robust.variance
    impute.betas[,j] <- imp.fit$coefficients
  }
  Vg <- Vg/n.impute
  beta.impute.mean <- rowMeans(impute.betas)
  names(beta.impute.mean) <- names(imp.fit$coefficients)

  Vb <- matrix(0,n.vars,n.vars)
  mean.centered.betas <- impute.betas - beta.impute.mean
  #colnames(mean.centered.betas) <- names(imp.fit$coefficients)

  rownames(mean.centered.betas) <- names(imp.fit$coefficients)
  colnames(Vb) <- names(imp.fit$coefficients)
  rownames(Vb) <- names(imp.fit$coefficients)
  colnames(Vg) <- names(imp.fit$coefficients)
  rownames(Vg) <- names(imp.fit$coefficients)
  for(j in 1:n.impute){
    Vb <- Vb + outer(mean.centered.betas[,j], mean.centered.betas[,j], "*")
  }
  Vb <- (1/(n.impute - 1))*Vb

  rubin.var <- Vg + (1 + 1/n.impute)*Vb
  impute.varfrac <- max(diag((((1/n.impute)*Vb)/rubin.var)))
  print(paste0("Residual Fraction of Variance due to imputation: ", round(impute.varfrac, 3)))
  z.scores <- beta.impute.mean/sqrt(diag(rubin.var))
  p.vals.imp <- 2*pnorm(-abs(z.scores))
  cc.var <- naive.fit$robust.variance
  cc.beta <- naive.fit$coefficients
  cc.z.scores <- cc.beta/sqrt(diag(cc.var))
  cc.p.vals <- 2*pnorm(-abs(cc.z.scores))
  out.list <- list("coefficients" = beta.impute.mean,
                   "variance" = rubin.var,
                   "z-scores" = z.scores,
                   "p-values" = p.vals.imp,
                   "cc-coefficients" = cc.beta,
                   "cc-variance" = cc.var,
                   "cc-z-scores" = cc.z.scores,
                   "cc-p-values" = cc.p.vals,
                   "impute.variance.fraction" = impute.varfrac)

  return(out.list)
}





