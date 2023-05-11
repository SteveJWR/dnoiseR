

#' Cross validation for regression and regularization parameters
#'
#' @param data Data used to estimate the model
#' @param cond Conditional Model
#' @param folds Number of folds
#' @param outcome Which column is the outcome
#' @param X.cols Which columns are covariates used
#' @param ker.params Parameters to input to scale each parameter
#' @param reg.kers Kernels used for smoothing across covariates
#' @param mu.params Parameter vector for the regularization term
#' @param threshold Threshold for EM algorithm likelihood change
#' @param max.iter Maximum number of iterations for EM algorithm
#'
#' @return List of cross validation likelihood, and optimal parameters
#' @export
#'
cv_regression_gridsearch <- function(data,cond,folds = 5, outcome = 1, X.cols = seq(2,ncol(data)), ker.params, reg.kers, mu.params, threshold = 10**(-5), max.iter = 50){
  if(length(ker.params) != length(X.cols)){
    stop("there must be a kernel parameter set which matches the ")
  }
  K = length(ker.params)
  all.params = ker.params
  all.params[[K + 1]] = mu.params
  hyperparam.grid = expand.grid(all.params)
  J = nrow(hyperparam.grid)
  n = nrow(data)
  fold.idx <- caret::createFolds(seq(n), k = folds)
  res = rep(NA,J)
  for(j in seq(J)){
    ker.set = list()
    for(k in seq(K)){
      ker.set[[k]] = scale_kernel(reg.kers[[k]], hyperparam.grid[j,k])
    }
    mu = hyperparam.grid[j,K + 1]
    cat(paste0("Parameter combination: ", j, "/", J), end = "\r")
    tmp = cv_regression_eval(data, cond, outcome, X.cols, fold.idx, ker.set, mu)
    res[j] = tmp
  }
  j.max = which.max(res)
  out.list = list("cv.lik" = res, "parameters" = hyperparam.grid, "opt.lik" = res[j.max], "opt.parameters" =hyperparam.grid[j.max,] )
  return(out.list)
}

#' @export
cv_regression_eval <- function(data, cond, outcome = 1, X.cols = seq(2,ncol(data)),
                               fold.idx, ker.set, mu, verbose = F){
  K = length(fold.idx)
  X.ref <- data[,X.cols]
  y <- data[,outcome]
  n = nrow(data)
  if (length(X.cols) == 1) {
    X.ref <- matrix(X.ref, ncol = 1)
    X.un <- data.frame(tmp = unique((data[, X.cols])))
  }
  else {
    X.un <- unique_X(data[, X.cols])
  }
  if(all(is.character(X.cols))){
    colnames(X.un) = X.cols
  }

  if (length(X.cols) != length(ker.set)) {
    stop("X.cols and ker.set must be the same length")
  }

  if (missing(cond)) {
    stop("must specify conditional distributions")
  }
  if (missing(mu)) {
    stop("must specify regularization parameters")
  }
  if (missing(data)) {
    stop("must specify training data")
  }

  N <- length(cond.y(0.5)) - 1
  A.matrix<- compute_A_matrix_2(R.bins, cond)


  res.vec = rep(NA, K)

  for(k in seq(K)){
    test.idx <- fold.idx[[k]]
    train.idx <- setdiff(seq(n), test.idx)
    y.train <- as.matrix(data[train.idx, outcome])
    X.train <- as.data.frame(data[train.idx, X.cols])

    mixture.block <- matrix(NA, nrow = nrow(X.un), ncol = R.bins)
    for (i in seq(nrow(X.un))) {
      if (verbose) {
        cat(paste("Unique mixture estimate:", i, "/",
                  nrow(X.un)), end = "\r")
      }
      x <- as.numeric(X.un[i, ])
      weights <- weight_vec(x, X.train, ker.set)
      p.hat.train <- compute_edf(y.train, Ny, weights)

      mix.model <- estimate_mixing_npem(p.hat.train, A.matrix,
                                        mu, threshold = threshold,
                                        max.iter = max.iter)
      mixture <- mix.model$latent

      mixture.block[i,] = mixture
    }
    cat("")
    ce = 0

    for (i in seq(nrow(X.un))) {
      if (verbose) {
        cat(paste("Unique mixture estimate:", i, "/",
                  nrow(X.un)), end = "\r")
      }
      x <- as.numeric(X.un[i, ])


      mixture = mixture.block[i,]

      if (length(X.cols) > 1) {
        match.id <- which(apply(data[,X.cols], 1,
                                function(z) return(all(z == x))))
      }
      else {
        match.id <- which(data[,X.cols] == x)
      }
      counts = rep(0,N + 1)
      match.id <- intersect(match.id,test.idx)
      for(id in match.id){
        counts[y[id] + 1] = counts[y[id] + 1] + 1
      }
      p.ma <- as.numeric(A.matrix %*% mixture)
      ce = ce + sum(log(p.ma)*counts)
    }
    res.vec[k] = ce
  }
  return(mean(res.vec))
}
