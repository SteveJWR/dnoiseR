### Basic Math functions which will be used in the package.  We need not document theselogit <- function(x){


logit <- function(x){
  if ({
    any(x < 0)
  } || {
    any(x > 1)
  })
    stop("x must be in [0,1].")
  out = log(x/(1-x))
  return(out)
}

logistic <- function(x){
  return(1/(1 + exp(-x)))
}



kl_divergence <- function(p,q){
  if(!(all(dim(p) == dim(q)))){
    stop("p,q must be of the same dimensions" )
  }
  out <- p*log(p/q)
  out[is.nan(out)] <- 0
  return(sum(out))
}


# TODO: Is this necessary
tv_norm <- function(x,y){
  n <- min(min(x), min(y))
  m <- max(max(x), max(y))

  Nx <- length(x)
  Ny <- length(y)

  series_term <- 0
  for(k in n:m){
    series_term <- series_term + abs((sum(x == k)/Nx) - (sum(y == k)/Ny))
  }
  out <- (1/2)*series_term
  return(out)
}



# Function: compute empirical distribution function
# Input: N, (support size {0,...,N}); integer
#        x, sample of observations; vector of integers
# Output: Empirical distribution function for x
compute_edf <- function(x,N){
  p.hat <- c()
  for(y in 0:(N)){
    prop <- mean(x == y)
    p.hat <- c(p.hat, prop)
  }
  return(p.hat)
}




# all of these are simple functions which can be used to define a kernel function

gaussian_kernel <- function(x){
  return(exp(-x^2/2))
}

exponential_kernel <- function(x){
  return(exp(-abs(x)))
}

logistic_kernel <- function(x){
  return(1/(exp(x) + 2 + exp(-x)))
}


## Kernel support  |x| <= 1
triangle_kernel <- function(x){
  out <- (1 - abs(x))*(abs(x) <= 1)
  return(out)
}


epanechnikov_kernel <- function(x){
  out <- (3/4)*(1 - x^2)*(abs(x) <= 1)
  return(out)
}

### Not Continuous which may cause problems
uniform_kernel <- function(x){
  out <- 1*(abs(x) <= 1)
  return(out)
}

quartic_kernel <- function(x){
  out <- (15/16)*(1 - x^2)^2*(abs(x) <= 1)
  return(out)
}

triweight_kernel <- function(x){
  out <- (35/32)*(1 - x^2)^3*(abs(x) <= 1)
  return(out)
}

tricube_kernel <- function(x){
  out <- (70/81)*(1 - abs(x)^3)^3*(abs(x) <= 1)
  return(out)
}

cosine_kernel <- function(x){
  out <- (pi/4)*(cos(pi*x/2))^3*(abs(x) <= 1)
  return(out)
}


