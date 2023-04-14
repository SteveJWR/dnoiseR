


# function for creating latent joint distribution.
JointDistribution <- function(gamma,zeta){
  ng = length(gamma)
  nz = length(zeta)
  grid1 <- seq(0,1, length.out = ng)
  grid2 <- seq(0,1, length.out = nz)
  cost = abs(outer(grid1, grid2, "-"))
  trns.res <- transport(gamma,zeta, cost)
  p.gz <- sparseMatrix(i = trns.res$from, j = trns.res$to, x = trns.res$mass)
  return(p.gz) # returns a sparse Matrix data type
}
