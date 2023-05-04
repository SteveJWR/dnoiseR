library(CVXR)
library(transport)
library(Matrix)
n1 = 1000
n2 = 2000
gamma <- runif(n1)
zeta <- runif(n2)
gamma = gamma/sum(gamma)
zeta = zeta/sum(zeta)

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





library(reshape2)
library(ggplot2)
p.gz.mat <- as.matrix(p.gz)
longData<-melt(p.gz.mat)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))



library(dnoiseR)
#library(CVXR)
N = 31
p.hat <- seq(0,1, length.out = N + 1)
p.hat <- rep(c(0,1), (N + 1)/2)
#p.hat <- rep(1,length.out = N + 1)
#p.hat <- c(p.hat,0)
p.hat <- p.hat/sum(p.hat)
plot(p.hat)
#cond = dnoiseR::generate_cond_binomial(N)
cond = dnoiseR::conditional_mkm(N,gaussian_kernel,h = 0.1)
R_bins = 1000
mu = 0.001
inner.grid = 25
threshold = 0.0005


A.matrix <- compute_A_matrix(R_bins, cond, numeric.points = 3)


time1 <- Sys.time()
res1 <- dnoiseR::estimate_mixing_numeric(p.hat, A.matrix, mu, cvx.solver = "SCS")
time2 <- Sys.time()
res2 <- estimate_mixing_npem(p.hat, A.matrix, mu, threshold = threshold)
time3 <- Sys.time()

time3 - time2
time2 - time1

plot(res1$latent)
plot(res2$latent)


like.new <- -kl_divergence(p.hat.1,p.ma.new1)  -kl_divergence(p.hat.2,p.ma.new2) - mu* kl_divergence(uniform.latent, latent.trait.new)


#library(CVXR)
N = 31
obs.set = c(rep(1,50), rep(2,50))
cond = dnoiseR::generate_cond_binomial(N)
gamma.set <- rep(c(0.3,0.7),50)
Y <- simulate_test_cond(obs.set, cond, gamma.set)
#p.hat <- rep(1,length.out = N + 1)
#p.hat <- c(p.hat,0)


#cond = dnoiseR::conditional_mkm(N,gaussian_kernel,h = 0.1)
R_bins = 500
mu = 0.1
inner.grid = 25
threshold = 10**(-6)


A.matrix <- compute_A_matrix(R_bins, cond, numeric.points = 3)
A.tensor <- compute_A_tensor(R_bins, cond, numeric.points = 3)

time1 <- Sys.time()
res1 <- estimate_mixing_numeric_2(Y, A.matrix, A.tensor, mu, cvx.solver = "SCS")
time2 <- Sys.time()
res2 <- estimate_mixing_npem_2(Y, A.matrix, A.tensor,  weights = rep(1,nrow(Y)), mu = mu, threshold = threshold, verbose = T, momentum = 0)
time3 <- Sys.time()

plot(res1$latent)
plot(res2$latent)

time3 - time2
time2 - time1




rm(list = ls())





