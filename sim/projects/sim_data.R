rm(list = ls())
library(MASS)

#------------------------- generate simulation data ---------------------------#
set.seed(12)
N <- 200
coords <- cbind(runif(N), runif(N))
X <- NULL
beta <- NULL
# X <- as.matrix(cbind(1, rnorm(N)))
# beta <- as.matrix(c(1, 2))
sigma.sq <- 1
tau.sq <- 0.5
#phi <- 3 / 0.5
phi <- 3

D <- as.matrix(dist(coords))
R <- exp(- phi * D)
w <- mvrnorm(n = 1, mu = rep(0, N), Sigma = sigma.sq * R)
if(is.null(X)){
  y <- rnorm(N, w, sqrt(tau.sq))
}else{
  y <- rnorm(N, X %*% beta + w, sqrt(tau.sq))
}


save(list = c("X", "y", "w", "beta", "sigma.sq", "tau.sq", "phi", 
              "N", "coords"), 
     file = "./sim/data/simdata.RData", 
     envir = .GlobalEnv)

## --- plot the location --- ##
plot(coords, xlim = c(0, 1), ylim = c(0, 1), xlab = " ", ylab = " ")

######################################################################
########################## no Xbeta ##################################
######################################################################


set.seed(12)
N <- 250
coords <- cbind(runif(N), runif(N))
X <- NULL
beta <- NULL
sigma.sq <- 1
tau.sq <- 1
#phi <- 3 / 0.5
phi <- 7.1

D <- as.matrix(dist(coords))
R <- exp(- phi * D)
w <- mvrnorm(n = 1, mu = rep(0, N), Sigma = sigma.sq * R)
y <- rnorm(N, w, sqrt(tau.sq))

save(list = c("X", "y", "w", "beta", "sigma.sq", "tau.sq", "phi", 
              "N", "coords"), 
     file = "./sim/data/simdata.RData", 
     envir = .GlobalEnv)

