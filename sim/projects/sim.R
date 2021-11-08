rm(list = ls())
source("utils.R")
load("./sim/data/simdata.RData")


## fit conjugate spatial model ##
alpha_grid = c(0.25, 0.5, 0.75)
D_max = sqrt(2)
phi_grid = 3/(sqrt(2)*c(0.25, 0.5, 0.75))
priors <- list(mu_beta = rep(0, length(beta)),
               inv_V_beta = 1/4 * diag(length(beta)),
               a_sigma = 2,
               b_sigma = 1)
L = 100
