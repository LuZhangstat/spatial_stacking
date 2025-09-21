rm(list = ls())
setwd("./Stacking")
library(ggplot2)
library(MASS)
library(spBayes)
library(geoR)
library(rbenchmark)
library("gridExtra")
library(fields)
source("utils.R") 

args <- commandArgs(trailingOnly = TRUE)
input_id <- as.numeric(args[1])  # Assuming input_id is the first argument

cat("id = ", input_id)


# set parameters for simulations
beta <- as.matrix(c(1, 2))
sigma.sq <- 1
tau.sq <- 0.3
nu <- 0.5
phi <- 2

priors <- list(mu_beta = rep(0, 2),
               inv_V_beta = 1/4 * diag(2),
               a_sigma = 2,
               b_sigma = 2)

samplesize_ls  = seq(200, 900, 100)

## set the candidate values for hyper-pars in stacking ##
# decide the bound of the candidate values for phi
range(c(1 / Matern.cor.to.range(0.6 * sqrt(2), 0.5, cor.target=.05),
        1 / Matern.cor.to.range(0.1 * sqrt(2), 0.5, cor.target=.05),
        1 / Matern.cor.to.range(0.6 * sqrt(2), 1.75, cor.target=.05),
        1 / Matern.cor.to.range(0.1 * sqrt(2), 1.75, cor.target=.05)))

nu_grid = c(0.5, 1, 1.5, 1.75)
phi_grid = c(3, 14, 25, 36)   #3.5 to 35.8 #old: c(3, 9, 15, 31) 
deltasq_grid <- pick_deltasq(E_sigmasq = sigma.sq, 
                             E_tausq = tau.sq, b = max(sigma.sq, tau.sq),
                             p_ls = c(0.05, 0.35, 0.65, 0.95))
deltasq_grid

K_fold = 10 # K-fold cross-validation

# pre-allocation for results #
N_list = length(samplesize_ls)
weights_M_LSE = matrix(0, length(deltasq_grid) * length(phi_grid) * 
                         length(nu_grid), N_list)
weights_M_LP = matrix(0, length(deltasq_grid) * length(phi_grid) * 
                        length(nu_grid), N_list)
raw_data <- list() # record raw data
expect_w <- list() # save the weighted latent process 
expect_y <- list() # save the weighted prediction 
DIV_matrix <- matrix(NA, nrow = N_list, ncol = 18)
colnames(DIV_matrix) <- c("SPE_stack_LSE", "SPE_stack_LP", "SPE_M0", 
                          "SPE_MCMC", "SPE_stack_LSE_P", "SPE_stack_LP_P",
                          "ELPD_stack_LSE", "ELPD_stack_LP", "ELPD_M0", 
                          "ELPD_MCMC", "ELPD_stack_LSE_P", "ELPD_stack_LP_P",
                          "SPE_w_stack_LP", "SPE_w_stack_LSE", "SPE_w_M0", 
                          "SPE_w_MCMC", "SPE_w_stack_LP_P", "SPE_w_stack_LSE_P")
rownames(DIV_matrix) <- paste(samplesize_ls) # check
run_time <- matrix(0, 6, ncol = N_list)
MCMC_par <- list() # record the thinned MCMC chains for hyperparameters

for(r in 1:N_list){ # repeat
  cat("\n", "samplesize:", samplesize_ls[r], "\t")
  seed = samplesize_ls[r] + input_id
  set.seed(seed)
  N <- samplesize_ls[r]
  coords <- cbind(runif(N), runif(N))
  X <- as.matrix(cbind(1, rnorm(N)))
  
  D <- as.matrix(dist(coords))
  R <- geoR::matern(D, phi = 1/phi, kappa = nu)
  w <- mvrnorm(n = 1, mu = rep(0, N), Sigma = sigma.sq * R)
  if(is.null(X)){
    y <- rnorm(N, w, sqrt(tau.sq))
  }else{
    y <- rnorm(N, X %*% beta + w, sqrt(tau.sq))
  }
  N_ho = 100
  ind_mod = 1:(N - N_ho)
  
  raw_data[[r]] <- list(X = X, y = y, w = w, coords = coords, beta = beta,
                        phi = phi, nu = nu, tau.sq = tau.sq, 
                        sigma.sq = sigma.sq, ind_mod = ind_mod)
  
  ####################################################################
  ## stacking with prediction (LSE) and log predictive density (LP) ##
  ####################################################################
  CV_fit_LSE <- sp_stacking_K_fold(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
    deltasq_grid = deltasq_grid, phi_grid = phi_grid,
    nu_grid = nu_grid, priors = priors, K_fold = K_fold,
    seed = seed, label = "LSE")
  weights_M_LSE[, r] <- CV_fit_LSE$wts
  run_time[1, r] <- CV_fit_LSE$time[3]
  
  CV_fit_LP <- sp_stacking_K_fold(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
    deltasq_grid = deltasq_grid, phi_grid = phi_grid,
    nu_grid = nu_grid, priors = priors, K_fold = K_fold,
    seed = seed, label = "LP", MC = FALSE)
  weights_M_LP[, r] <- CV_fit_LP$wts
  run_time[2, r] <- CV_fit_LP$time[3]
  
  
  ## stacking mean squared prediction error ##
  y_pred_grid <- matrix(0, nrow = N_ho, ncol = nrow(CV_fit_LSE$grid_all))
  w_expect_grid <- matrix(0, nrow = N, ncol = nrow(CV_fit_LSE$grid_all))
  
  for (i in 1:nrow(CV_fit_LSE$grid_all)){
    if( (CV_fit_LSE$wts[i]>0) | (CV_fit_LP$wts[i]>0)){
      pred_grid <- Conj_predict(X.mod = X[ind_mod, ], y.mod = y[ind_mod],
                                coords.mod = coords[ind_mod, ],
                                deltasq_pick = CV_fit_LSE$grid_all$deltasq[i],
                                phi_pick = CV_fit_LSE$grid_all$phi[i], 
                                nu_pick = CV_fit_LSE$grid_all$nu[i],
                                priors,
                                X.ho = X[-ind_mod, ], 
                                coords.ho = coords[-ind_mod, ])
      y_pred_grid[, i] <- pred_grid$y_expect
      w_expect_grid[, i] <- pred_grid$w_expect
    }
  }
  y_pred_stack_LSE = y_pred_grid %*% CV_fit_LSE$wts
  DIV_matrix[r, "SPE_stack_LSE"] <- mean((y_pred_stack_LSE - y[-ind_mod])^2)
  y_pred_stack_LP = y_pred_grid %*% CV_fit_LP$wts
  DIV_matrix[r, "SPE_stack_LP"] <- mean((y_pred_stack_LP - y[-ind_mod])^2)
  
  w_expect_stack_LSE = w_expect_grid %*% CV_fit_LSE$wts
  DIV_matrix[r, "SPE_w_stack_LSE"] <- mean((w_expect_stack_LSE - w)^2)
  w_expect_stack_LP = w_expect_grid %*% CV_fit_LP$wts
  DIV_matrix[r, "SPE_w_stack_LP"] <- mean((w_expect_stack_LP - w)^2)
  
  
  ## stacking Expected log pointwise predictive density ##
  lp_pred_grid <- matrix(0, nrow = N_ho, ncol = nrow(CV_fit_LSE$grid_all))
  for (i in 1:nrow(CV_fit_LSE$grid_all)){
    if((CV_fit_LSE$wts[i] > 0) | (CV_fit_LP$wts[i] > 0)){
      lp_pred_grid[, i] <- Conj_lpd(X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                                    coords.mod = coords[ind_mod, ], 
                                    deltasq_pick = CV_fit_LSE$grid_all$deltasq[i],
                                    phi_pick = CV_fit_LSE$grid_all$phi[i], 
                                    nu_pick = CV_fit_LSE$grid_all$nu[i],
                                    priors, X.ho = X[-ind_mod, ], 
                                    y.ho = y[-ind_mod], 
                                    coords.ho = coords[-ind_mod, ], MC = FALSE)
    }
  }
  DIV_matrix[r, "ELPD_stack_LSE"] = mean(log(exp(lp_pred_grid) %*% CV_fit_LSE$wts))
  DIV_matrix[r, "ELPD_stack_LP"] = mean(log(exp(lp_pred_grid) %*% CV_fit_LP$wts))
  
  
  ##################################
  ## predict with the exact model ##
  ##################################
  t0 <- proc.time()
  pred_M0 <- Conj_predict(X.mod = X[ind_mod, ], y.mod = y[ind_mod],
                          coords.mod = coords[ind_mod, ],
                          deltasq_pick = tau.sq / sigma.sq,
                          phi_pick = phi, nu_pick = nu, priors,
                          X.ho = X[-ind_mod, ],
                          coords.ho = coords[-ind_mod, ])
  t1 <- proc.time() - t0
  run_time[3, r] = t1[3]
  DIV_matrix[r, "SPE_M0"] <- mean((pred_M0$y_expect - y[-ind_mod])^2)
  DIV_matrix[r, "SPE_w_M0"] <- mean((pred_M0$w_expect - w)^2)
  
  ## exact model Expected log pointwise predictive density ##
  lp_pred_M0 <- Conj_lpd(X.mod = X[ind_mod, ], y.mod = y[ind_mod],
                         coords.mod = coords[ind_mod, ],
                         deltasq_pick = tau.sq / sigma.sq,
                         phi_pick = phi, nu_pick = nu,
                         priors, X.ho = X[-ind_mod, ],
                         y.ho = y[-ind_mod],
                         coords.ho = coords[-ind_mod, ])
  DIV_matrix[r, "ELPD_M0"] <- mean(lp_pred_M0)
  
  
  #######################
  ## predict with MCMC ##
  #######################
  
  ### fit with spBayes ###
  n.samples <- 20000
  starting <- list("phi"=3/0.5, "sigma.sq"=1, "tau.sq"=1, "nu" = 0.5)
  tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1, "nu" = 0.1)
  priors.1 <- list("beta.Norm"=list(rep(0, ncol(X)), solve(priors$inv_V_beta)),
                   "phi.Unif"=c(3, 36), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 2), "nu.unif" = c(0.25, 2))
  cov.model <- "matern"
  n.report <- 5000
  verbose <- TRUE
  m.1 <- spLM(y[ind_mod]~X[ind_mod, ]-1, coords=coords[ind_mod, ],
              starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=n.samples, verbose=verbose, n.report=n.report)
  
  ## recover beta ##
  t0 <- proc.time()
  r.1 <- spRecover(m.1, get.w = FALSE, start = 0.5*n.samples, thin = 10,
                   n.report =  500)
  t1 <- proc.time() - t0
  
  run_time[4, r] <- m.1$run.time[3]
  run_time[5, r] <- t1[3]
  ## recover latent process on all locations ##
  ## compute expected response and latent process ##
  MCMC_out <- expects_MCMC(theta.recover = r.1$p.theta.recover.samples,
                           beta.recover = r.1$p.beta.recover.samples,
                           y.mod = y[ind_mod], X.mod = X[ind_mod, ],
                           coords.mod = coords[ind_mod, ],
                           X.ho = X[-ind_mod, ], y.ho = y[-ind_mod],
                           coords.ho = coords[-ind_mod, ])
  run_time[6, r] <- MCMC_out$time[3]
  MCMC_par[[r]] <- r.1$p.theta.recover.samples
  
  DIV_matrix[r, "SPE_MCMC"] <- mean((MCMC_out$y_expect_MCMC - y[-ind_mod])^2)
  DIV_matrix[r, "SPE_w_MCMC"] <- mean((MCMC_out$w_expect_MCMC - w)^2)
  DIV_matrix[r, "ELPD_MCMC"] <- mean(MCMC_out$lp_expect_MCMC)
  
  expect_y[[r]] <- cbind(y_pred_stack_LSE, y_pred_stack_LP, pred_M0$y_expect,
                         MCMC_out$y_expect_MCMC)
  colnames(expect_y[[r]]) <- c("LSE", "LP", "M0", "MCMC")
  expect_w[[r]] <- cbind(w_expect_stack_LSE, w_expect_stack_LP, pred_M0$w_expect,
                         MCMC_out$w_expect_MCMC)
  colnames(expect_w[[r]]) <- c("LSE", "LP", "M0", "MCMC")
  
  ####################################################################
  ## select prefixed value based on marginal posterior distribution ##
  ####################################################################
  pick_ind <- seq(300, 1000, by = 11)
  all_prefix_ls <- cbind(MCMC_par[[r]][pick_ind, "tau.sq"] / 
                           MCMC_par[[r]][pick_ind, "sigma.sq"],
                         MCMC_par[[r]][pick_ind, c("phi", "nu")]) 
  colnames(all_prefix_ls) <- c("delatsq", "phi", "nu")
  all_prefix_ls <- as.data.frame(all_prefix_ls)
  
  CV_fit_LSE_P <- sp_stacking_K_fold3(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
    all_prefix_ls = all_prefix_ls, priors = priors, 
    K_fold = K_fold, seed = seed, label = "LSE")
  
  CV_fit_LP_P <- sp_stacking_K_fold3(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
    all_prefix_ls = all_prefix_ls, 
    priors = priors, K_fold = K_fold,
    seed = seed, label = "LP", MC = FALSE)
  
  ## stacking mean squared prediction error ##
  y_pred_grid_P <- matrix(0, nrow = N_ho, ncol = nrow(CV_fit_LSE_P$grid_all))
  w_expect_grid_P <- matrix(0, nrow = N, ncol = nrow(CV_fit_LSE_P$grid_all))
  
  for (i in 1:nrow(CV_fit_LSE_P$grid_all)){
    if( (CV_fit_LSE_P$wts[i]>0) | (CV_fit_LP_P$wts[i]>0)){
      pred_grid_P <- Conj_predict(X.mod = X[ind_mod, ], y.mod = y[ind_mod],
                                  coords.mod = coords[ind_mod, ],
                                  deltasq_pick = CV_fit_LSE_P$grid_all$deltasq[i],
                                  phi_pick = CV_fit_LSE_P$grid_all$phi[i], 
                                  nu_pick = CV_fit_LSE_P$grid_all$nu[i],
                                  priors,
                                  X.ho = X[-ind_mod, ], 
                                  coords.ho = coords[-ind_mod, ])
      y_pred_grid_P[, i] <- pred_grid_P$y_expect
      w_expect_grid_P[, i] <- pred_grid_P$w_expect
    }
  }
  y_pred_stack_LSE_P = y_pred_grid_P %*% CV_fit_LSE_P$wts
  DIV_matrix[r, "SPE_stack_LSE_P"] <- mean((y_pred_stack_LSE_P - y[-ind_mod])^2)
  y_pred_stack_LP_P = y_pred_grid_P %*% CV_fit_LP_P$wts
  DIV_matrix[r, "SPE_stack_LP_P"] <- mean((y_pred_stack_LP_P - y[-ind_mod])^2)
  
  w_expect_stack_LSE_P = w_expect_grid_P %*% CV_fit_LSE_P$wts
  DIV_matrix[r, "SPE_w_stack_LSE_P"] <- mean((w_expect_stack_LSE_P - w)^2)
  w_expect_stack_LP_P = w_expect_grid_P %*% CV_fit_LP_P$wts
  DIV_matrix[r, "SPE_w_stack_LP_P"] <- mean((w_expect_stack_LP_P - w)^2)
  
  
  ## stacking Expected log pointwise predictive density ##
  lp_pred_grid_P <- matrix(0, nrow = N_ho, ncol = nrow(CV_fit_LSE_P$grid_all))
  for (i in 1:nrow(CV_fit_LSE_P$grid_all)){
    if((CV_fit_LSE_P$wts[i] > 0) | (CV_fit_LP_P$wts[i] > 0)){
      lp_pred_grid_P[, i] <- Conj_lpd(X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                                      coords.mod = coords[ind_mod, ], 
                                      deltasq_pick = CV_fit_LSE_P$grid_all$deltasq[i],
                                      phi_pick = CV_fit_LSE_P$grid_all$phi[i], 
                                      nu_pick = CV_fit_LSE_P$grid_all$nu[i],
                                      priors, X.ho = X[-ind_mod, ], 
                                      y.ho = y[-ind_mod], 
                                      coords.ho = coords[-ind_mod, ], MC = FALSE)
    }
  }
  DIV_matrix[r, "ELPD_stack_LSE_P"] = mean(log(exp(lp_pred_grid_P) %*% 
                                                 CV_fit_LSE_P$wts))
  DIV_matrix[r, "ELPD_stack_LP_P"] = mean(log(exp(lp_pred_grid_P) %*% 
                                                CV_fit_LP_P$wts))
  
}
summary(DIV_matrix)

output_filename <- paste0("./sim/results/sim3_", input_id, ".Rdata")
save(raw_data, DIV_matrix, run_time, samplesize_ls, weights_M_LSE,
     weights_M_LP, expect_w, expect_y, MCMC_par, file = output_filename)



