rm(list = ls())
#setwd("../../")
library(ggplot2)
library(MASS)
library(spBayes)
library(geoR)
library(rbenchmark)
library("gridExtra")
library(fields)
source("utils2.R") # utils2.R is the testing code #

args <- commandArgs(trailingOnly = TRUE)
input_id <- as.numeric(args[1])  # Assuming input_id is the first argument

cat("id = ", input_id)


## test for the choice of prefixed parameter ##
# test 1: grid #
# decide the bound of the candidate values for phi
range(c(1 / Matern.cor.to.range(0.6 * sqrt(2), 0.5, cor.target=.05),
        1 / Matern.cor.to.range(0.1 * sqrt(2), 0.5, cor.target=.05),
        1 / Matern.cor.to.range(0.6 * sqrt(2), 1.75, cor.target=.05),
        1 / Matern.cor.to.range(0.1 * sqrt(2), 1.75, cor.target=.05)))

deltasq_grid <- c(0.1, 0.5, 1, 2)
phi_grid = c(3, 14, 25, 36)   #3.5 to 35.8 #old: c(3, 9, 15, 31) 
nu_grid = c(0.5, 1, 1.5, 1.75)

# test 2: select delatsq through beta #
deltasq_grid2 <- pick_deltasq(E_sigmasq = 1, E_tausq = 0.3, b = 2,
                              p_ls = c(0.2, 0.4, 0.6, 0.8))
deltasq_grid2

# test 3: posterior sample from marginal distribution #
data_filename <- paste0("./sim_hoffman2/results/sim2_", input_id, ".Rdata")
load(data_filename)

# Simulation #
priors <- list(mu_beta = rep(0, 2),
               inv_V_beta = 1/4 * diag(2),
               a_sigma = 2,
               b_sigma = 2)

K_fold = 10
samplesize_ls  = seq(200, 900, 100)
#samplesize_ls  = seq(200, 400, 100)


N_list = length(samplesize_ls)

weights_M_LSE = matrix(0, length(deltasq_grid) * length(phi_grid) * 
                         length(nu_grid), N_list)
weights_M_LP = matrix(0, length(deltasq_grid) * length(phi_grid) * 
                        length(nu_grid), N_list)
#raw_data <- list() # record raw data
expect_w <- list() # save the weighted latent process 
expect_y <- list() # save the weighted prediction 
DIV_matrix2 <- matrix(NA, nrow = N_list, ncol = 12)
#G: grid; E: empirical; P: posterior
colnames(DIV_matrix2) <- c("SPE_stack_LSE_2", "SPE_stack_LP_2",
                           "SPE_stack_LSE_P", "SPE_stack_LP_P",
                           "ELPD_stack_LSE_2", "ELPD_stack_LP_2",
                           "ELPD_stack_LSE_P", "ELPD_stack_LP_P",
                           "SPE_w_stack_LSE_2", "SPE_w_stack_LP_2",
                           "SPE_w_stack_LSE_P", "SPE_w_stack_LP_P")
rownames(DIV_matrix2) <- paste(samplesize_ls) # check
run_time2 <- matrix(0, 4, ncol = N_list)
rownames(run_time2) <- c("Stack_LSE_2", "Stack_LP_2",
                         "Stack_LSE_P", "Stack_LP_P")
#MCMC_par <- list() # record the thinned MCMC chains for hyperparameters

t <- proc.time()
for(r in 1:N_list){ # repeat
  cat("\n", "samplesize:", samplesize_ls[r], "\t")
  seed = samplesize_ls[r] + input_id
  set.seed(seed)
  
  ind_mod = raw_data[[r]]$ind_mod
  X <- raw_data[[r]]$X
  y <- raw_data[[r]]$y
  w <- raw_data[[r]]$w
  coords <- raw_data[[r]]$coords
  N <- samplesize_ls[r]
  N_ho <- N - length(raw_data[[r]]$ind_mod)
  ##########################################################
  ## select prefixed value based on empirical variogram   ##
  ##########################################################
  CV_fit_LSE <- sp_stacking_K_fold(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ], 
    deltasq_grid = deltasq_grid2, phi_grid = phi_grid, 
    nu_grid = nu_grid, priors = priors, 
    K_fold = K_fold, seed = seed, label = "LSE")
  weights_M_LSE[, r] <- CV_fit_LSE$wts
  run_time2["Stack_LSE_2", r] <- CV_fit_LSE$time[3]
  
  CV_fit_LP <- sp_stacking_K_fold(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
    deltasq_grid = deltasq_grid2, phi_grid = phi_grid, 
    nu_grid = nu_grid, priors = priors, K_fold = K_fold,
    seed = seed, label = "LP", MC = FALSE)
  weights_M_LP[, r] <- CV_fit_LP$wts
  run_time2["Stack_LP_2", r] <- CV_fit_LP$time[3]
  
  
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
  DIV_matrix2[r, "SPE_stack_LSE_2"] <- mean((y_pred_stack_LSE - y[-ind_mod])^2)
  y_pred_stack_LP = y_pred_grid %*% CV_fit_LP$wts
  DIV_matrix2[r, "SPE_stack_LP_2"] <- mean((y_pred_stack_LP - y[-ind_mod])^2)
  
  w_expect_stack_LSE = w_expect_grid %*% CV_fit_LSE$wts
  DIV_matrix2[r, "SPE_w_stack_LSE_2"] <- mean((w_expect_stack_LSE - w)^2)
  w_expect_stack_LP = w_expect_grid %*% CV_fit_LP$wts
  DIV_matrix2[r, "SPE_w_stack_LP_2"] <- mean((w_expect_stack_LP - w)^2)
  
  
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
  DIV_matrix2[r, "ELPD_stack_LSE_2"] = 
    mean(log(exp(lp_pred_grid) %*% CV_fit_LSE$wts))
  DIV_matrix2[r, "ELPD_stack_LP_2"] = 
    mean(log(exp(lp_pred_grid) %*% CV_fit_LP$wts))
  
  
  ####################################################################
  ## select prefixed value based on marginal posterior distribution ##
  ####################################################################
  pick_ind <- seq(300, 1000, by = 11)
  all_prefix_ls <- cbind(MCMC_par[[r]][pick_ind, "tau.sq"] / 
                           MCMC_par[[r]][pick_ind, "sigma.sq"],
                         MCMC_par[[r]][pick_ind, c("phi", "nu")]) 
  colnames(all_prefix_ls) <- c("delatsq", "phi", "nu")
  all_prefix_ls <- as.data.frame(all_prefix_ls)
  
  CV_fit_LSE <- sp_stacking_K_fold3(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
    all_prefix_ls = all_prefix_ls, priors = priors, 
    K_fold = K_fold, seed = seed, label = "LSE")
  weights_M_LSE[, r] <- CV_fit_LSE$wts
  run_time2["Stack_LSE_P", r] <- CV_fit_LSE$time[3]
  
  CV_fit_LP <- sp_stacking_K_fold3(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
    all_prefix_ls = all_prefix_ls, 
    priors = priors, K_fold = K_fold,
    seed = seed, label = "LP", MC = FALSE)
  weights_M_LP[, r] <- CV_fit_LP$wts
  run_time2["Stack_LP_P", r] <- CV_fit_LP$time[3]
  
  
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
  DIV_matrix2[r, "SPE_stack_LSE_P"] <- mean((y_pred_stack_LSE - y[-ind_mod])^2)
  y_pred_stack_LP = y_pred_grid %*% CV_fit_LP$wts
  DIV_matrix2[r, "SPE_stack_LP_P"] <- mean((y_pred_stack_LP - y[-ind_mod])^2)
  
  w_expect_stack_LSE = w_expect_grid %*% CV_fit_LSE$wts
  DIV_matrix2[r, "SPE_w_stack_LSE_P"] <- mean((w_expect_stack_LSE - w)^2)
  w_expect_stack_LP = w_expect_grid %*% CV_fit_LP$wts
  DIV_matrix2[r, "SPE_w_stack_LP_P"] <- mean((w_expect_stack_LP - w)^2)
  
  
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
  DIV_matrix2[r, "ELPD_stack_LSE_P"] = mean(log(exp(lp_pred_grid) %*% CV_fit_LSE$wts))
  DIV_matrix2[r, "ELPD_stack_LP_P"] = mean(log(exp(lp_pred_grid) %*% CV_fit_LP$wts))
  
}
proc.time() - t
summary(DIV_matrix2)

output_filename <- paste0("./sim_carc/results/sim2_", input_id, "_prefix2.Rdata")
save(DIV_matrix2, file = output_filename)


