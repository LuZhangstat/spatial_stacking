## check the distribution of inference ##
rm(list = ls())
library(ggplot2)
library(MASS)
library(spBayes)
library(ggplot2)
library("gridExtra")
library("bayesplot")
library("coda")

# colorblind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim_ind = 2 # simulation index 1 or 2
load(paste0("./sim_hoffman2/results/sim", sim_ind, "_1.Rdata"))
priors <- list(mu_beta = rep(0, 2),
               inv_V_beta = 1/4 * diag(2),
               a_sigma = 2,
               b_sigma = 2)

### fit with spBayes ###
r = 3 # rerun the simulation with N = 400
set.seed(r+1)
N_ho = 100
ind_mod = 1:(length(raw_data[[r]]$y) - N_ho)
n.samples <- 40000

fit_flag = FALSE
if(fit_flag){
  starting <- list("phi"=3/0.5, "sigma.sq"=1, "tau.sq"=1, "nu" = 0.5)
  tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1, "nu" = 0.1)
  priors.1 <- list("beta.Norm"=list(rep(0, ncol(raw_data[[r]]$X)), 
                                    solve(priors$inv_V_beta)),
                   "phi.Unif"=c(3, 21), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 2), "nu.unif" = c(0.25, 2))
  cov.model <- "matern"
  n.report <- 5000
  verbose <- TRUE
  t <- proc.time()
  m.1 <- spLM(raw_data[[r]]$y[ind_mod]~raw_data[[r]]$X[ind_mod, ]-1, 
              coords=raw_data[[r]]$coords[ind_mod, ],
              starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=n.samples, verbose=verbose, n.report=n.report)
  proc.time() - t
  save_filename <- 
    paste0("./sim_hoffman2/results/sim", sim_ind, "_MCMCfit.Rdata")
  save(m.1, file = save_filename)
}else{
  load(paste0("./sim_hoffman2/results/sim", sim_ind, "_MCMCfit.Rdata"))
}

## load the results of the simulations ##
K_fold = 10
N_list = length(samplesize_ls)
N_sim = 60
## recover the estimated hyper-parameters ##
deltasq_grid <- c(0.1, 0.5, 1, 2)
phi_grid = c(3, 9, 15, 21)   #3/(0.6*sqrt(2)) to 3/(0.1*sqrt(2)) phi_grid = c(3, 13, 23, 33)
nu_grid = c(0.5, 1, 1.5, 1.75)
grid_all <- expand.grid(deltasq_grid, phi_grid, nu_grid)
colnames(grid_all) <- c("deltasq", "phi", "nu")
weights_M_LSE_all = array(0, dim = c(length(deltasq_grid) * length(phi_grid) * 
                                       length(nu_grid), N_list, N_sim))
weights_M_LP_all = array(0, dim = c(length(deltasq_grid) * length(phi_grid) * 
                                      length(nu_grid), N_list, N_sim))


for(i in 1:N_sim){
  filename <- paste0("./sim_hoffman2/results/sim", sim_ind, "_", i, ".Rdata")
  load(filename)
  weights_M_LSE_all[, , i] <- weights_M_LSE
  weights_M_LP_all[, , i] <- weights_M_LP
}

## check the inference of hyper-parameters ##
expect_w_phi_LSE <- matrix(0, nrow = length(phi_grid), ncol = N_sim * N_list)
expect_w_phi_LP <- matrix(0, nrow = length(phi_grid), ncol = N_sim * N_list)
expect_w_nu_LSE <- matrix(0, nrow = length(nu_grid), ncol = N_sim * N_list)
expect_w_nu_LP <- matrix(0, nrow = length(nu_grid), ncol = N_sim * N_list)
expect_w_deltasq_LSE <- matrix(0, nrow = length(deltasq_grid), ncol = N_sim * N_list)
expect_w_deltasq_LP <- matrix(0, nrow = length(deltasq_grid), ncol = N_sim * N_list)


for(i in 1:length(phi_grid)){
  expect_w_phi_LSE[i, ] <- c(apply(
    weights_M_LSE_all[which(grid_all$phi == phi_grid[i]), , ], 2:3, sum))
  expect_w_phi_LP[i, ] <- c(apply(
    weights_M_LP_all[which(grid_all$phi == phi_grid[i]), , ], 2:3, sum))
  expect_w_nu_LSE[i, ] <- c(apply(
    weights_M_LSE_all[which(grid_all$nu == nu_grid[i]), , ], 2:3, sum))
  expect_w_nu_LP[i, ] <- c(apply(
    weights_M_LP_all[which(grid_all$nu == nu_grid[i]), , ], 2:3, sum))
  expect_w_deltasq_LSE[i, ] <- c(apply(
    weights_M_LSE_all[which(grid_all$deltasq == deltasq_grid[i]), , ], 
    2:3, sum))
  expect_w_deltasq_LP[i, ] <- c(apply(
    weights_M_LP_all[which(grid_all$deltasq == deltasq_grid[i]), , ], 
    2:3, sum))
}

est_phi_LSE <- c(t(expect_w_phi_LSE) %*% phi_grid)
est_phi_LP <- c(t(expect_w_phi_LP) %*% phi_grid)
est_nu_LSE <- c(t(expect_w_nu_LSE) %*% nu_grid)
est_nu_LP <- c(t(expect_w_nu_LP) %*% nu_grid)
est_deltasq_LSE <- c(t(expect_w_deltasq_LSE) %*% deltasq_grid)
est_deltasq_LP <- c(t(expect_w_deltasq_LP) %*% deltasq_grid)

pick_MCMC_ind <- seq(10001, n.samples, by = 10)
str(m.1$p.theta.samples)
mcmc_trace(m.1$p.theta.samples[pick_MCMC_ind,], pars = "phi") 
effectiveSize(m.1$p.theta.samples[pick_MCMC_ind, "phi"]) 
mcmc_trace(m.1$p.theta.samples[pick_MCMC_ind,], pars = "nu") 
effectiveSize(m.1$p.theta.samples[pick_MCMC_ind, "nu"]) 
mcmc_trace(m.1$p.theta.samples[pick_MCMC_ind,], pars = "sigma.sq") 
effectiveSize(m.1$p.theta.samples[pick_MCMC_ind, "sigma.sq"]) 
mcmc_trace(m.1$p.theta.samples, pars = "tau.sq") 
effectiveSize(m.1$p.theta.samples[pick_MCMC_ind, "tau.sq"]) 
median(m.1$p.theta.samples[pick_MCMC_ind, "phi"])


# compare the histograms 
pick_est_ind <- seq(r, 472 + r, by = 8)
infer_compar = data.frame(
  value = c(m.1$p.theta.samples[pick_MCMC_ind, "phi"], est_phi_LSE[pick_est_ind], 
            est_phi_LP[pick_est_ind], m.1$p.theta.samples[pick_MCMC_ind, "nu"], 
            est_nu_LSE[pick_est_ind], est_nu_LP[pick_est_ind], 
            (m.1$p.theta.samples[pick_MCMC_ind, "tau.sq"] / 
               m.1$p.theta.samples[pick_MCMC_ind, "sigma.sq"]), 
            est_deltasq_LSE[pick_est_ind], est_deltasq_LP[pick_est_ind]),
  type = rep(c(rep(1, length(pick_MCMC_ind)), rep(2, 60), rep(3, 60)), 3),
  par = rep(1:3, each = (60*2 + length(pick_MCMC_ind))))
infer_compar$type <- factor(infer_compar$type, levels = c(1, 2, 3), 
                            labels = c("MCMC", "stacking of means", 
                                       "stacking of predictive densities"))
infer_compar$par <- factor(infer_compar$par, levels = c(1, 2, 3), 
                           labels = c(expression(phi), expression(nu),
                                      expression(delta^2)))
vline.data <- data.frame(z = c(raw_data[[r]]$phi, raw_data[[r]]$nu, 
                               raw_data[[r]]$tau.sq / raw_data[[r]]$sigma.sq),
                         par = levels(infer_compar$par))

inf_compar_plot =  ggplot(infer_compar, aes(x=value, fill=type)) + 
  geom_density(color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080", "#E69F00")) +
  #theme_bw(base_size = 10) + 
  geom_vline(aes(xintercept = z), data=vline.data, colour = "red", 
             linetype="dashed") +
  facet_wrap(~par, nrow = 1, scales = "free", labeller = label_parsed) +
  labs(fill="") + xlab("") + 
  # theme(legend.position = c(0.7, 0.8), 
  #       legend.title = element_blank(), legend.background = element_blank())
  theme(legend.position="bottom")

inf_compar_plot

ggsave(paste0("./sim_hoffman2/pics/inf_compar_sim", sim_ind, ".png"), 
       plot = inf_compar_plot, 
       width = 6.5, height = 3, units = "in", dpi = 600)

