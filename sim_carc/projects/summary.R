rm(list = ls())
library(ggplot2)
library(MASS)
library(spBayes)
library(ggplot2)
library("gridExtra")
library("coda")
source("utils.R")
# colorblind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim_ind = 4 # simulation index 1, 2, 3 or 4


load(paste0("./sim_carc/results/sim", sim_ind, "_1.Rdata"))
K_fold = 10
N_list = length(samplesize_ls)
N_sim = 60
SPE_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_M0 = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_M0 = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_M0 = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)

## recover the estimated hyper-parameters ##
deltasq_grid <- pick_deltasq(E_sigmasq = raw_data[[1]]$sigma.sq, 
                             E_tausq = raw_data[[1]]$tau.sq, 
                             b = max(raw_data[[1]]$sigma.sq, 
                                     raw_data[[1]]$tau.sq),
                              p_ls = c(0.05, 0.35, 0.65, 0.95))
deltasq_grid
phi_grid = c(3, 14, 25, 36)   #c(3, 9, 15, 21) 3/(0.6*sqrt(2)) to 3/(0.1*sqrt(2)) phi_grid = c(3, 13, 23, 33)
nu_grid = c(0.5, 1, 1.5, 1.75)
grid_all <- expand.grid(deltasq_grid, phi_grid, nu_grid)
colnames(grid_all) <- c("deltasq", "phi", "nu")
weights_M_LSE_all = array(0, dim = c(length(deltasq_grid) * length(phi_grid) * 
                                       length(nu_grid), N_list, N_sim))
weights_M_LP_all = array(0, dim = c(length(deltasq_grid) * length(phi_grid) * 
                                      length(nu_grid), N_list, N_sim))

## check the effective sample size of all MCMC samples
ESS_MCMC_M <- array(NA, dim = c(N_sim, N_list, 3))
dimnames(ESS_MCMC_M)[[3]] <- c("phi", "nu", "deltasq")


for(i in (1:N_sim)){ 
  filename <- paste0("./sim_carc/results/sim", sim_ind, "_", i, ".Rdata")
  load(filename)
  #print(DIV_matrix)
  
  SPE_stack_LP[i, ] = sqrt(DIV_matrix[, "SPE_stack_LP"])
  SPE_stack_LSE[i, ] = sqrt(DIV_matrix[, "SPE_stack_LSE"])
  SPE_M0[i, ] = sqrt(DIV_matrix[, "SPE_M0"])
  SPE_MCMC[i, ] = sqrt(DIV_matrix[, "SPE_MCMC"])
  ELPD_stack_LSE[i, ] = DIV_matrix[, "ELPD_stack_LSE"]
  ELPD_stack_LP[i, ] = DIV_matrix[, "ELPD_stack_LP"]
  ELPD_M0[i, ] = DIV_matrix[, "ELPD_M0"]
  ELPD_MCMC[i, ] = DIV_matrix[, "ELPD_MCMC"]
  SPE_w_stack_LP[i, ] = sqrt(DIV_matrix[, "SPE_w_stack_LP"])
  SPE_w_stack_LSE[i, ] = sqrt(DIV_matrix[, "SPE_w_stack_LSE"])
  SPE_w_M0[i, ] = sqrt(DIV_matrix[, "SPE_w_M0"])
  SPE_w_MCMC[i, ] = sqrt(DIV_matrix[, "SPE_w_MCMC"])
  
  
  weights_M_LSE_all[, , i] <- weights_M_LSE
  weights_M_LP_all[, , i] <- weights_M_LP
  
  for(j in 1:N_list){
    ESS_MCMC_M[i, j, "phi"] <- effectiveSize(MCMC_par[[j]][, "phi"]) 
    ESS_MCMC_M[i, j, "nu"] <- effectiveSize(MCMC_par[[j]][, "nu"])
    ESS_MCMC_M[i, j, "deltasq"] <- 
      effectiveSize(MCMC_par[[j]][, "tau.sq"] / MCMC_par[[j]][, "sigma.sq"])
  }
}


# Figure 1 & S3: Distributions of the diagnostic metrics for prediction performance 
type = c("M0", "stacking of means", "stacking of predictive densities", 
         "MCMC")

test = c("RMSPE", "RMSEZ", "MLPD")


dat_check <- data.frame(N_sample = rep(rep(rep(paste(samplesize_ls), 
                                               each = N_sim), 
                                       length(type)), length(test)),
                        value = c(c(c(SPE_M0), c(SPE_stack_LSE), 
                                    c(SPE_stack_LP), c(SPE_MCMC)),
                                c(c(SPE_w_M0), c(SPE_w_stack_LSE), 
                                  c(SPE_w_stack_LP), c(SPE_w_MCMC)),
                                c(c(ELPD_M0), c(ELPD_stack_LSE), 
                                 c(ELPD_stack_LP), c(ELPD_MCMC))),
                        label = rep(rep(c(1:4), each = N_list * N_sim), 
                                    length(test)),
                        test = rep(1:3, each = N_list * N_sim * length(type)))

dat_check$label <- factor(dat_check$label, levels = 1:4,
                          labels = type)

dat_check$test <- factor(dat_check$test, levels = 1:3,
                          labels = test)

p_summary <- ggplot(dat_check, aes(x = N_sample, y = value, color = label)) +
  geom_violin(draw_quantiles = c(0.5)) +theme_bw() + xlab("sample size") +
  facet_wrap(~ test, ncol = 1, scales = "free_y", strip.position="right") +
  theme(legend.position="top", legend.title = element_blank()) + ylab(" ") +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"))

p_summary

ggsave(paste0("./sim_carc/pics/CVexperiment_sim", sim_ind, ".png"), 
       plot = p_summary, 
       width = 6.5, height = 4.5, units = "in", dpi = 600)


# Figure 2: Distributions of the counts of nonzero weights 

# On average, only 3.5 out of 64 models have no-zero weights
weights_nonzero_LSE = (weights_M_LSE_all > 0.001)

sum(weights_nonzero_LSE) / (64 * 8) # 3.55 in sim1; 3.87 in sim2; 3.9 in sim3; 3.0 in sim4
weights_nonzero_LP = (weights_M_LP_all > 0.001)
sum(weights_nonzero_LP) / (64 * 8) # 4.13 in sim1; 4.74 in sim2; 4.87 in sim3; 3.65 in sim4

weight_data <- data.frame(
  nonzero_count = c(c(apply(weights_nonzero_LSE, 3:2, sum)), 
                    c(apply(weights_nonzero_LP, 3:2, sum))),
  N_sample = rep(rep(paste(samplesize_ls), each = N_sim), 2),
  label = rep(1:2, each = N_list * N_sim)
)

weight_data$label <- factor(weight_data$label, levels = 1:2,
                          labels = c("stacking of means", 
                                     "stacking of predictive densities"))

p_nonzero_counts <- 
  ggplot(weight_data, aes(x = N_sample, y = nonzero_count, color = label)) +
  geom_violin(draw_quantiles = c(0.5))  + theme_bw() + ylim(c(0, 11)) +
  theme(legend.position = c(0.45, 0.9), legend.title = element_blank(), 
        legend.background = element_blank()) +
  xlab("sample size") + ylab("No. of nonzero weights") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"))
p_nonzero_counts 

## Obviously, stacking based on LP has more nonzero weights than stacking based on LSE
ggsave(paste0("./sim_carc/pics/nonzero_check_sim", sim_ind, ".png"), 
       plot = p_nonzero_counts, 
       width = 6.5, height = 3, units = "in", dpi = 600)

# Historgram of the stacking weights
# Reshape the weight arrays to create a data frame for plotting
# For each sample size and simulation, we have 64 weights

# Create data frame for LSE weights
weights_LSE_long <- data.frame(
  weight_value = c(weights_M_LSE_all),  # Flatten the 64*8*60 array
  N_sample = rep(rep(paste(samplesize_ls), each = 64), N_sim),  # 64 weights per sample size per sim
  sim_rep = rep(1:N_sim, each = 64 * N_list),  # Simulation replicate identifier
  method = "stacking of means"
)

# Create data frame for LP weights  
weights_LP_long <- data.frame(
  weight_value = c(weights_M_LP_all),  # Flatten the 64*8*60 array
  N_sample = rep(rep(paste(samplesize_ls), each = 64), N_sim),  # 64 weights per sample size per sim
  sim_rep = rep(1:N_sim, each = 64 * N_list),  # Simulation replicate identifier
  method = "stacking of predictive densities"
)

# Combine the data frames
weights_all_long <- rbind(weights_LSE_long, weights_LP_long)

# Create violin plots for weight distributions (including zeros with binned small weights)
# Create binned weights for the full dataset: weights < 0.001 are set to 0.0005
weights_all_binned <- weights_all_long
weights_all_binned$weight_value[weights_all_binned$weight_value < 0.001] <- 0.0005

p_weights_dist <- 
  ggplot(weights_all_binned, aes(x = N_sample, y = weight_value, color = method)) +
  geom_violin(alpha = 0.7) + 
  theme_bw() + 
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_blank()) +
  xlab("sample size") + 
  ylab("Stacking weight value (binned small weights)") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  scale_y_log10(limits = c(0.0005, 1), 
                breaks = c(0.0005, 0.001, 0.01, 0.1, 1),
                labels = c("<0.001", "0.001", "0.01", "0.1", "1"))

p_weights_dist

ggsave(paste0("./sim_carc/pics/weights_distribution_sim", sim_ind, ".png"), 
       plot = p_weights_dist, 
       width = 8, height = 5, units = "in", dpi = 600)

# Figure S4: Distributions of stacking weights
# Create violin plots for weight distributions (>0.001 only)
weights_nonzero_LSE <- weights_LSE_long[weights_LSE_long$weight_value > 0.001, ]
weights_nonzero_LP <- weights_LP_long[weights_LP_long$weight_value > 0.001, ]
weights_nonzero_all <- rbind(weights_nonzero_LSE, weights_nonzero_LP)

p_weights_dist_nonzero <- 
  ggplot(weights_nonzero_all, aes(x = N_sample, y = weight_value, color = method)) +
  geom_violin(draw_quantiles = c(0.5), alpha = 0.7) + 
  theme_bw() + 
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_blank()) +
  xlab("sample size") + 
  ylab("Stacking weight value (>0.001 only)") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
  scale_y_continuous(limits = c(0.001, 1), 
                    breaks = c(0.001, 0.25, 0.5, 0.75, 1),
                    labels = c("0.001", "0.25", "0.5", "0.75", "1"))

p_weights_dist_nonzero

ggsave(paste0("./sim_carc/pics/weights_distribution_nonzero_sim", sim_ind, ".png"), 
       plot = p_weights_dist_nonzero, 
       width = 8, height = 5, units = "in", dpi = 600)


## check the inference of hyper-parameters ##
# check phi: #
# Figure S13: Distributions of the estimated $\phi$#
expect_w_phi_LSE <- matrix(0, nrow = length(phi_grid), ncol = N_sim * N_list)
expect_w_phi_LP <- matrix(0, nrow = length(phi_grid), ncol = N_sim * N_list)

for(i in 1:length(phi_grid)){
  expect_w_phi_LSE[i, ] <- c(apply(
    weights_M_LSE_all[which(grid_all$phi == phi_grid[i]), , ], 2:3, sum))
  expect_w_phi_LP[i, ] <- c(apply(
    weights_M_LP_all[which(grid_all$phi == phi_grid[i]), , ], 2:3, sum))
}
est_phi_LSE <- c(t(expect_w_phi_LSE) %*% phi_grid)
est_phi_LP <- c(t(expect_w_phi_LP) %*% phi_grid)
phi_dat <- data.frame(est_phi = c(est_phi_LSE, est_phi_LP), 
                      N_sample = rep(rep(paste(samplesize_ls), N_sim), 2),
                      label = rep(1:2, each = N_list * N_sim))
phi_dat$label <- factor(phi_dat$label, levels = 1:2,
                        labels = c("stacking of means", 
                                   "stacking of predictive densities"))
# adjust the position of the legend for different simulation
if((sim_ind == 1)| (sim_ind == 4)){
  leg_pos = c(0.8, 0.85)
}else if (sim_ind == 2){
  leg_pos = c(0.8, 0.2)
}else if (sim_ind == 3){
  leg_pos = c(0.8, 0.88)
}
p_est_phi <- 
  ggplot(phi_dat, aes(x = N_sample, y = est_phi, color = label)) +
  geom_violin(draw_quantiles = c(0.5))  + theme_bw() + ylim(c(0, 24)) +
  theme(legend.position = leg_pos, legend.title = element_blank(), 
        legend.background = element_blank()) +
  xlab("sample size") + ylab(expression("estimated "*phi)) + 
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) + 
  geom_hline(yintercept = raw_data[[1]]$phi, linetype="dashed", color = "red")
p_est_phi
# the estimation of phi is not reliable
ggsave(paste0("./sim_carc/pics/est_phi_sim", sim_ind, ".png"), 
       plot = p_est_phi, 
       width = 6.5, height = 3, units = "in", dpi = 600)

# check nu #
# Figure S14: Distributions of the estimated $\nu$#

expect_w_nu_LSE <- matrix(0, nrow = length(nu_grid), ncol = N_sim * N_list)
expect_w_nu_LP <- matrix(0, nrow = length(nu_grid), ncol = N_sim * N_list)

for(i in 1:length(nu_grid)){
  expect_w_nu_LSE[i, ] <- c(apply(
    weights_M_LSE_all[which(grid_all$nu == nu_grid[i]), , ], 2:3, sum))
  expect_w_nu_LP[i, ] <- c(apply(
    weights_M_LP_all[which(grid_all$nu == nu_grid[i]), , ], 2:3, sum))
}
est_nu_LSE <- c(t(expect_w_nu_LSE) %*% nu_grid)
est_nu_LP <- c(t(expect_w_nu_LP) %*% nu_grid)
nu_dat <- data.frame(est_nu = c(est_nu_LSE, est_nu_LP), 
                      N_sample = rep(rep(paste(samplesize_ls), N_sim), 2),
                      label = rep(1:2, each = N_list * N_sim))
nu_dat$label <- factor(nu_dat$label, levels = 1:2,
                        labels = c("stacking of means", 
                                   "stacking of predictive densities"))
# adjust the position of the legend for different simulation
if((sim_ind == 1) | (sim_ind == 4) ){
  leg_pos_nu = c(0.8, 0.15)
}else if (sim_ind == 2){
  leg_pos_nu = c(0.8, 0.91)
}else if (sim_ind == 3){
  leg_pos_nu = c(0.8, 0.12)
}

p_est_nu <- 
  ggplot(nu_dat, aes(x = N_sample, y = est_nu, color = label)) +
  geom_violin(draw_quantiles = c(0.5))  + theme_bw() + ylim(c(0.2, 2)) +
  theme(legend.position = leg_pos_nu, legend.title = element_blank(), 
        legend.background = element_blank()) +
  xlab("sample size") + ylab(expression("estimated "*nu)) + 
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) + 
  geom_hline(yintercept = raw_data[[1]]$nu, linetype="dashed", color = "red")
p_est_nu
# the estimation of phi is not reliable
ggsave(paste0("./sim_carc/pics/est_nu_sim", sim_ind, ".png"), 
       plot = p_est_nu, 
       width = 6.5, height = 3, units = "in", dpi = 600)

# check deltasq #
# Figure S15: Distributions of the estimated $\delta^2$#
expect_w_deltasq_LSE <- matrix(0, nrow = length(deltasq_grid), ncol = N_sim * N_list)
expect_w_deltasq_LP <- matrix(0, nrow = length(deltasq_grid), ncol = N_sim * N_list)

for(i in 1:length(deltasq_grid)){
  expect_w_deltasq_LSE[i, ] <- c(apply(
    weights_M_LSE_all[which(grid_all$deltasq == deltasq_grid[i]), , ], 
    2:3, sum))
  expect_w_deltasq_LP[i, ] <- c(apply(
    weights_M_LP_all[which(grid_all$deltasq == deltasq_grid[i]), , ], 
    2:3, sum))
}
est_deltasq_LSE <- c(t(expect_w_deltasq_LSE) %*% deltasq_grid)
est_deltasq_LP <- c(t(expect_w_deltasq_LP) %*% deltasq_grid)
deltasq_dat <- data.frame(est_deltasq = c(est_deltasq_LSE, est_deltasq_LP), 
                     N_sample = rep(rep(paste(samplesize_ls), N_sim), 2),
                     label = rep(1:2, each = N_list * N_sim))
deltasq_dat$label <- factor(deltasq_dat$label, levels = 1:2,
                       labels = c("stacking of means", 
                                  "stacking of predictive densities"))
# adjust the position of the legend for different simulation
if( (sim_ind == 1) | (sim_ind == 2) ){
  leg_pos_deltasq = "none"
}else if (sim_ind == 3){
  leg_pos_deltasq = c(0.8, 0.88)
}else if (sim_ind == 4){
  leg_pos_deltasq = c(0.8, 0.1) #c(0.8, 0.88)
}
p_est_deltasq <- 
  ggplot(deltasq_dat, aes(x = N_sample, y = est_deltasq, color = label)) +
  geom_violin(draw_quantiles = c(0.5))  + theme_bw() + ylim(c(0, 2.1)) +
  theme(legend.position = leg_pos_deltasq, legend.title = element_blank(), 
        legend.background = element_blank()) +
  xlab("sample size") + ylab(expression("estimated "*delta^2)) + 
  scale_colour_manual(values=c("#E69F00", "#56B4E9")) + 
  geom_hline(yintercept = (raw_data[[1]]$tau.sq / raw_data[[1]]$sigma.sq), 
                           linetype="dashed", color = "red")
p_est_deltasq
# the estimation of phi is not reliable
ggsave(paste0("./sim_carc/pics/est_deltasq_sim", sim_ind, ".png"), 
       plot = p_est_deltasq, 
       width = 6.5, height = 3, units = "in", dpi = 600)


## Check the effective sample size of MCMC ##
ESS_MCMC_dat <- 
  data.frame(ESS = c(ESS_MCMC_M), 
             N_sample = rep(rep(paste(samplesize_ls), each = N_sim), 3),
             label = rep(1:3, each = N_list * N_sim))
ESS_MCMC_dat$label <- factor(ESS_MCMC_dat$label, levels = 1:3,
                             labels = dimnames(ESS_MCMC_M)[[3]])

if((sim_ind == 1) |(sim_ind == 3) ){
  leg_pos_ESS = "none"
}else if (sim_ind == 2){
  leg_pos_ESS = c(0.8, 0.9)
}else if (sim_ind == 4){
  leg_pos_ESS = c(0.8, 0.9)
}
p_ESS_MCMC<- 
  ggplot(ESS_MCMC_dat, aes(x = N_sample, y = ESS, color = label)) +
  geom_violin(draw_quantiles = c(0.5))  + theme_bw() + #ylim(c(0, 900)) +
  theme(legend.position = leg_pos_ESS, legend.title = element_blank(), 
        legend.background = element_blank(),
        legend.direction="horizontal") + 
  scale_y_log10(breaks = c(20, 50, seq(0, 600, by = 100), 800, 1000)) +
  xlab("sample size") + ylab("MCMC ESS") + 
  scale_colour_manual(values=c("#009E73", "#F0E442", "#0072B2")) 
p_ESS_MCMC
ggsave(paste0("./sim_carc/pics/ESS_MCMC_sim", sim_ind, ".png"), 
       plot = p_ESS_MCMC, 
       width = 6.5, height = 3, units = "in", dpi = 600)

