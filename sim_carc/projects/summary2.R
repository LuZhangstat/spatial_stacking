rm(list = ls())
library(ggplot2)
library(MASS)
library(spBayes)
library(ggplot2)
library("gridExtra")
library("coda")

#' Figure  S12: Distributions of the diagnostic metrics for prediction performance 

# colorblind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim_ind = 4 # simulation index 1, 2, 3 or 4


load(paste0("./sim_carc/results/sim", sim_ind, "_1.Rdata"))
K_fold = 10
N_list = length(samplesize_ls)
N_sim = 60
SPE_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_stack_LP_P = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_stack_LSE_P = matrix(NA, nrow = N_sim, ncol = N_list)

ELPD_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_stack_LSE_P = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_stack_LP_P = matrix(NA, nrow = N_sim, ncol = N_list)

SPE_w_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_stack_LP_P = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_stack_LSE_P = matrix(NA, nrow = N_sim, ncol = N_list)

for(i in 1:N_sim){   
  filename <- paste0("./sim_carc/results/sim", sim_ind, "_", i, ".Rdata")
  load(filename)
  SPE_MCMC[i, ] = sqrt(DIV_matrix[, "SPE_MCMC"])
  SPE_stack_LP[i, ] = sqrt(DIV_matrix[, "SPE_stack_LP"])
  SPE_stack_LSE[i, ] = sqrt(DIV_matrix[, "SPE_stack_LSE"])
  SPE_stack_LP_P[i, ] = sqrt(DIV_matrix[, "SPE_stack_LP_P"])
  SPE_stack_LSE_P[i, ] = sqrt(DIV_matrix[, "SPE_stack_LSE_P"])
  
  ELPD_MCMC[i, ] = DIV_matrix[, "ELPD_MCMC"]
  ELPD_stack_LSE[i, ] = DIV_matrix[, "ELPD_stack_LSE"]
  ELPD_stack_LP[i, ] = DIV_matrix[, "ELPD_stack_LP"]
  ELPD_stack_LSE_P[i, ] = DIV_matrix[, "ELPD_stack_LSE_P"]
  ELPD_stack_LP_P[i, ] = DIV_matrix[, "ELPD_stack_LP_P"]
  
  SPE_w_MCMC[i, ] = sqrt(DIV_matrix[, "SPE_w_MCMC"])
  SPE_w_stack_LP[i, ] = sqrt(DIV_matrix[, "SPE_w_stack_LP"])
  SPE_w_stack_LSE[i, ] = sqrt(DIV_matrix[, "SPE_w_stack_LSE"])
  SPE_w_stack_LP_P[i, ] = sqrt(DIV_matrix[, "SPE_w_stack_LP_P"])
  SPE_w_stack_LSE_P[i, ] = sqrt(DIV_matrix[, "SPE_w_stack_LSE_P"])
}


type = c("MCMC LSE", "stacking of means", "MCMC LP", 
         "stacking of predictive densities", "MCMC")

test = c("RMSPE", "RMSEZ", "MLPD")


dat_check <- data.frame(N_sample = rep(rep(rep(paste(samplesize_ls), 
                                               each = N_sim), 
                                       length(type)), length(test)),
                        value = c(c(c(SPE_stack_LSE_P), c(SPE_stack_LSE), 
                                    c(SPE_stack_LP_P), c(SPE_stack_LP),
                                    c(SPE_MCMC)),
                                c(c(SPE_w_stack_LSE_P), c(SPE_w_stack_LSE), 
                                  c(SPE_w_stack_LP_P), c(SPE_w_stack_LP), 
                                  c(SPE_w_MCMC)),
                                c(c(ELPD_stack_LSE_P), c(ELPD_stack_LSE), 
                                  c(ELPD_stack_LP_P), c(ELPD_stack_LP),
                                  c(ELPD_MCMC))),
                        label = rep(rep(c(1:length(type)), 
                                        each = N_list * N_sim), 
                                    length(test)),
                        test = rep(1:length(test), 
                                   each = N_list * N_sim * length(type)))

dat_check$label <- factor(dat_check$label, levels = 1:length(type),
                          labels = type)

dat_check$test <- factor(dat_check$test, levels = 1:length(test),
                          labels = test)

p_summary <- ggplot(dat_check, aes(x = N_sample, y = value, color = label)) +
  geom_violin(draw_quantiles = c(0.5)) +theme_bw() + xlab("sample size") +
  facet_wrap(~ test, ncol = 1, scales = "free_y", strip.position="right") +
  theme(legend.position="top", legend.title = element_blank()) + ylab(" ") +
  scale_colour_manual(values=c("#D55E00", "#E69F00", "#F0E442", #"#999999", 
                               "#0072B2", "#56B4E9"#, "#009E73"#, "#CC79A7"
                               ))
p_summary

ggsave(paste0("./sim_carc/pics/CVexperiment_prefix_sim", sim_ind, ".png"), 
       plot = p_summary, 
       width = 6.5, height = 4.5, units = "in", dpi = 600)

