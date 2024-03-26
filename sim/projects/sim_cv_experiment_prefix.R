rm(list = ls())
library(ggplot2)
library(MASS)
library(spBayes)
library(ggplot2)
library(cmdstanr)
library(geoR)
library(rbenchmark)
library("gridExtra")
library(fields)
library(INLA)
source("utils.R") # utils2.R is the testing code #
#options(mc.cores = parallel::detectCores())

############################################################
######## the impact of candidate value for sigma^2 #########
############################################################
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

load("./sim_carc/results/sim1_1.RData")
# default method #
sim_ind = 1
phi_grid = c(3, 14, 25, 36)   #3.5 to 35.8 #old: c(3, 9, 15, 31) 
nu_grid = c(0.5, 1, 1.5, 1.75)
deltasq_grid <- pick_deltasq(E_sigmasq = raw_data[[1]]$sigma.sq, 
                             E_tausq = raw_data[[1]]$tau.sq, 
                             b = max(raw_data[[1]]$tau.sq, raw_data[[1]]$sigma.sq),
                             p_ls = c(0.05, 0.35, 0.65, 0.95))
deltasq_grid
seed = 123

r = 8 # r = 2,8 for sim1 r = 4 for sim2
ind_mod = raw_data[[r]]$ind_mod
X <- raw_data[[r]]$X
y <- raw_data[[r]]$y
w <- raw_data[[r]]$w
coords <- raw_data[[r]]$coords
N <- samplesize_ls[r]
N_ho <- N - length(raw_data[[r]]$ind_mod)

priors <- list(mu_beta = rep(0, 2),
               inv_V_beta = 1/4 * diag(2),
               a_sigma = 2,
               b_sigma = 2)
K_fold = 10 # K-fold cross-validation


CV_fit_LSE <- sp_stacking_K_fold(
  X = X[ind_mod, ], y = y[ind_mod], 
  coords = coords[ind_mod, ],
  deltasq_grid = deltasq_grid, phi_grid = phi_grid,
  nu_grid = nu_grid, priors = priors, 
  K_fold = K_fold, seed = seed, label = "LSE")
cbind(CV_fit_LSE$grid_all[CV_fit_LSE$wts>0.00001, ], 
      CV_fit_LSE$wts[CV_fit_LSE$wts>0.00001])

pos_sam_LSE <- 
  stacking_pos_sample(Stack_fit = CV_fit_LSE, L1 = 900, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 1,
                      recover_w_obs = T, recover_w_U = T)


CV_fit_LP <- sp_stacking_K_fold(
  X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
  deltasq_grid = deltasq_grid, phi_grid = phi_grid, 
  nu_grid = nu_grid,
  priors = priors, K_fold = K_fold,
  seed = seed, label = "LP", MC = FALSE)
cbind(CV_fit_LP$grid_all[CV_fit_LP$wts>0.00001, ], 
      CV_fit_LP$wts[CV_fit_LP$wts>0.00001])

pos_sam_LP <- 
  stacking_pos_sample(Stack_fit = CV_fit_LP, L1 = 900, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 2,
                      recover_w_obs = T, recover_w_U = T)

cat("95%CI coverage:")
y_U_CI_stack_LP <- 
  apply(pos_sam_LP$pred_y_U_stack_sam, 1, 
        function(x){quantile(x, probs = c(0.025, 0.975))})
sum((y_U_CI_stack_LP[1, ] < y[-ind_mod]) & 
      (y_U_CI_stack_LP[2, ] > y[-ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
# 0.96 for sim1 r = 8; 0.88 for sim2 r = 4 

cat("98%CI coverage:")
y_U_CI_stack_LP <- 
  apply(pos_sam_LP$pred_y_U_stack_sam, 1, 
        function(x){quantile(x, probs = c(0.01, 0.99))})
sum((y_U_CI_stack_LP[1, ] < y[-ind_mod]) & 
      (y_U_CI_stack_LP[2, ] > y[-ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
#0.98 for sim1 r = 8; 0.93 for sim 2 r = 4

# for w #
cat("obs w 95%CI coverage:")
w_obs_CI_stack_LP <- 
  apply(pos_sam_LP$pred_w_obs_stack_sam, 1, 
        function(x){quantile(x, probs = c(0.025, 0.975))})
sum((w_obs_CI_stack_LP[1, ] < w[ind_mod]) & 
      (w_obs_CI_stack_LP[2, ] > w[ind_mod]))/length(w[raw_data[[r]]$ind_mod])
# 1 for sim1 r=8; 0.935 for sim2 r = 4 

cat("unobs w 95%CI coverage:")
w_U_CI_stack_LP <- 
  apply(pos_sam_LP$pred_w_U_stack_sam, 1, 
        function(x){quantile(x, probs = c(0.025, 0.975))})
sum((w_U_CI_stack_LP[1, ] < w[-ind_mod]) & 
      (w_U_CI_stack_LP[2, ] > w[-ind_mod]))/length(w[-raw_data[[r]]$ind_mod])
# 1 for sim1 r = 8; 0.68 for sim2 r = 4 

# INLA #
df = data.frame(y=c(raw_data[[r]]$y[raw_data[[r]]$ind_mod], rep(NA, 100)), 
                locx=raw_data[[r]]$coords[, 1], 
                locy=raw_data[[r]]$coords[, 2], x = raw_data[[r]]$X[, 2])
fake.locations = matrix(c(0,0,1,1, 0, 1, 1, 0), nrow = 4, byrow = T)
n <- length(raw_data[[r]]$coords[, 1]) - 100
mesh = inla.mesh.2d(loc = fake.locations, max.edge=c(2/sqrt(n), 10/sqrt(n)))
mesh$n

A = inla.spde.make.A(mesh=mesh, loc=data.matrix(df[ , c('locx', 'locy')]))
prior.median.sd = 1; prior.median.range = 1/raw_data[[r]]$phi

a1 <- raw_data[[r]]$nu +1
spde1 = inla.spde2.pcmatern(mesh, alpha = a1,
                            prior.range = c(prior.median.range, .5), 
                            prior.sigma = c(prior.median.sd, .5))
stack1 = inla.stack(tag='est',
                    # - Name (nametag) of the stack
                    # - Here: est for estimating
                    data=list(y=df$y),
                    effects=list(
                      # - The Model Components
                      s=1:spde1$n.spde, 
                      # - The first is 's' (for spatial)
                      data.frame(intercept=1, x=df$x)),
                    # - The second is all fixed effects
                    A=list(A, 1)
)

family = "gaussian"
prior.median.sd.g = 1 # prior median for sigma.epsilon
control.family = list(hyper = list(prec = list(
  prior = "pc.prec", param = c(prior.median.sd.g,0.5))))

initial.theta = NULL #c(2.35, 0.79, 0.46)

res1 = inla(y~x + f(s, model=spde1), data=inla.stack.data(stack1),
            family = family,
            control.family = control.family,
            control.predictor=list(A = inla.stack.A(stack1)),
            quantiles=c(0.01, 0.025, 0.1, 0.5, 0.9, 0.975, 0.99),
            control.mode = list(restart = T, theta = initial.theta))

cat("95%CI coverage:")
sum((res1$summary.fitted.values[
  inla.stack.index(stack1, "est")$data[-raw_data[[r]]$ind_mod],
  "0.025quant"] < y[-raw_data[[r]]$ind_mod]) &
    (res1$summary.fitted.values[
      inla.stack.index(stack1, "est")$data[-raw_data[[r]]$ind_mod],
      "0.975quant"] > y[-raw_data[[r]]$ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
# 0.52 for sim1 r = 8; 0.55 for sim2 with r = 4

cat("98%CI coverage:")
sum((res1$summary.fitted.values[
  inla.stack.index(stack1, "est")$data[-raw_data[[r]]$ind_mod],
  "0.01quant"] < y[-raw_data[[r]]$ind_mod]) &
    (res1$summary.fitted.values[
      inla.stack.index(stack1, "est")$data[-raw_data[[r]]$ind_mod],
      "0.99quant"] > y[-raw_data[[r]]$ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
# 0.58 for sim1 r = 8; 0.62 for sim2 r = 4

# test 2: posterior sample from marginal distribution #
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
cbind(CV_fit_LSE_P$grid_all[CV_fit_LSE_P$wts>0.00001, ], 
      CV_fit_LSE_P$wts[CV_fit_LSE_P$wts>0.00001])

pos_sam_LSE_P <- 
  stacking_pos_sample(Stack_fit = CV_fit_LSE_P, L1 = 900, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 5,
                      recover_w_obs = T, recover_w_U = T)


CV_fit_LP_P <- sp_stacking_K_fold3(
  X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
  all_prefix_ls = all_prefix_ls, 
  priors = priors, K_fold = K_fold,
  seed = seed, label = "LP", MC = FALSE)
cbind(CV_fit_LP_P$grid_all[CV_fit_LP_P$wts>0.00001, ], 
      CV_fit_LP_P$wts[CV_fit_LP_P$wts>0.00001])

pos_sam_LP_P <- 
  stacking_pos_sample(Stack_fit = CV_fit_LP_P, L1 = 900, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 6,
                      recover_w_obs = T, recover_w_U = T)

cat("95%CI coverage:")
y_U_CI_stack_LP_P <- 
  apply(pos_sam_LP_P$pred_y_U_stack_sam, 1, 
        function(x){quantile(x, probs = c(0.025, 0.975))})
sum((y_U_CI_stack_LP_P[1, ] < y[-ind_mod]) & 
      (y_U_CI_stack_LP_P[2, ] > y[-ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
#0.96 for sim1 r=8; 0.88 for sim2 r = 4 

cat("98%CI coverage:")
y_U_CI_stack_LP_P <- 
  apply(pos_sam_LP_P$pred_y_U_stack_sam, 1, 
        function(x){quantile(x, probs = c(0.01, 0.99))})
sum((y_U_CI_stack_LP_P[1, ] < y[-ind_mod]) & 
      (y_U_CI_stack_LP_P[2, ] > y[-ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
#0.98 for sim1 r=8; 0.93 for sim 2 r = 4


# compare sigmasq #
draws_ls_sigmasq <- c()
draws_ls_sigmasq[[1]] <- MCMC_par[[r]][101:1000, "sigma.sq"]
draws_ls_sigmasq[[2]] <- pos_sam_LSE$sigmasq_sam
draws_ls_sigmasq[[3]] <- pos_sam_LSE_P$sigmasq_sam
draws_ls_sigmasq[[4]] <- MCMC_par[[r]][101:1000, "sigma.sq"]
draws_ls_sigmasq[[5]] <- pos_sam_LP$sigmasq_sam
draws_ls_sigmasq[[6]] <- pos_sam_LP_P$sigmasq_sam
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sigmasq_compar <- 
  hist_compar(draws_ls = draws_ls_sigmasq, 
              type_colors = c("#999999", "#69b3a2", "#404080"),  #"#56B4E9",
              type_names = c("MCMC", "default",  "MCMC+stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = raw_data[[r]]$sigma.sq,
              yname = expression(sigma^2), bins = 45)

sigmasq_compar
ggsave(paste0("./sim/pics/sigmasq_prefix_compar", sim_ind, "_r", r, ".png"),
       plot = sigmasq_compar,
       width = 6.5, height = 3.5, units = "in", dpi = 600)


# compare tausq #
draws_ls_tausq <- c()
draws_ls_tausq[[1]] <- MCMC_par[[r]][101:1000, "tau.sq"]
draws_ls_tausq[[2]] <- pos_sam_LSE$tausq_sam
draws_ls_tausq[[3]] <- pos_sam_LSE_P$tausq_sam
draws_ls_tausq[[4]] <- MCMC_par[[r]][101:1000, "tau.sq"]
draws_ls_tausq[[5]] <- pos_sam_LP$tausq_sam
draws_ls_tausq[[6]] <- pos_sam_LP_P$tausq_sam
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
tausq_compar <- 
  hist_compar(draws_ls = draws_ls_tausq, 
              type_colors = c("#999999", "#69b3a2", "#404080"),  #"#56B4E9",
              type_names = c("MCMC", "default",  "MCMC+stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = raw_data[[r]]$tau.sq,
              yname = expression(tau^2), bins = 45)

tausq_compar
ggsave(paste0("./sim/pics/tausq_prefix_compar", sim_ind, "_r", r, ".png"),
       plot = tausq_compar,
       width = 6.5, height = 3.5, units = "in", dpi = 600)



# check histograms for individuals (pick the 50th and 100th point)#
### fit with spBayes ###
fit_flag <- FALSE
if(fit_flag){
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
  
  ## recover posterior samples for y and w ##
  pos_wy <- recover_MCMC(theta.recover = r.1$p.theta.recover.samples,
                         beta.recover = r.1$p.beta.recover.samples,
                         y.mod = y[ind_mod], X.mod = X[ind_mod, ],
                         coords.mod = coords[ind_mod, ],
                         X.ho = X[-ind_mod, ], y.ho = y[-ind_mod],
                         coords.ho = coords[-ind_mod, ])
  save(m.1, r.1, pos_wy, 
       file = paste0("./sim/results/indi_compar", sim_ind, "_r", r, ".Rdata"))
}else{
  load(paste0("./sim/results/indi_compar", sim_ind, "_r", r, ".Rdata"))
}

cat("95%CI coverage:")
pos_y_U_CI_P <- 
  apply(pos_wy$y.ho.sample, 1, 
        function(x){quantile(x, probs = c(0.025, 0.975))})
sum((pos_y_U_CI_P[1, ] < y[-ind_mod]) & 
      (pos_y_U_CI_P[2, ] > y[-ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
# 0.96 for sim1 r=8;  0.89 for sim2 r = 4 

cat("98%CI coverage:")
pos_y_U_CI_P <- 
  apply(pos_wy$y.ho.sample, 1, 
        function(x){quantile(x, probs = c(0.01, 0.99))})
sum((pos_y_U_CI_P[1, ] < y[-ind_mod]) & 
      (pos_y_U_CI_P[2, ] > y[-ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
#0.99 for sim1 r=8; 0.94 for sim 2 r = 4


# compare beta1 #
draws_ls_b1 <- c()
draws_ls_b1[[1]] <- r.1$p.beta.recover.samples[101:1000, 1]
draws_ls_b1[[2]] <- pos_sam_LSE$pred_beta_stack_sam[1, ]
draws_ls_b1[[3]] <- pos_sam_LSE_P$pred_beta_stack_sam[1, ]
draws_ls_b1[[4]] <- r.1$p.beta.recover.samples[101:1000, 1]
draws_ls_b1[[5]] <- pos_sam_LP$pred_beta_stack_sam[1, ]
draws_ls_b1[[6]] <- pos_sam_LP_P$pred_beta_stack_sam[1, ]
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
b1_compar <- 
  hist_compar(draws_ls = draws_ls_b1, 
              type_colors = c("#999999", "#69b3a2", "#404080"),  #"#56B4E9",
              type_names = c("MCMC", "default", "MCMC+stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = raw_data[[r]]$beta[1],
              yname = expression(beta[1]), bins = 45)

b1_compar
ggsave(paste0("./sim/pics/beta1_prefix_compar", sim_ind, "_r", r, ".png"),
       plot = b1_compar,
       width = 6.5, height = 3.5, units = "in", dpi = 600)


# compare beta2 #
draws_ls_b2 <- c()
draws_ls_b2[[1]] <- r.1$p.beta.recover.samples[101:1000, 2]
draws_ls_b2[[2]] <- pos_sam_LSE$pred_beta_stack_sam[2, ]
draws_ls_b2[[3]] <- pos_sam_LSE_P$pred_beta_stack_sam[2, ]
draws_ls_b2[[4]] <- r.1$p.beta.recover.samples[101:1000, 2]
draws_ls_b2[[5]] <- pos_sam_LP$pred_beta_stack_sam[2, ]
draws_ls_b2[[6]] <- pos_sam_LP_P$pred_beta_stack_sam[2, ]
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
b2_compar <- 
  hist_compar(draws_ls = draws_ls_b2, 
              type_colors = c("#999999", "#69b3a2", "#404080"),  #"#56B4E9",
              type_names = c("MCMC", "default", "MCMC+stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = raw_data[[r]]$beta[2],
              yname = expression(beta[2]), bins = 45)

b2_compar
ggsave(paste0("./sim/pics/beta2_prefix_compar", sim_ind, "_r", r, ".png"),
       plot = b2_compar,
       width = 6.5, height = 3.5, units = "in", dpi = 600)


# check y #
pick_indi <- c(50, 90)
draws_ls1 <- c()
draws_ls1[[1]] <- pos_wy$y.ho.sample[pick_indi[1], 101:1000]
draws_ls1[[2]] <- pos_sam_LSE$pred_y_U_stack_sam[pick_indi[1], ]
draws_ls1[[3]] <- pos_sam_LSE_P$pred_y_U_stack_sam[pick_indi[1], ]
draws_ls1[[4]] <- pos_wy$y.ho.sample[pick_indi[1], 101:1000]
draws_ls1[[5]] <- pos_sam_LP$pred_y_U_stack_sam[pick_indi[1], ]
draws_ls1[[6]] <- pos_sam_LP_P$pred_y_U_stack_sam[pick_indi[1], ]

individual_y_compar1 <- 
  hist_compar(draws_ls = draws_ls1, 
              type_colors = c("#999999", "#69b3a2", "#404080"), 
              type_names = c("MCMC", "default", "MCMC+Stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = y[-ind_mod][pick_indi[1]],
              yname = expression(y(s[50]))
  )

draws_ls2 <- c()
draws_ls2[[1]] <- pos_wy$y.ho.sample[pick_indi[2], 101:1000]
draws_ls2[[2]] <- pos_sam_LSE$pred_y_U_stack_sam[pick_indi[2], ]
draws_ls2[[3]] <- pos_sam_LSE_P$pred_y_U_stack_sam[pick_indi[2], ]
draws_ls2[[4]] <- pos_wy$y.ho.sample[pick_indi[2], 101:1000]
draws_ls2[[5]] <- pos_sam_LP$pred_y_U_stack_sam[pick_indi[2], ]
draws_ls2[[6]] <- pos_sam_LP_P$pred_y_U_stack_sam[pick_indi[2], ]

individual_y_compar2 <- 
  hist_compar(draws_ls = draws_ls2, 
              type_colors = c("#999999", "#69b3a2", "#404080"), 
              type_names = c("MCMC", "default", "MCMC+stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = y[-ind_mod][pick_indi[2]],
              yname = expression(y(s[90]))
  )

print(individual_y_compar1)
print(individual_y_compar2)
ggsave(paste0("./sim/pics/indi50_prefix_compar", sim_ind, "ICI_r", r, ".png"),
       plot = individual_y_compar1,
       width = 6.5, height = 3.5, units = "in", dpi = 600)
ggsave(paste0("./sim/pics/indi90_prefix_compar", sim_ind, "ICI_r", r, ".png"),
       plot = individual_y_compar2,
       width = 6.5, height = 3.5, units = "in", dpi = 600)

# check latent process z #
pick_indi <- c(50, 90)
draws_ls1 <- c()
draws_ls1[[1]] <- pos_wy$w.recover.sample[r*100+pick_indi[1], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls1[[2]] <- pos_sam_LSE$pred_w_U_stack_sam[pick_indi[1], ] +
  pos_sam_LSE$pred_beta_stack_sam[1, ]
draws_ls1[[3]] <- pos_sam_LSE_P$pred_w_U_stack_sam[pick_indi[1], ]+
  pos_sam_LSE_P$pred_beta_stack_sam[1, ]
draws_ls1[[4]] <- pos_wy$w.recover.sample[r*100+pick_indi[1], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls1[[5]] <- pos_sam_LP$pred_w_U_stack_sam[pick_indi[1], ] +
  pos_sam_LP$pred_beta_stack_sam[1, ]
draws_ls1[[6]] <- pos_sam_LP_P$pred_w_U_stack_sam[pick_indi[1], ] +
  pos_sam_LP_P$pred_beta_stack_sam[1, ]

individual_w_compar1 <- 
  hist_compar(draws_ls = draws_ls1, 
              type_colors = c("#999999", "#69b3a2", "#404080"), 
              type_names = c("MCMC", "default", "MCMC+Stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = w[-ind_mod][pick_indi[1]]+raw_data[[r]]$beta[1],
              yname = expression(z(s[50])+beta[1])
  )

draws_ls2 <- c()
draws_ls2[[1]] <- pos_wy$w.recover.sample[r*100+pick_indi[2], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls2[[2]] <- pos_sam_LSE$pred_w_U_stack_sam[pick_indi[2], ] +
  pos_sam_LSE$pred_beta_stack_sam[1, ]
draws_ls2[[3]] <- pos_sam_LSE_P$pred_w_U_stack_sam[pick_indi[2], ]+
  pos_sam_LSE_P$pred_beta_stack_sam[1, ]
draws_ls2[[4]] <- pos_wy$w.recover.sample[r*100+pick_indi[2], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls2[[5]] <- pos_sam_LP$pred_w_U_stack_sam[pick_indi[2], ] +
  pos_sam_LP$pred_beta_stack_sam[1, ]
draws_ls2[[6]] <- pos_sam_LP_P$pred_w_U_stack_sam[pick_indi[2], ] +
  pos_sam_LP_P$pred_beta_stack_sam[1, ]

individual_w_compar2 <- 
  hist_compar(draws_ls = draws_ls2, 
              type_colors = c("#999999", "#69b3a2", "#404080"), 
              type_names = c("posterior", "default", "MCMC+stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = w[-ind_mod][pick_indi[2]]+raw_data[[r]]$beta[1],
              yname = expression(z(s[90])+beta[1])
  )

print(individual_w_compar1)
print(individual_w_compar2)
ggsave(paste0("./sim/pics/indi50_prefix_incpw_compar", sim_ind, "_r", r, ".png"),
       plot = individual_w_compar1,
       width = 6.5, height = 3.5, units = "in", dpi = 600)
ggsave(paste0("./sim/pics/indi90_prefix_incpw_compar", sim_ind, "_r", r, ".png"),
       plot = individual_w_compar2,
       width = 6.5, height = 3.5, units = "in", dpi = 600)

## check the latent process at observed locations ##
# check w #
pick_indi <- c(50, 100)
draws_ls1 <- c()
draws_ls1[[1]] <- pos_wy$w.recover.sample[pick_indi[1], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls1[[2]] <- pos_sam_LSE$pred_w_obs_stack_sam[pick_indi[1], ] +
  pos_sam_LSE$pred_beta_stack_sam[1, ]
draws_ls1[[3]] <- pos_sam_LSE_P$pred_w_obs_stack_sam[pick_indi[1], ]+
  pos_sam_LSE_P$pred_beta_stack_sam[1, ]
draws_ls1[[4]] <- pos_wy$w.recover.sample[pick_indi[1], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls1[[5]] <- pos_sam_LP$pred_w_obs_stack_sam[pick_indi[1], ] +
  pos_sam_LP$pred_beta_stack_sam[1, ]
draws_ls1[[6]] <- pos_sam_LP_P$pred_w_obs_stack_sam[pick_indi[1], ] +
  pos_sam_LP_P$pred_beta_stack_sam[1, ]

individual_w_compar1 <- 
  hist_compar(draws_ls = draws_ls1, 
              type_colors = c("#999999", "#69b3a2", "#404080"), 
              type_names = c("MCMC", "default", "MCMC+Stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = w[pick_indi[1]]+raw_data[[r]]$beta[1],
              yname = expression(z(s[50])+beta[1])
  )

draws_ls2 <- c()
draws_ls2[[1]] <- pos_wy$w.recover.sample[pick_indi[2], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls2[[2]] <- pos_sam_LSE$pred_w_obs_stack_sam[pick_indi[2], ] +
  pos_sam_LSE$pred_beta_stack_sam[1, ]
draws_ls2[[3]] <- pos_sam_LSE_P$pred_w_obs_stack_sam[pick_indi[2], ]+
  pos_sam_LSE_P$pred_beta_stack_sam[1, ]
draws_ls2[[4]] <- pos_wy$w.recover.sample[pick_indi[2], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls2[[5]] <- pos_sam_LP$pred_w_obs_stack_sam[pick_indi[2], ] +
  pos_sam_LP$pred_beta_stack_sam[1, ]
draws_ls2[[6]] <- pos_sam_LP_P$pred_w_obs_stack_sam[pick_indi[2], ] +
  pos_sam_LP_P$pred_beta_stack_sam[1, ]

individual_w_compar2 <- 
  hist_compar(draws_ls = draws_ls2, 
              type_colors = c("#999999", "#69b3a2", "#404080"), 
              type_names = c("posterior", "default", "MCMC+stacking"), 
              test_names = c("stacking of means", "stacking of pds"), 
              true_value = w[pick_indi[2]]+raw_data[[r]]$beta[1],
              yname = expression(z(s[90])+beta[1])
  )

print(individual_w_compar1)
print(individual_w_compar2)
ggsave(paste0("./sim/pics/indi50_prefix_incpw_obs_compar", sim_ind, "_r", r, ".png"),
       plot = individual_w_compar1,
       width = 6.5, height = 3.5, units = "in", dpi = 600)
ggsave(paste0("./sim/pics/indi90_prefix_incpw_obs_compar", sim_ind, "_r", r, ".png"),
       plot = individual_w_compar2,
       width = 6.5, height = 3.5, units = "in", dpi = 600)

# check all y and w #
library(ggplot2)
library(dplyr)
# Assuming y, ind_mod, pos_sam_LP$pred_y_U_stack_sam, and pos_sam_LSE$pred_y_U_stack_sam are defined
# Function to calculate means and CIs
calculate_means_and_cis <- function(pos_sam) {
  row_means <- rowMeans(pos_sam)
  CI <- apply(pos_sam, 1, 
              function(x){quantile(x, probs = c(0.025, 0.975))})
  data.frame(
    mean = row_means,
    lower = CI[1, ],
    upper = CI[2, ]
  )
}

# check y prediction # 
lp_data <- calculate_means_and_cis(pos_sam_LP$pred_y_U_stack_sam)
lp_data$x <- y[-ind_mod]
cover_lp <- round(sum(lp_data$lower<lp_data$x & lp_data$upper>lp_data$x) / 
  length(lp_data$x) *100, 1) 
label_lp <- paste0('stacking of pds ', cover_lp,"% coverage")
lp_data$source <- label_lp
lse_data <- calculate_means_and_cis(pos_sam_LSE$pred_y_U_stack_sam)
lse_data$x <-  y[-ind_mod]
cover_lse <- round(sum(lse_data$lower<lse_data$x & lse_data$upper>lse_data$x) / 
                    length(lse_data$x) *100, 1) 
label_lse <- paste0('stacking of means ', cover_lse, "% coverage")
lse_data$source <- label_lse
mcmc_data <- calculate_means_and_cis(pos_wy$y.ho.sample)
mcmc_data$x <-  y[-ind_mod]
cover_mcmc <- round(sum(mcmc_data$lower<mcmc_data$x & 
                          mcmc_data$upper>mcmc_data$x) / 
                     length(mcmc_data$x) *100, 1) 
label_mcmc <- paste0('MCMC ', cover_mcmc, "% coverage")
mcmc_data$source <- label_mcmc


# Combine data
combined_data <- rbind(lse_data, lp_data, mcmc_data)
combined_data$source <- factor(combined_data$source, 
                               levels = c(label_lse, label_lp, label_mcmc))

# Determine the common range for x and y if not already known
x_range <- range(combined_data$x, na.rm = TRUE)
y_range <- range(combined_data$mean, combined_data$lower, combined_data$upper, na.rm = TRUE)
# Plot
pts_y <- ggplot(combined_data, aes(x = x, y = mean)) +
  geom_point(size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ source) +  # Use 'free_x' to unify y range
  scale_x_continuous(limits = x_range) +
  scale_y_continuous(limits = y_range) +
  theme_bw() +
  labs(x = "y", y = "mean and 95%CI")
pts_y
ggsave(paste0("./sim/pics/y_U_95CIsim", sim_ind, "_r", r, ".png"),
       plot = pts_y,
       width = 8, height = 3, units = "in", dpi = 600)

# check observed w prediction # 
lp_data <- calculate_means_and_cis(
  pos_sam_LP$pred_w_obs_stack_sam + 
    rep(1,r*100) %*% t(pos_sam_LP$pred_beta_stack_sam[1,]))
lp_data$x <- w[ind_mod] + raw_data[[r]]$beta[1]
cover_lp <- round(sum(lp_data$lower<lp_data$x & lp_data$upper>lp_data$x) / 
                    length(lp_data$x) *100, 1) 
label_lp <- paste0('stacking of pds ', cover_lp,"% coverage")
lp_data$source <- label_lp
lse_data <- calculate_means_and_cis(
  pos_sam_LSE$pred_w_obs_stack_sam + 
    rep(1,r*100) %*% t(pos_sam_LSE$pred_beta_stack_sam[1,]))
lse_data$x <-  w[ind_mod] + raw_data[[r]]$beta[1]
cover_lse <- round(sum(lse_data$lower<lse_data$x & lse_data$upper>lse_data$x) / 
                     length(lse_data$x) *100, 1) 
label_lse <- paste0('stacking of means ', cover_lse, "% coverage")
lse_data$source <- label_lse
mcmc_data <- calculate_means_and_cis(
  pos_wy$w.recover.sample[ind_mod, ] +
    rep(1,length(ind_mod)) %*% t(r.1$p.beta.recover.samples[, 1]))
mcmc_data$x <-  w[ind_mod] + raw_data[[r]]$beta[1]
cover_mcmc <- round(sum(mcmc_data$lower<mcmc_data$x & 
                          mcmc_data$upper>mcmc_data$x) / 
                      length(mcmc_data$x) *100, 1) 
label_mcmc <- paste0('MCMC ', cover_mcmc, "% coverage")
mcmc_data$source <- label_mcmc

# Combine data
combined_data <- rbind(lse_data, lp_data, mcmc_data)
combined_data$source <- factor(combined_data$source, 
                               levels = c(label_lse, label_lp, label_mcmc))

# Determine the common range for x and y if not already known
x_range <- range(combined_data$x, na.rm = TRUE)
y_range <- range(combined_data$mean, combined_data$lower, combined_data$upper, 
                 na.rm = TRUE)
# Plot
pts_w_obs <- ggplot(combined_data, aes(x = x, y = mean)) +
  geom_point(size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ source, scales = 'free_x') +  # Use 'free_x' to unify y range
  scale_x_continuous(limits = x_range) +
  scale_y_continuous(limits = y_range) +
  theme_bw() +
  labs(x = "z + intercept on observed locations", y = "mean and 95%CI")
pts_w_obs
ggsave(paste0("./sim/pics/w_obs_95CIsim", sim_ind, "_r", r, ".png"),
       plot = pts_w_obs,
       width = 8, height = 3, units = "in", dpi = 600)


# check unobserved w prediction # 
lp_data <- calculate_means_and_cis(
  pos_sam_LP$pred_w_U_stack_sam +
    rep(1, 100) %*% t(pos_sam_LP$pred_beta_stack_sam[1,]))
lp_data$x <- w[-ind_mod] + raw_data[[r]]$beta[1]
cover_lp <- round(sum(lp_data$lower<lp_data$x & lp_data$upper>lp_data$x) / 
                    length(lp_data$x) *100, 1) 
label_lp <- paste0('stacking of pds ', cover_lp,"% coverage")
lp_data$source <- label_lp
lse_data <- calculate_means_and_cis(
  pos_sam_LSE$pred_w_U_stack_sam +
    rep(1, 100) %*% t(pos_sam_LSE$pred_beta_stack_sam[1,]))
lse_data$x <-  w[-ind_mod] + raw_data[[r]]$beta[1]
cover_lse <- round(sum(lse_data$lower<lse_data$x & lse_data$upper>lse_data$x) / 
                     length(lse_data$x) *100, 1) 
label_lse <- paste0('stacking of means ', cover_lse, "% coverage")
lse_data$source <- label_lse
mcmc_data <- calculate_means_and_cis(
  pos_wy$w.recover.sample[-ind_mod, ] + 
    rep(1, 100) %*% t(r.1$p.beta.recover.samples[, 1]))
mcmc_data$x <-  w[-ind_mod] + raw_data[[r]]$beta[1]
cover_mcmc <- round(sum(mcmc_data$lower<mcmc_data$x & 
                          mcmc_data$upper>mcmc_data$x) / 
                      length(mcmc_data$x) *100, 1) 
label_mcmc <- paste0('MCMC ', cover_mcmc, "% coverage")
mcmc_data$source <- label_mcmc

# Combine data
combined_data <- rbind(lse_data, lp_data, mcmc_data)
combined_data$source <- factor(combined_data$source, 
                               levels = c(label_lse, label_lp, label_mcmc))

# Determine the common range for x and y if not already known
x_range <- range(combined_data$x, na.rm = TRUE)
y_range <- range(combined_data$mean, combined_data$lower, combined_data$upper, 
                 na.rm = TRUE)
# Plot
pts_w_U <- ggplot(combined_data, aes(x = x, y = mean)) +
  geom_point(size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ source, scales = 'free_x') +  # Use 'free_x' to unify y range
  scale_x_continuous(limits = x_range) +
  scale_y_continuous(limits = y_range) +
  theme_bw() +
  labs(x = "z + intercept on unobserved locations", y = "mean and 95%CI")
pts_w_U
ggsave(paste0("./sim/pics/w_U_95CIsim", sim_ind, "_r", r, ".png"),
       plot = pts_w_U,
       width = 8, height = 3, units = "in", dpi = 600)


# check unobserved w prediction for fun MCMC vs lp # 
mcmc_data <- calculate_means_and_cis(
  pos_wy$w.recover.sample[-ind_mod, ] + 
    rep(1, 100) %*% t(r.1$p.beta.recover.samples[, 1]))
lp_data <- calculate_means_and_cis(
  pos_sam_LP$pred_w_U_stack_sam+
    rep(1, 100) %*% t(pos_sam_LP$pred_beta_stack_sam[1,]))
lp_data$x <- rowMeans(pos_wy$w.recover.sample[-ind_mod, ]+ 
                        rep(1, 100) %*% t(r.1$p.beta.recover.samples[, 1]))
cover_lp <- round(sum(lp_data$lower<lp_data$x & lp_data$upper>lp_data$x) / 
                    length(lp_data$x) *100, 1) 
label_lp <- paste0("stacking of pds")
lp_data$source <- label_lp
lp_data$mcmc_lower <- mcmc_data$lower
lp_data$mcmc_upper <- mcmc_data$upper
lse_data <- calculate_means_and_cis(
  pos_sam_LSE$pred_w_U_stack_sam + 
    rep(1, 100) %*% t(pos_sam_LSE$pred_beta_stack_sam[1,]))
lse_data$x <- lp_data$x
cover_lse <- round(sum(lse_data$lower<lse_data$x & lse_data$upper>lse_data$x) / 
                     length(lse_data$x) *100, 1) 
label_lse <- paste0("stacking of means")
lse_data$source <- label_lse
lse_data$mcmc_lower <- mcmc_data$lower
lse_data$mcmc_upper <- mcmc_data$upper

# Combine lp_data, lse_data, and mcmc_data
combined_data <- rbind(lp_data, lse_data)

# Determine the common range for x and y
x_range <- range(combined_data$x, combined_data$mcmc_lower, 
                 combined_data$mcmc_upper, na.rm = TRUE)
y_range <- range(combined_data$mean, combined_data$lower, 
                 combined_data$upper, na.rm = TRUE)

# Plot
pts_w_U_mcmc_compar <- ggplot(combined_data, aes(x = x, y = mean)) +
  geom_point(size = 0.2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, alpha = 0.2) +
  geom_errorbarh(aes(xmin = mcmc_lower, xmax = mcmc_upper), height = 0.05, 
                 alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ source) + coord_fixed(ratio = 1) +
  scale_x_continuous(limits = x_range) +
  scale_y_continuous(limits = y_range) +
  theme_bw() + 
  labs(x = "MCMC mean and 95%CI for z + intercept on unobserved locations", 
       y = "mean and 95%CI")

pts_w_U_mcmc_compar

ggsave(paste0("./sim/pics/w_U_95CI_MCMC_vs_stacking_sim", sim_ind, "_r", r, ".png"),
       plot = pts_w_U_mcmc_compar,
       width = 5.5, height = 3, units = "in", dpi = 600)








