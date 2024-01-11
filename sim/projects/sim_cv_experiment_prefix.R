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
source("utils.R") # utils2.R is the testing code #
#options(mc.cores = parallel::detectCores())

############################################################
######## the impact of candidate value for sigma^2 #########
############################################################
load("./sim_hoffman2/results/sim1_1.RData")
# default method #
sim_ind = 1
deltasq_grid <- c(0.1, 0.5, 1, 2) # c(0.1, 0.5, 1, 2)
phi_grid = c(3, 14, 25, 36)   #3.5 to 35.8 #old: c(3, 9, 15, 31) 
nu_grid = c(0.5, 1, 1.5, 1.75)
deltasq_grid <- pick_deltasq(E_sigmasq = raw_data[[1]]$sigma.sq, 
                              E_tausq = raw_data[[1]]$tau.sq, 
                              b = max(raw_data[[1]]$tau.sq, raw_data[[1]]$sigma.sq),
                              p_ls = c(0.2, 0.4, 0.6, 0.8))
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


CV_fit_LSE <- sp_stacking_K_fold(
  X = X[ind_mod, ], y = y[ind_mod], 
  coords = coords[ind_mod, ],
  deltasq_grid = deltasq_grid, phi_grid = phi_grid,
  nu_grid = nu_grid, priors = priors, 
  K_fold = K_fold, seed = seed, label = "LSE")
cbind(CV_fit_LSE$grid_all[CV_fit_LSE$wts>0.00001, ], 
      CV_fit_LSE$wts[CV_fit_LSE$wts>0.00001])

pos_sam_LSE <- 
  stacking_pos_sample(Stack_fit = CV_fit_LSE, L1 = 300, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 1)

CV_fit_LP <- sp_stacking_K_fold(
  X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
  deltasq_grid = deltasq_grid, phi_grid = phi_grid, 
  nu_grid = nu_grid,
  priors = priors, K_fold = K_fold,
  seed = seed, label = "LP", MC = FALSE)
cbind(CV_fit_LP$grid_all[CV_fit_LP$wts>0.00001, ], 
      CV_fit_LP$wts[CV_fit_LP$wts>0.00001])

pos_sam_LP <- 
  stacking_pos_sample(Stack_fit = CV_fit_LP, L1 = 300, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 2)

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

a1 <- 3/2
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

sum((res1$summary.fitted.values[
  inla.stack.index(stack1, "est")$data[-raw_data[[r]]$ind_mod], 
  "0.01quant"] < y[-raw_data[[r]]$ind_mod]) & 
      (res1$summary.fitted.values[
        inla.stack.index(stack1, "est")$data[-raw_data[[r]]$ind_mod], 
        "0.99quant"] > y[-raw_data[[r]]$ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
# 0.55 for sim1 r = 8 0.66 for sim2 r = 4

a2 <- 2
spde2 = inla.spde2.pcmatern(mesh, alpha = a2,
                            prior.range = c(prior.median.range, .5), 
                            prior.sigma = c(prior.median.sd, .5))
stack2 = inla.stack(tag='est',
                    # - Name (nametag) of the stack
                    # - Here: est for estimating
                    data=list(y=df$y),
                    effects=list(
                      # - The Model Components
                      s=1:spde2$n.spde, 
                      # - The first is 's' (for spatial)
                      data.frame(intercept=1, x=df$x)),
                    # - The second is all fixed effects
                    A=list(A, 1)
)

res2 = inla(y~x + f(s, model=spde2), data=inla.stack.data(stack2),
            family = family,
            control.family = control.family,
            control.predictor=list(A = inla.stack.A(stack2)),
            quantiles=c(0.01, 0.025, 0.1, 0.5, 0.9, 0.975, 0.99),
            control.mode = list(restart = T, theta = initial.theta))
sum((res2$summary.fitted.values[
  inla.stack.index(stack2, "est")$data[-raw_data[[r]]$ind_mod], 
  "0.01quant"] < y[-raw_data[[r]]$ind_mod]) & 
    (res2$summary.fitted.values[
      inla.stack.index(stack2, "est")$data[-raw_data[[r]]$ind_mod], 
      "0.99quant"] > y[-raw_data[[r]]$ind_mod]))/length(y[-raw_data[[r]]$ind_mod])
# 0.55 for sim1 r = 8 0.66 for sim2 r = 4
all_prefix_ls <- cbind(
  c(1/(exp(res1$joint.hyper$`Log precision for the Gaussian observations`)*
         exp(res1$joint.hyper$`log(Stdev) for s`)^2),
    1/(exp(res2$joint.hyper$`Log precision for the Gaussian observations`)*
         exp(res2$joint.hyper$`log(Stdev) for s`)^2)) ,
  c(1/exp(res1$joint.hyper$`log(Range) for s`), 
    1/exp(res2$joint.hyper$`log(Range) for s`)),
  c(rep(a1 - 1, length(res1$joint.hyper$`log(Range) for s`)), 
    rep(a2 - 1, length(res2$joint.hyper$`log(Range) for s`))))
colnames(all_prefix_ls) <- c("delatsq", "phi", "nu")
all_prefix_ls <- as.data.frame(all_prefix_ls)

CV_fit_LSE_I <- sp_stacking_K_fold3(
  X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
  all_prefix_ls = all_prefix_ls, priors = priors, 
  K_fold = K_fold, seed = seed, label = "LSE")
cbind(CV_fit_LSE_I$grid_all[CV_fit_LSE_I$wts>0.00001, ], 
      CV_fit_LSE_I$wts[CV_fit_LSE_I$wts>0.00001])

pos_sam_LSE_I <- 
  stacking_pos_sample(Stack_fit = CV_fit_LSE_I, L1 = 300, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 3)


CV_fit_LP_I <- sp_stacking_K_fold3(
  X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
  all_prefix_ls = all_prefix_ls, 
  priors = priors, K_fold = K_fold,
  seed = seed, label = "LP", MC = FALSE)
cbind(CV_fit_LP_I$grid_all[CV_fit_LP_I$wts>0.00001, ], 
      CV_fit_LP_I$wts[CV_fit_LP_I$wts>0.00001])

pos_sam_LP_I <- 
  stacking_pos_sample(Stack_fit = CV_fit_LP_I, L1 = 300, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 4)

pos_y_LP_I_CI <- 
  apply(pos_sam_LP_I$pred_y_U_stack_sam, 1, 
        function(x){quantile(x, probs = c(0.025, 0.975))})

sum((pos_y_LP_I_CI[1, ] < y[-ind_mod]) & 
      (pos_y_LP_I_CI[2, ] > y[-ind_mod]))/100

# 0.96 for sim1 r = 8 and 0.86 for sim2 r = 4 
# I tested several example and CI coverage improved with Stacking

# test 3: posterior sample from marginal distribution #
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
  stacking_pos_sample(Stack_fit = CV_fit_LSE_P, L1 = 300, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 5)

CV_fit_LP_P <- sp_stacking_K_fold3(
  X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
  all_prefix_ls = all_prefix_ls, 
  priors = priors, K_fold = K_fold,
  seed = seed, label = "LP", MC = FALSE)
cbind(CV_fit_LP_P$grid_all[CV_fit_LP_P$wts>0.00001, ], 
      CV_fit_LP_P$wts[CV_fit_LP_P$wts>0.00001])

pos_sam_LP_P <- 
  stacking_pos_sample(Stack_fit = CV_fit_LP_P, L1 = 300, L2 = 900, 
                      X.mod = X[ind_mod, ], y.mod = y[ind_mod], 
                      coords.mod = coords[ind_mod, ], priors = priors,
                      X.ho = X[-ind_mod, ], 
                      coords.ho = coords[-ind_mod, ], seed = 6)

# compare sigmasq #
draws_ls_sigmasq <- c()
draws_ls_sigmasq[[1]] <- MCMC_par[[r]][101:1000, "sigma.sq"]
draws_ls_sigmasq[[2]] <- pos_sam_LSE$sigmasq_sam
draws_ls_sigmasq[[3]] <- pos_sam_LSE_I$sigmasq_sam
draws_ls_sigmasq[[4]] <- pos_sam_LSE_P$sigmasq_sam
draws_ls_sigmasq[[5]] <- MCMC_par[[r]][101:1000, "sigma.sq"]
draws_ls_sigmasq[[6]] <- pos_sam_LP$sigmasq_sam
draws_ls_sigmasq[[7]] <- pos_sam_LP_I$sigmasq_sam
draws_ls_sigmasq[[8]] <- pos_sam_LP_P$sigmasq_sam
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
sigmasq_compar <- 
  hist_compar(draws_ls = draws_ls_sigmasq, 
              type_colors = c("#999999", "#69b3a2", "#404080", "#E69F00"),  #"#56B4E9",
              type_names = c("posterior", "default", "INLA+stacking", 
                             "MCMC+stacking"), 
              test_names = c("LSE", "LP"), 
              true_value = raw_data[[r]]$sigma.sq,
              yname = "sigmasq", bins = 45)

sigmasq_compar
ggsave(paste0("./sim/pics/sigmasq_prefix_compar", sim_ind, "_r", r, ".png"),
       plot = sigmasq_compar,
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

# compare beta2 #
draws_ls_b2 <- c()
draws_ls_b2[[1]] <- r.1$p.beta.recover.samples[101:1000, 2]
draws_ls_b2[[2]] <- pos_sam_LSE$pred_beta_stack_sam[2, ]
draws_ls_b2[[3]] <- pos_sam_LSE_I$pred_beta_stack_sam[2, ]
draws_ls_b2[[4]] <- pos_sam_LSE_P$pred_beta_stack_sam[2, ]
draws_ls_b2[[5]] <- r.1$p.beta.recover.samples[101:1000, 2]
draws_ls_b2[[6]] <- pos_sam_LP$pred_beta_stack_sam[2, ]
draws_ls_b2[[7]] <- pos_sam_LP_I$pred_beta_stack_sam[2, ]
draws_ls_b2[[8]] <- pos_sam_LP_P$pred_beta_stack_sam[2, ]
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
b2_compar <- 
  hist_compar(draws_ls = draws_ls_b2, 
              type_colors = c("#999999", "#69b3a2", "#404080", "#E69F00"),  #"#56B4E9",
              type_names = c("posterior", "default", "INLA+stacking", 
                             "MCMC+stacking"), 
              test_names = c("LSE", "LP"), 
              true_value = raw_data[[r]]$beta[2],
              yname = "beta 2", bins = 45)

b2_compar
ggsave(paste0("./sim/pics/beta2_prefix_compar", sim_ind, "_r", r, ".png"),
       plot = sigmasq_compar,
       width = 6.5, height = 3.5, units = "in", dpi = 600)


# check y #
pick_indi <- c(50, 90)
draws_ls1 <- c()
draws_ls1[[1]] <- pos_wy$y.ho.sample[pick_indi[1], 101:1000]
draws_ls1[[2]] <- pos_sam_LSE$pred_y_U_stack_sam[pick_indi[1], ]
draws_ls1[[3]] <- pos_sam_LSE_I$pred_y_U_stack_sam[pick_indi[1], ]
draws_ls1[[4]] <- pos_sam_LSE_P$pred_y_U_stack_sam[pick_indi[1], ]
draws_ls1[[5]] <- pos_wy$y.ho.sample[pick_indi[1], 101:1000]
draws_ls1[[6]] <- pos_sam_LP$pred_y_U_stack_sam[pick_indi[1], ]
draws_ls1[[7]] <- pos_sam_LP_I$pred_y_U_stack_sam[pick_indi[1], ]
draws_ls1[[8]] <- pos_sam_LP_P$pred_y_U_stack_sam[pick_indi[1], ]

individual_y_compar1 <- 
  hist_compar(draws_ls = draws_ls1, 
              type_colors = c("#999999", "#69b3a2", "#404080", "#E69F00"), 
              type_names = c("posterior", "default", "INLA+Stacking", 
                             "MCMC+Stacking"), 
              test_names = c("LSE", "LP"), 
              true_value = y[-ind_mod][pick_indi[1]],
              yname = "y",
              INLA_CI = c(unlist(res1$summary.fitted.values[
                inla.stack.index(stack1, "est")$
                  data[-raw_data[[r]]$ind_mod][pick_indi[1]], 
                c("0.01quant", "0.99quant")]))
              # INLA_CI = c(unlist(res2$summary.fitted.values[
              #   inla.stack.index(stack2, "est")$
              #     data[-raw_data[[r]]$ind_mod][pick_indi[1]], 
              #   c("0.01quant", "0.99quant")]))
              )

res2$summary.fitted.values[r*100+pick_indi[2], c("0.01quant", "0.99quant")]
draws_ls2 <- c()
draws_ls2[[1]] <- pos_wy$y.ho.sample[pick_indi[2], 101:1000]
draws_ls2[[2]] <- pos_sam_LSE$pred_y_U_stack_sam[pick_indi[2], ]
draws_ls2[[3]] <- pos_sam_LSE_I$pred_y_U_stack_sam[pick_indi[2], ]
draws_ls2[[4]] <- pos_sam_LSE_P$pred_y_U_stack_sam[pick_indi[2], ]
draws_ls2[[5]] <- pos_wy$y.ho.sample[pick_indi[2], 101:1000]
draws_ls2[[6]] <- pos_sam_LP$pred_y_U_stack_sam[pick_indi[2], ]
draws_ls2[[7]] <- pos_sam_LP_I$pred_y_U_stack_sam[pick_indi[2], ]
draws_ls2[[8]] <- pos_sam_LP_P$pred_y_U_stack_sam[pick_indi[2], ]

individual_y_compar2 <- 
  hist_compar(draws_ls = draws_ls2, 
              type_colors = c("#999999", "#69b3a2", "#404080", "#E69F00"), 
              type_names = c("posterior", "default", "INLA+stacking", 
                             "MCMC+stacking"), 
              test_names = c("LSE", "LP"), 
              true_value = y[-ind_mod][pick_indi[2]],
              yname = "y",
              INLA_CI = c(unlist(res1$summary.fitted.values[
                inla.stack.index(stack1, "est")$
                  data[-raw_data[[r]]$ind_mod][pick_indi[2]], 
                c("0.01quant", "0.99quant")]))
              # INLA_CI = c(unlist(res2$summary.fitted.values[
              #   inla.stack.index(stack2, "est")$
              #     data[-raw_data[[r]]$ind_mod][pick_indi[2]], 
              #   c("0.01quant", "0.99quant")]))
              )

print(individual_y_compar1)
print(individual_y_compar2)
ggsave(paste0("./sim/pics/indi50_prefix_compar", sim_ind, "ICI_r", r, ".png"),
       plot = individual_y_compar1,
       width = 6.5, height = 3.5, units = "in", dpi = 600)
ggsave(paste0("./sim/pics/indi90_prefix_compar", sim_ind, "ICI_r", r, ".png"),
       plot = individual_y_compar2,
       width = 6.5, height = 3.5, units = "in", dpi = 600)

# check w #
pick_indi <- c(50, 90)
draws_ls1 <- c()
draws_ls1[[1]] <- pos_wy$w.recover.sample[r*100+pick_indi[1], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls1[[2]] <- pos_sam_LSE$pred_w_U_stack_sam[pick_indi[1], ] +
  pos_sam_LSE$pred_beta_stack_sam[1, ]
draws_ls1[[3]] <- pos_sam_LSE_I$pred_w_U_stack_sam[pick_indi[1], ] +
  pos_sam_LSE_I$pred_beta_stack_sam[1, ]
draws_ls1[[4]] <- pos_sam_LSE_P$pred_w_U_stack_sam[pick_indi[1], ]+
  pos_sam_LSE_P$pred_beta_stack_sam[1, ]
draws_ls1[[5]] <- pos_wy$w.recover.sample[r*100+pick_indi[1], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls1[[6]] <- pos_sam_LP$pred_w_U_stack_sam[pick_indi[1], ] +
  pos_sam_LP$pred_beta_stack_sam[1, ]
draws_ls1[[7]] <- pos_sam_LP_I$pred_w_U_stack_sam[pick_indi[1], ] +
  pos_sam_LP_I$pred_beta_stack_sam[1, ]
draws_ls1[[8]] <- pos_sam_LP_P$pred_w_U_stack_sam[pick_indi[1], ] +
  pos_sam_LP_P$pred_beta_stack_sam[1, ]

individual_w_compar1 <- 
  hist_compar(draws_ls = draws_ls1, 
              type_colors = c("#999999", "#69b3a2", "#404080", "#E69F00"), 
              type_names = c("posterior", "default", "INLA+Stacking", 
                             "MCMC+Stacking"), 
              test_names = c("LSE", "LP"), 
              true_value = w[-ind_mod][pick_indi[1]]+raw_data[[r]]$beta[1],
              yname = "w+beta1"
              # INLA_CI = c(unlist(res2$summary.fitted.values[
              #   inla.stack.index(stack2, "est")$
              #     data[-raw_data[[r]]$ind_mod][pick_indi[1]], 
              #   c("0.01quant", "0.99quant")]))
  )

res2$summary.fitted.values[r*100+pick_indi[2], c("0.01quant", "0.99quant")]
draws_ls2 <- c()
draws_ls2[[1]] <- pos_wy$w.recover.sample[r*100+pick_indi[2], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls2[[2]] <- pos_sam_LSE$pred_w_U_stack_sam[pick_indi[2], ] +
  pos_sam_LSE$pred_beta_stack_sam[1, ]
draws_ls2[[3]] <- pos_sam_LSE_I$pred_w_U_stack_sam[pick_indi[2], ] +
  pos_sam_LSE_I$pred_beta_stack_sam[1, ]
draws_ls2[[4]] <- pos_sam_LSE_P$pred_w_U_stack_sam[pick_indi[2], ]+
  pos_sam_LSE_P$pred_beta_stack_sam[1, ]
draws_ls2[[5]] <- pos_wy$w.recover.sample[r*100+pick_indi[2], 101:1000]+
  r.1$p.beta.recover.samples[101:1000, 1]
draws_ls2[[6]] <- pos_sam_LP$pred_w_U_stack_sam[pick_indi[2], ] +
  pos_sam_LP$pred_beta_stack_sam[1, ]
draws_ls2[[7]] <- pos_sam_LP_I$pred_w_U_stack_sam[pick_indi[2], ] +
  pos_sam_LP_I$pred_beta_stack_sam[1, ]
draws_ls2[[8]] <- pos_sam_LP_P$pred_w_U_stack_sam[pick_indi[2], ] +
  pos_sam_LP_P$pred_beta_stack_sam[1, ]

individual_w_compar2 <- 
  hist_compar(draws_ls = draws_ls2, 
              type_colors = c("#999999", "#69b3a2", "#404080", "#E69F00"), 
              type_names = c("posterior", "default", "INLA+stacking", 
                             "MCMC+stacking"), 
              test_names = c("LSE", "LP"), 
              true_value = w[-ind_mod][pick_indi[2]]+raw_data[[r]]$beta[1],
              yname = "w+beta1"
              # INLA_CI = c(unlist(res2$summary.fitted.values[
              #   inla.stack.index(stack2, "est")$
              #     data[-raw_data[[r]]$ind_mod][pick_indi[2]], 
              #   c("0.01quant", "0.99quant")]))
  )

print(individual_w_compar1)
print(individual_w_compar2)
ggsave(paste0("./sim/pics/indi50_prefix_incpw_compar", sim_ind, "_r", r, ".png"),
       plot = individual_w_compar1,
       width = 6.5, height = 3.5, units = "in", dpi = 600)
ggsave(paste0("./sim/pics/indi90_prefix_incpw_compar", sim_ind, "_r", r, ".png"),
       plot = individual_w_compar2,
       width = 6.5, height = 3.5, units = "in", dpi = 600)



