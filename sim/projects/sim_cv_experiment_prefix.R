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
source("utils2.R") # utils2.R is the testing code #
#options(mc.cores = parallel::detectCores())

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

# test 2: Empirical method: semivariogram + indep deltasq #
load("./sim_hoffman2/results/sim1_1.RData")
## Variogram ##
library(geoR)
library(fields)
### variogram of raw data and residuals ###
r = 2
coords_train <- raw_data[[r]]$coords[raw_data[[r]]$ind_mod, ]
lm_fit <- lm(raw_data[[r]]$y[raw_data[[r]]$ind_mod]~
               raw_data[[r]]$X[raw_data[[r]]$ind_mod, 2])
res <- residuals(lm_fit)
max.dist=0.6*max(rdist(coords_train))
bins=20
vario.resid <- variog(coords=coords_train,
                          data=res,
                          uvec=(seq(0.01, max.dist, length=bins)))
plot(vario.resid)
vfit_wls=variofit(vario.resid, 
                  ini.cov.pars=c(1, 0.4), 
                  nugget=1, 
                  fix.nugget=FALSE, fix.kappa = FALSE,
                  cov.model='matern', 
                  weights='cressie')
vfit_wls
plot(vario.resid)
lines(vfit_wls, col="purple", lwd=1.5)


eff_range <- c(0.2, 0.4, 0.6, 0.8)
phi_ls <- decay_est(eff_range, nu_grid)
phi_nu_ls <- cbind(c(phi_ls), rep(nu_grid, each = length(eff_range))) # put all phi and nu candidate value here
colnames(phi_nu_ls) = c("phi", "nu")


# test 3: posterior sample from marginal distribution #
input_id = 1
load("./sim_hoffman2/results/sim1_1.RData")

# test 4: INLA #
## fit INLA ##
library(INLA)
r = 8
df = data.frame(y=c(raw_data[[r]]$y[raw_data[[r]]$ind_mod], rep(NA, 100)), 
                locx=raw_data[[r]]$coords[, 1], 
                locy=raw_data[[r]]$coords[, 2], x = raw_data[[r]]$X[, 2])
summary(df)

fake.locations = matrix(c(0,0,1,1, 0, 1, 1, 0), nrow = 4, byrow = T)
mesh = inla.mesh.2d(loc = fake.locations, max.edge=c(0.08, 1))
mesh$n

A = inla.spde.make.A(mesh=mesh, loc=data.matrix(df[ , c('locx', 'locy')]))
dim(A);

par(mfrow=c(1,1))
plot(mesh, asp=1)
points(df[ , c('locx', 'locy')], col='red', lwd=.1)

a <- 3/2
prior.median.sd = 1; prior.median.range = 1/7
spde = inla.spde2.pcmatern(mesh, alpha = a,
                           prior.range = c(prior.median.range, .5), 
                           prior.sigma = c(prior.median.sd, .5))
stack = inla.stack(tag='est',
                   # - Name (nametag) of the stack
                   # - Here: est for estimating
                   data=list(y=df$y),
                   effects=list(
                     # - The Model Components
                     s=1:spde$n.spde, 
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

res = inla(y~x + f(s, model=spde), data=inla.stack.data(stack),
           family = family,
           control.family = control.family,
           control.predictor=list(A = inla.stack.A(stack)),
           quantiles=c(0.01, 0.1, 0.5, 0.9, 0.99),
           control.mode = list(restart = T, theta = initial.theta))

summary(res)
all_prefix_ls <- cbind(
  1/(exp(res$joint.hyper$`Log precision for the Gaussian observations`)*
       exp(res$joint.hyper$`log(Stdev) for s`)^2),
  1/exp(res$joint.hyper$`log(Range) for s`),
  rep(a - 1, length(res$joint.hyper$`log(Range) for s`)))
colnames(all_prefix_ls) <- c("delatsq", "phi", "nu")
all_prefix_ls <- as.data.frame(all_prefix_ls)


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
colnames(DIV_matrix2) <- c("SPE_stack_LSE_E", "SPE_stack_LP_E",
                          "SPE_stack_LSE_P", "SPE_stack_LP_P",
                          "ELPD_stack_LSE_E", "ELPD_stack_LP_E",
                          "ELPD_stack_LSE_P", "ELPD_stack_LP_P",
                          "SPE_w_stack_LSE_E", "SPE_w_stack_LP_E",
                          "SPE_w_stack_LSE_P", "SPE_w_stack_LP_P")
rownames(DIV_matrix2) <- paste(samplesize_ls) # check
run_time2 <- matrix(0, 4, ncol = N_list)
rownames(run_time2) <- c("Stack_LSE_E", "Stack_LP_E",
                        "Stack_LSE_P", "Stack_LP_P")
#MCMC_par <- list() # record the thinned MCMC chains for hyperparameters

t <- proc.time()
for(r in 1:N_list){ # repeat
  cat("\n", "samplesize:", samplesize_ls[r], "\t")
  seed = samplesize_ls[r] + input_id
  
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
  CV_fit_LSE <- sp_stacking_K_fold2(
    X = X[ind_mod, ], y = y[ind_mod], 
    coords = coords[ind_mod, ],
    deltasq_grid = deltasq_grid, phi_nu_ls = phi_nu_ls, priors = priors, 
    K_fold = K_fold, seed = seed, label = "LSE")
  weights_M_LSE[, r] <- CV_fit_LSE$wts
  run_time2["Stack_LSE_E", r] <- CV_fit_LSE$time[3]
  
  CV_fit_LP <- sp_stacking_K_fold2(
    X = X[ind_mod, ], y = y[ind_mod], coords = coords[ind_mod, ],
    deltasq_grid = deltasq_grid, phi_nu_ls = phi_nu_ls, 
    priors = priors, K_fold = K_fold,
    seed = seed, label = "LP", MC = FALSE)
  weights_M_LP[, r] <- CV_fit_LP$wts
  run_time2["Stack_LP_E", r] <- CV_fit_LP$time[3]
  
  
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
  DIV_matrix2[r, "SPE_stack_LSE_E"] <- mean((y_pred_stack_LSE - y[-ind_mod])^2)
  y_pred_stack_LP = y_pred_grid %*% CV_fit_LP$wts
  DIV_matrix2[r, "SPE_stack_LP_E"] <- mean((y_pred_stack_LP - y[-ind_mod])^2)
  
  w_expect_stack_LSE = w_expect_grid %*% CV_fit_LSE$wts
  DIV_matrix2[r, "SPE_w_stack_LSE_E"] <- mean((w_expect_stack_LSE - w)^2)
  w_expect_stack_LP = w_expect_grid %*% CV_fit_LP$wts
  DIV_matrix2[r, "SPE_w_stack_LP_E"] <- mean((w_expect_stack_LP - w)^2)
  
  
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
  DIV_matrix2[r, "ELPD_stack_LSE_E"] = mean(log(exp(lp_pred_grid) %*% CV_fit_LSE$wts))
  DIV_matrix2[r, "ELPD_stack_LP_E"] = mean(log(exp(lp_pred_grid) %*% CV_fit_LP$wts))
  
  
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

type = c("stacking LSE E", "stacking LP E", "stacking LSE P", "stacking LP P")

dat_check <- data.frame(N_sample = rep(samplesize_ls, length(type)),
                        SPE = c(DIV_matrix[, c("SPE_stack_LSE", "SPE_stack_LP", 
                                               "SPE_M0", "SPE_MCMC")]),
                        SPE_w = c(DIV_matrix[, c("SPE_w_stack_LSE", 
                                                 "SPE_w_stack_LP", "SPE_w_M0", 
                                                 "SPE_w_MCMC")]),
                        ELPD = c(DIV_matrix[, c("ELPD_stack_LSE", 
                                                "ELPD_stack_LP", "ELPD_M0", 
                                                "ELPD_MCMC")]),
                        label = rep(type, each = N_list))

p_SPE <- ggplot(dat_check) +
  geom_line(aes(N_sample, SPE, group = label,
                color = label)) +theme_bw()
p_SPE

p_SPE_w <- ggplot(dat_check) +
  geom_line(aes(N_sample, SPE_w, group = label,
                color = label))+theme_bw()
p_SPE_w 
p_ELPD <- ggplot(dat_check) +
  geom_line(aes(N_sample, ELPD, group = label,
                color = label))+theme_bw()
p_ELPD

p_summary <- grid.arrange(p_SPE,  p_SPE_w, p_ELPD,
                          ncol = 1, nrow = 3)
ggsave("./sim/pics/CVexperiment_lit.png", plot = p_summary, 
       width = 6, height = 4, units = "in")
save(dat_check, samplesize_ls, weights_M_LSE,
     weights_M_LP, file = "./sim/results/CVexperiment_lit.Rdata")
#load("./sim/results/CVexperiment_lit.Rdata")
weight_stack <- 
  data.frame(weights = c(weights_M_LSE), 
             model = c(rep(paste0("deltasq:", 
                                  round(CV_fit_LSE$grid_all[, 1], digits = 3), 
                                  " phi: ",
                                  round(CV_fit_LSE$grid_all[, 2], digits = 3)),
                           N_list)))

p1 <- ggplot(data = weight_stack, aes(x = weights)) + 
  geom_histogram() + facet_wrap(~model)
p1

summary(colSums(weights_M_LSE))

expect_w_phi <- matrix(0, nrow = length(phi_grid), ncol = N_list)
for(i in 1:length(phi_grid)){
  expect_w_phi[i, ] <- colSums(weights_M_LSE[which(CV_fit_LSE$grid_all$phi == phi_grid[i]), ])
}
summary(c(t(expect_w_phi) %*%phi_grid))
hist(t(expect_w_phi) %*%phi_grid )
plot(phi_grid, rowMeans(expect_w_phi))

expect_w_deltasq <- matrix(0, nrow = length(deltasq_grid), ncol = N_list)
for(i in 1:length(deltasq_grid)){
  expect_w_deltasq[i, ] <- 
    colSums(weights_M[which(CV_fit_LSE$grid_all$deltasq == deltasq_grid[i]),])
}
summary(c(t(expect_w_deltasq) %*% deltasq_grid))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1000  0.5164  0.7113  0.8436  1.0581  2.0000
hist(t(expect_w_deltasq) %*% deltasq_grid)
plot(deltasq_grid, rowMeans(expect_w_deltasq))

p1 <- mcmc_trace(fit_mcmc$draws("lp__")) 
print(p1)


save(dat_check, samplesize_ls, weights_M_LSE,
     weights_M_LP, file = "./sim/results/CVexperiment_lit.Rdata")


# save(weights_M_LSE, weights_M_LP, raw_data, expect_w, expect_y, DIV_matrix,
#      run_time, file = "./sim/results/1.Rdata")

## check the plots of latent process ##
library(coda)
library(spBayes)
library(MBA)
library(fields)
library(classInt)
library(RColorBrewer)
library(sp)

r = 3
h <- 12
surf.raw <- mba.surf(cbind(raw_data[[r]]$coords, raw_data[[r]]$w), no.X = 300, 
                     no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.LSE <- mba.surf(cbind(raw_data[[r]]$coords, expect_w[[r]][, "LSE"]), 
                     no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.LP <- mba.surf(cbind(raw_data[[r]]$coords, expect_w[[r]][, "LP"]), 
                    no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.M0 <- mba.surf(cbind(raw_data[[r]]$coords, expect_w[[r]][, "M0"]), 
                    no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.MCMC <- mba.surf(cbind(raw_data[[r]]$coords, expect_w[[r]][, "MCMC"]), 
                      no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est

surf.brks <- classIntervals(surf.raw$z, 500, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1.13)
zlim <- range(c(surf.raw[["z"]], surf.LSE[["z"]], surf.LP[["z"]], 
                surf.M0[["z"]], surf.MCMC[["z"]]))

# size for the mapping of w               
width <- 360
height <- 360
pointsize <- 16

# setEPS()
# postscript("./pic/map-w-true.eps")
par(mfrow = c(2, 3))
i <- as.image.SpatialGridDataFrame(surf.raw)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "raw")
axis(2, las=1)
axis(1)
image.plot(i, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
# dev.off()

i1 <- as.image.SpatialGridDataFrame(surf.LSE)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "stacking-LSE") 
axis(2, las=1)
axis(1)
image.plot(i1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i2 <- as.image.SpatialGridDataFrame(surf.LP)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "stacking-LP") 
axis(2, las=1)
axis(1)
image.plot(i2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i3 <- as.image.SpatialGridDataFrame(surf.M0)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "M0") 
axis(2, las=1)
axis(1)
image.plot(i3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i4 <- as.image.SpatialGridDataFrame(surf.MCMC)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "MCMC") 
axis(2, las=1)
axis(1)
image.plot(i4, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

# check predicted w
surf.raw <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                           raw_data[[r]]$w[-raw_data[[r]]$ind_mod]), no.X = 300, 
                     no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.LSE <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                           expect_w[[r]][-raw_data[[r]]$ind_mod, "LSE"]), 
                     no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.LP <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ],
                          expect_w[[r]][-raw_data[[r]]$ind_mod, "LP"]), 
                    no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.M0 <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                          expect_w[[r]][-raw_data[[r]]$ind_mod, "M0"]), 
                    no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.MCMC <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                            expect_w[[r]][-raw_data[[r]]$ind_mod, "MCMC"]), 
                      no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est

surf.brks <- classIntervals(surf.raw$z, 500, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1.13)
zlim <- range(c(surf.raw[["z"]], surf.LSE[["z"]], surf.LP[["z"]], 
                surf.M0[["z"]], surf.MCMC[["z"]]))

par(mfrow = c(2, 3))
i <- as.image.SpatialGridDataFrame(surf.raw)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "raw")
axis(2, las=1)
axis(1)
image.plot(i, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
# dev.off()

i1 <- as.image.SpatialGridDataFrame(surf.LSE)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "stacking-LSE") 
axis(2, las=1)
axis(1)
image.plot(i1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i2 <- as.image.SpatialGridDataFrame(surf.LP)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "stacking-LP") 
axis(2, las=1)
axis(1)
image.plot(i2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i3 <- as.image.SpatialGridDataFrame(surf.M0)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "M0") 
axis(2, las=1)
axis(1)
image.plot(i3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i4 <- as.image.SpatialGridDataFrame(surf.MCMC)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "MCMC") 
axis(2, las=1)
axis(1)
image.plot(i4, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)


### check the predicted y ##
#r = 2
surf.raw <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                           raw_data[[r]]$y[-raw_data[[r]]$ind_mod]), no.X = 300, 
                     no.Y = 300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.LSE <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                           expect_y[[r]][, "LSE"]), 
                     no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.LP <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                          expect_y[[r]][, "LP"]), 
                    no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.M0 <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                          expect_y[[r]][, "M0"]), 
                    no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est
surf.MCMC <- mba.surf(cbind(raw_data[[r]]$coords[-raw_data[[r]]$ind_mod, ], 
                            expect_y[[r]][, "MCMC"]), 
                      no.X=300, no.Y=300, exten = TRUE, sp = TRUE, h = h)$xyz.est

surf.brks <- classIntervals(surf.raw$z, 500, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0, 1.13)
zlim <- range(c(surf.raw[["z"]], surf.LSE[["z"]], surf.LP[["z"]], 
                surf.M0[["z"]], surf.MCMC[["z"]]))

# size for the mapping of w               
width <- 360
height <- 360
pointsize <- 16

# setEPS()
# postscript("./pic/map-w-true.eps")
par(mfrow = c(2, 3))
i <- as.image.SpatialGridDataFrame(surf.raw)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "raw")
axis(2, las=1)
axis(1)
image.plot(i, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
# dev.off()

i1 <- as.image.SpatialGridDataFrame(surf.LSE)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "stacking-LSE") 
axis(2, las=1)
axis(1)
image.plot(i1, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i2 <- as.image.SpatialGridDataFrame(surf.LP)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "stacking-LP") 
axis(2, las=1)
axis(1)
image.plot(i2, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i3 <- as.image.SpatialGridDataFrame(surf.M0)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "M0") 
axis(2, las=1)
axis(1)
image.plot(i3, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

i4 <- as.image.SpatialGridDataFrame(surf.MCMC)
plot(raw_data[[r]]$coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", 
     xlab="x", main = "MCMC") 
axis(2, las=1)
axis(1)
image.plot(i4, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)



### old code ###
### fit with stan ###
# #######################
# ## predict with MCMC ##
# #######################
# 
# mod <- cmdstan_model(stan_file = "./sim/projects/spatial_fit.stan")
# 
# ### setup ###
# dat = list(N_mod = (N - N_ho), P = ncol(X), N_ho = N_ho, N = N, X = X,
#            y = y, coords = coords,
#            mu_beta = priors$mu_beta,
#            V_beta = chol2inv(chol(priors$inv_V_beta)),
#            a_sigma = priors$a_sigma, b_sigma = priors$b_sigma,
#            a_tau = priors$a_sigma, b_tau = priors$b_sigma,
#            a_phi = min(phi_grid), b_phi = max(phi_grid))
# 
# init_fun <- function() list(phi = phi, tau = sqrt(tau.sq),
#                             sigma = sqrt(sigma.sq), beta = beta)
# 
# fit_mcmc <- mod$sample(data = dat,
#                        seed = 123,
#                        chains = 2,
#                        refresh = 500,
#                        init = init_fun,
#                        parallel_chains = 2,
#                        iter_warmup = 1000,
#                        iter_sampling = 100)
# posterior <- as.array(fit_mcmc$draws())
# dim(posterior)
# mcmc_trace(posterior, pars = c("tau", "sigma", "phi"))
# fit_mcmc
# 
# y_pred_MCMC <- fit_mcmc$summary("y_pred", "mean")$mean
# w_expect_MCMC <- fit_mcmc$summary("w", "mean")$mean
# SPE_MCMC[r] <- mean((y_pred_MCMC - y[(N - N_ho + 1):N])^2)
# SPE_w_MCMC[r] <- mean((w_expect_MCMC - w)^2)
# ELPD_MCMC[r] <- mean(fit_mcmc$draws("lp_pred"))

### fit with spBayes ###
n.samples <- 4000
starting <- list("phi"=3/0.5, "sigma.sq"=1, "tau.sq"=1, "nu" = 0.5)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1, "nu" = 0.1)
priors.1 <- list("beta.Norm"=list(rep(0, ncol(X)), solve(priors$inv_V_beta)),
                 "phi.Unif"=c(3, 33), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 2), "nu.unif" = c(0.25, 2))
# priors.2 <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
#                  "sigma.sq.IG"=c(2, 2), "tau.sq.IG"=c(2, 0.1))
cov.model <- "matern"
n.report <- 5000
verbose <- TRUE
m.1 <- spLM(y[ind_mod]~X[ind_mod, ]-1, coords=coords[ind_mod, ],
            starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)


## check MCMC trace plot ##
library(bayesplot)
color_scheme_set("blue")
mcmc_trace(m.1$p.theta.samples, n_warmup = 0.5*n.samples)

# ## recover beta ##
# r.1 <- spRecover(m.1, get.w = FALSE, start = 0.5*n.samples, thin = 10,
#                  n.report =  500)
# 
# ## recover latent process on all locations ##
# ## compute expected response and latent process ##
# MCMC_out <- expects_MCMC(theta.recover = r.1$p.theta.recover.samples, 
#         beta.recover = r.1$p.beta.recover.samples, 
#         y.mod = y[ind_mod], X.mod = X[ind_mod, ], coords.mod = coords[ind_mod, ], 
#         X.ho = X[-ind_mod, ], y.ho = y[-ind_mod], coords.ho = coords[-ind_mod, ])

# test solvers, MOSEK is installed in julia#
#RMOSEKDIR = "/home/luzhang/.julia/packages/Mosek/au3Cq/deps/src/mosek/9.3/tools/platform/linux64x86/rmosek"
#source(paste0(RMOSEKDIR, "/builder.R"))
#attachbuilder(what_mosek_bindir="/home/luzhang/.julia/packages/Mosek/au3Cq/deps/src/mosek/9.3/tools/platform/linux64x86/bin", pos=2L, name="Rmosek:builder", warn.conflicts=TRUE)
#install.rmosek()

# library(Rmosek)
# result <- psolve(prob, solver = "ECOS_BB", verbose = TRUE)
# round(result$getValue(w)[,1], digits = 3)
# sum(log(exp_lpd_point%*%result$getValue(w)))
# #[1] 43.18918
# 
# results2 <- psolve(prob, solver = "MOSEK", verbose = TRUE)
# round(results2$getValue(w)[,1], digits = 3)
# sum(log(exp_lpd_point%*%results2$getValue(w)))
# #[1] 43.18921  win
# 
# results3 <- psolve(prob, solver = "SCS", verbose = TRUE)
# round(results3$getValue(w)[,1], digits = 3)
# sum(log(exp_lpd_point%*%results3$getValue(w)))
# #[1] 42.40507
# 
# result_old <- stacking_weights_old(lp_expect)
# sum(log(exp_lpd_point%*%result_old))
# #[1] 39.34675
# sum(log(exp_lpd_point%*%rep(1/64, 64)))
# 
# result_test <- stacking_weights(lp_expect)
# sum(log(exp_lpd_point%*%result_test))

