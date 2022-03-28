rm(list = ls())
library(ggplot2)
library(MASS)
library(spBayes)
library(ggplot2)
library(cmdstanr)
library(geoR)
library(rbenchmark)
library("gridExtra")
source("utils.R")
#options(mc.cores = parallel::detectCores())

deltasq_grid <- c(0.1, 0.5, 1, 2)
phi_grid = c(3, 9, 15, 21)   #3/(0.6*sqrt(2)) to 3/(0.1*sqrt(2)) phi_grid = c(3, 13, 23, 33)
nu_grid = c(0.5, 1, 1.5, 1.75)

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
raw_data <- list() # record raw data
expect_w <- list() # save the weighted latent process 
expect_y <- list() # save the weighted prediction 
DIV_matrix <- matrix(NA, nrow = N_list, ncol = 12)
colnames(DIV_matrix) <- c("SPE_stack_LSE", "SPE_stack_LP", "SPE_M0", 
                          "SPE_MCMC", "ELPD_stack_LSE", "ELPD_stack_LP", 
                          "ELPD_M0", "ELPD_MCMC", "SPE_w_stack_LP", 
                          "SPE_w_stack_LSE", "SPE_w_M0", "SPE_w_MCMC")
rownames(DIV_matrix) <- paste(samplesize_ls) # check
run_time <- matrix(0, 6, ncol = N_list)
MCMC_par <- list() # record the thinned MCMC chains for hyperparameters


for(r in 1:N_list){ # repeat
  cat("\n", "samplesize:", samplesize_ls[r], "\t")
  seed = samplesize_ls[r]
  set.seed(seed)
  N <- samplesize_ls[r]
  coords <- cbind(runif(N), runif(N))
  X <- as.matrix(cbind(1, rnorm(N)))
  beta <- as.matrix(c(1, 2))
  sigma.sq <- 1
  tau.sq <- 0.3 #1
  nu <- 0.5 # 1.0
  phi <- 20 # 7
  
  D <- as.matrix(dist(coords))
  R <- matern(D, phi = 1/phi, kappa = nu)
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
    seed = seed, label = "LP")
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
                                    coords.ho = coords[-ind_mod, ])
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
                   "phi.Unif"=c(3, 21), "sigma.sq.IG"=c(2, 2),
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
}
summary(DIV_matrix)
(run_time[4, ] + run_time[5, ])/run_time[1, ]
(run_time[4, ] + run_time[5, ])/run_time[2, ]

type = c("stacking LSE", "stacking LP", "M0", "MCMC")

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

