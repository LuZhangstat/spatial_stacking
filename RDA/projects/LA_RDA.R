rm(list = ls())
library(reticulate)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)

# load raw data
np <- import("numpy")
LA_AOD <- np$load("./RDA/LAdemoData/LA_AOD_17000-997-zeroNotValid.npy")
LA_testmask <- np$load("./RDA/LAdemoData/LA_testmask_4146.npy")
LA_trainmask <- np$load("./RDA/LAdemoData/LA_trainmask_11857.npy")
LA_xfea3 <- np$load("./RDA/LAdemoData/LA_xfea3_17000-997.npy")
LA_xfea64 <- np$load("./RDA/LAdemoData/LA_xfea64_17000-997.npy")

## EDA ##
# transform the data from matrix into a dataframe
n_row <- nrow(LA_AOD)
n_col <- ncol(LA_AOD)
x_ind <- rep(1:n_col, each = n_row)
y_ind <- rep(n_row:1, time = n_col)    # y from high to low

df_aod <- data.frame(AOD = c(LA_AOD), x = x_ind, y = y_ind)
df_aod %>% glimpse()
summary(df_aod)
qqnorm(df_aod$AOD)        # the data has 0 for missing value
sum(df_aod$AOD == 0)      # there are 997 pixels have AOD values equal to 0. 
df_aod <- df_aod %>% filter(AOD > 0)   # Remove 0 from the data
df_aod %>% glimpse()
# our data has 16,003 pixles

# Add the 64 features into the data
covariates_matrix <- array(data = LA_xfea64, 
                           dim = c(dim(LA_xfea64)[1] * dim(LA_xfea64)[2], 
                                   dim(LA_xfea64)[3]))
covariates_df <- as.data.frame(covariates_matrix)
covariates_df$x <- x_ind
covariates_df$y <- y_ind

# Combine with df_aod
combined_data <- df_aod %>%
  left_join(covariates_df, by = c("x", "y"))

# Check the combined data frame
glimpse(combined_data)

# Convert the mask to a data frame
mask_df <- data.frame(mask = c(LA_testmask), x = x_ind, y = y_ind)
combined_data_train <- combined_data %>% 
  left_join(mask_df, by = c("x", "y")) %>%filter(mask == FALSE)
glimpse(combined_data_train)
combined_data_test <- combined_data %>% 
  left_join(mask_df, by = c("x", "y")) %>%filter(mask == TRUE)
glimpse(combined_data_test)
# training: 11,857; testing: 4,146

## Some simple data analysis ##
# linear regression:
lm_fit <- lm(AOD~x+y+., data = combined_data_train)
summary(lm_fit)   # Adjusted R-squared = 0.7387
lm_fit2 <- lm(log(AOD)~x+y+., data = combined_data_train)
summary(lm_fit2)   # Adjusted R-squared = 0.7336
# check residual
AOD_residual1 <- residuals(lm_fit); hist(AOD_residual1)
qqnorm(scale(AOD_residual1)); abline(a=0,b=1)
AOD_residual2 <- residuals(lm_fit2); hist(AOD_residual2)
qqnorm(scale(AOD_residual2)); abline(a=0,b=1)
#' The residual of log AOD behaves similar to Gaussian
#' We use Log_transform here in the model fitting.
combined_data_train$logAOD <- log(combined_data_train$AOD)


## Visualization ##
# The map of the raw AOD data in central LA ##
aod_p <- ggplot(combined_data, aes(x = x, y = y, fill = log(AOD))) + 
  geom_tile() + 
  scale_fill_gradientn(colors = viridis::viridis(6)) +
  theme_minimal() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())
aod_p

# check the plot of covariates
V_p <- ggplot(combined_data, aes(x = x, y = y, fill = V3)) + 
  geom_tile() + 
  scale_fill_gradientn(colors = viridis::viridis(6)) +
  theme_minimal() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())
V_p

# check the residual #
res_dat <- data.frame(residual = AOD_residual2, x = combined_data_train$x,
                      y = combined_data_train$y)
res_p <- ggplot(res_dat, aes(x = x, y = y, fill = residual)) + 
  geom_tile() + 
  scale_fill_gradientn(colors = viridis::viridis(6)) +
  theme_minimal() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())
res_p

# Plot the masked image #
aod_mask <- ggplot() + coord_fixed() + # Ensure the aspect ratio remains fixed
  geom_tile(data = combined_data,         # Add the image layer
            aes(x = x, y = y, fill = log(AOD))) +
  geom_tile(data = mask_df,        # Overlay the mask layer
            aes(x = x, y = y), fill = "white", 
            alpha = ifelse(mask_df$mask, 0.9, 0)) +
  scale_fill_gradientn(colors = viridis::viridis(6)) +  
  theme_minimal() +                # Adjust this for your image color scale
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
aod_mask
# The white part in the middle are cloud coverage



## Variogram ##
library(geoR)
library(fields)
### variogram of raw data and residuals ###
coords_train <- combined_data_train[, c("x", "y")]
max.dist=0.5*max(rdist(coords_train))
bins=20
vario.AOD <- variog(coords=coords_train, data=combined_data_train$logAOD, 
                    uvec=(seq(0.5, max.dist, length=bins)))
plot(vario.AOD)
max.dist=0.15*max(rdist(coords_train)) #0.1
bins=45 # 30
vario.AOD.resid <- variog(coords=coords_train, 
                          data=AOD_residual2, 
                          uvec=(seq(0.5, max.dist, length=bins)))
plot(vario.AOD.resid)

vario.AOD.resid2 <- variog(coords=coords_train, 
                           data=AOD_residual2, 
                           uvec=(seq(0.5, max.dist, length=bins)),
                           estimator.type="modulus")
plot(vario.AOD.resid2)
vfit_wls=variofit(vario.AOD.resid2, 
                  ini.cov.pars=c(0.02, 10), 
                  nugget=0.001, 
                  fix.nugget=FALSE, fix.kappa = FALSE,
                  cov.model='matern', 
                  weights='cressie')
vfit_wls
plot(vario.AOD.resid2)
lines(vfit_wls, col="purple", lwd=1.5)

## decide the grids for hyper-parameters ##
library(telefit)
library(GPBayes)
#' Based on empirical semivariogram, we select range to be 1 (almost independent) 
#' to 20 meters
#' effective range c(1, 5, 10, 20)
#' delatsq: 0.1, 0.5, 1, 2
#' 1. kappa = 0.5, effective range 1 -20. phi = 3-60

source("utils2.R") # utils2.R is the testing code #
eff_range <- c(1, 5, 10, 20)
nu_grid <- c(0.5, 1.0, 1.5)
phi_ls <- decay_est(eff_range, nu_grid)
phi_nu_ls <- cbind(c(phi_ls), rep(nu_grid, each = length(eff_range))) # put all phi and nu candidate value here
colnames(phi_nu_ls) = c("phi", "nu")
deltasq_grid <- c(0.01, 0.1, 0.5)


p = 64+2
priors <- list(mu_beta = rep(0, p),
               inv_V_beta = 1/4 * diag(p),
               a_sigma = 2,
               b_sigma = 0.05)
K_fold = 3
seed = 123

run_tag <- FALSE
if(run_tag){
CV_fit_LSE <- sp_stacking_K_fold2(
  X = as.matrix(combined_data_train[, c("x", "y", paste0("V", 1:64))]), 
  y = combined_data_train$logAOD, 
  coords = coords_train,
  deltasq_grid = deltasq_grid, phi_nu_ls = phi_nu_ls, priors = priors, 
  K_fold = K_fold, seed = seed, label = "LSE")
# let's test the R^2 and coverage 
CV_fit_LSE$time #4.3 hours
save(CV_fit_LSE, file = "./RDA/result/CV_fit_LSE.RData")
}
load("./RDA/result/CV_fit_LSE.RData")
weights_LSE <- CV_fit_LSE$wts

if(run_tag){
CV_fit_LP <- sp_stacking_K_fold2(
  X = as.matrix(combined_data_train[, c("x", "y", paste0("V", 1:64))]), 
  y = combined_data_train$logAOD, 
  coords = coords_train,
  deltasq_grid = deltasq_grid, phi_nu_ls = phi_nu_ls, 
  priors = priors, K_fold = K_fold,
  seed = seed, label = "LP", MC = FALSE)
save(CV_fit_LP, file = "./RDA/result/CV_fit_LP.RData")
}
load("./RDA/result/CV_fit_LP.RData")

weights_M_LP[, r] <- CV_fit_LP$wts
run_time2["Stack_LP_E", r] <- CV_fit_LP$time[3]

## stacking mean squared prediction error ##
N_ho <- length(combined_data_test$x)
N <- length(combined_data_train$x)
y_pred_grid <- matrix(0, nrow = N_ho, ncol = nrow(CV_fit_LSE$grid_all))
w_expect_grid <- matrix(0, nrow = N+N_ho, ncol = nrow(CV_fit_LSE$grid_all))

# TO DO: generate posterior samples and compute the 95% CI
# use function Conj_pos_sam


# obtain point estimator
t <- proc.time()
for (i in 1:nrow(CV_fit_LSE$grid_all)){
  if((CV_fit_LSE$wts[i]>0.00001)){ # | (CV_fit_LP$wts[i]>0)){
    pred_grid <- 
      Conj_predict(X.mod = as.matrix(
        combined_data_train[, c("x", "y", paste0("V", 1:64))]), 
        y.mod = combined_data_train$logAOD,
        coords.mod = coords_train,
        deltasq_pick = CV_fit_LSE$grid_all$deltasq[i],
        phi_pick = CV_fit_LSE$grid_all$phi[i], 
        nu_pick = CV_fit_LSE$grid_all$nu[i], priors,
        X.ho = as.matrix(
          combined_data_test[, c("x", "y", paste0("V", 1:64))]), 
        coords.ho = as.matrix(
          combined_data_test[, c("x", "y")]))
    y_pred_grid[, i] <- pred_grid$y_expect
    w_expect_grid[, i] <- pred_grid$w_expect
  }
}
proc.time() - t
y_pred_stack_LSE = y_pred_grid %*% CV_fit_LSE$wts

# test R^2
1 - (exp(y_pred_stack_LSE) - combined_data_test$AOD)^2/ 
  (combined_data_test$AOD - mean(combined_data_test$AOD))^2

#10*3*12*time







