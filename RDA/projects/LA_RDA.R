rm(list = ls())
#library(reticulate)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)

# load raw data
# np <- import("numpy")
# LA_AOD <- np$load("./RDA/LAdemoData/LA_AOD_17000-997-zeroNotValid.npy")
# LA_testmask <- np$load("./RDA/LAdemoData/LA_testmask_4146.npy")
# LA_trainmask <- np$load("./RDA/LAdemoData/LA_trainmask_11857.npy")
# LA_xfea3 <- np$load("./RDA/LAdemoData/LA_xfea3_17000-997.npy")
# LA_xfea64 <- np$load("./RDA/LAdemoData/LA_xfea64_17000-997.npy")
# save(LA_AOD, LA_testmask, LA_trainmask, LA_xfea3, LA_xfea64, 
#      file = "./RDA/LAdemoData/RDALA.RData")

load("./RDA/LAdemoData/RDALA.RData")

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
sum(df_aod$AOD == 0)      # there are 997 pixels having AOD values equal to 0. 
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
# the original 3 features #
covariates_fea3_matrix <- array(data = LA_xfea3, 
                                dim = c(dim(LA_xfea3)[1] * dim(LA_xfea3)[2], 
                                        dim(LA_xfea3)[3]))
fea3_df <- as.data.frame(covariates_fea3_matrix)
fea3_df$x <- x_ind
fea3_df$y <- y_ind
names(fea3_df)[1:3] = c("EVI", "IS", "RND")

# Combine with df_aod
combined_data <- df_aod %>%
  left_join(covariates_df, by = c("x", "y"))

combined_data <- combined_data %>%
  left_join(fea3_df, by = c("x", "y"))

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
N_ho <- length(combined_data_test$x)
N <- length(combined_data_train$x)
# training: 11,857; testing: 4,146

## Some simple data analysis ##
pick_var <- c("AOD", "x", "y", paste0("V", 1:64))
# linear regression:
lm_fit <- lm(AOD~x+y+., data = combined_data_train[, pick_var])
summary(lm_fit)   # Adjusted R-squared = 0.7387
lm_fit2 <- lm(log(AOD)~x+y+., data = combined_data_train[, pick_var])
summary(lm_fit2)   # Adjusted R-squared = 0.7336
# check residual
AOD_residual1 <- residuals(lm_fit); hist(AOD_residual1)
qqnorm(scale(AOD_residual1)); abline(a=0,b=1)
AOD_residual2 <- residuals(lm_fit2); hist(AOD_residual2)
qqnorm(scale(AOD_residual2)); abline(a=0,b=1)
#' The residual of log AOD behaves similar to Gaussian
#' We use Log_transform here in the model fitting.
combined_data_train$logAOD <- log(combined_data_train$AOD)
pred_lm <- predict(lm_fit2, combined_data_test[, c("x", "y", paste0("V", 1:64))])
1 - sum((exp(c(pred_lm)) - combined_data_test$AOD)^2)/ 
  sum((combined_data_test$AOD - mean(combined_data_test$AOD))^2)
#0.6877944

# consider the original 3 features EVI IS RND #
pick_var2 <- c("AOD", "x", "y", "EVI", "IS", "RND")
lm_fit_2 <- lm(AOD~x+y+., data = combined_data_train[, pick_var2])
summary(lm_fit_2)   # Adjusted R-squared = 0.5
lm_fit2_2 <- lm(log(AOD)~x+y+., data = combined_data_train[, pick_var2])
summary(lm_fit2_2)   # Adjusted R-squared = 0.478
# check residual
AOD_residual1_2 <- residuals(lm_fit_2); hist(AOD_residual1_2)
qqnorm(scale(AOD_residual1_2)); abline(a=0,b=1)
AOD_residual2_2 <- residuals(lm_fit2_2); hist(AOD_residual2_2)
qqnorm(scale(AOD_residual2_2)); abline(a=0,b=1)
pred_lm_2 <- predict(lm_fit2_2, 
                     combined_data_test[, c("x", "y", "EVI", "IS", "RND")])
1 - sum((exp(c(pred_lm_2)) - combined_data_test$AOD)^2)/ 
  sum((combined_data_test$AOD - mean(combined_data_test$AOD))^2)
#0.45



## Visualization ##
library("ggplot2")
library("viridis")
library("gridExtra")
library(gtable)

min_log_aod <- min(log(combined_data$AOD), na.rm = TRUE)
max_log_aod <- max(log(combined_data$AOD), na.rm = TRUE)

# Determine the range 
x_range <- range(combined_data$x, na.rm = TRUE)
y_range <- range(combined_data$y, na.rm = TRUE)
# Calculate aspect ratio (if necessary)
aspect_ratio <- diff(x_range) / diff(y_range)

# Figure 7: AOD and Predictors Visualization
# The map of the raw AOD data in central LA ##
aod_p <- ggplot() + 
  coord_fixed(#ratio = aspect_ratio, 
    xlim = x_range, ylim = y_range) +
  geom_tile(data = combined_data, 
            aes(x = x, y = y, fill = log(AOD))) + 
  scale_fill_gradientn(colors = viridis::viridis(6), 
                       limits = c(min_log_aod, max_log_aod),
                       name = "Log AOD") +
  theme_minimal() + labs(title = "Raw AOD") +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))
aod_p
aod_p_legend <- gtable_filter(ggplot_gtable(ggplot_build(aod_p)), "guide-box")
aod_p <- aod_p + guides(fill = "none")
# Plot the masked image #
aod_mask <- ggplot() + 
  coord_fixed(#ratio = aspect_ratio, 
    xlim = x_range, ylim = y_range) +
  geom_tile(data = combined_data,         # Add the image layer
            aes(x = x, y = y, fill = log(AOD))) +
  geom_tile(data = mask_df,        # Overlay the mask layer
            aes(x = x, y = y), fill = "white", 
            alpha = ifelse(mask_df$mask, 0.95, 0)) +
  scale_fill_gradientn(colors = viridis::viridis(6),
                       limits = c(min_log_aod, max_log_aod)) +  
  theme_minimal() +                # Adjust this for your image color scale
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + labs(title = "Training") +
  guides(fill = "none")#,
#legend.position = "none")
aod_mask
# The white part in the middle are cloud coverage

aod_test <- ggplot() + 
  coord_fixed(#ratio = aspect_ratio, 
    xlim = x_range, ylim = y_range) +
  geom_tile(data = combined_data,         # Add the image layer
            aes(x = x, y = y, fill = log(AOD))) +
  geom_tile(data = mask_df,        # Overlay the mask layer
            aes(x = x, y = y), fill = "white", 
            alpha = ifelse(mask_df$mask, 0, 0.95)) +
  scale_fill_gradientn(colors = viridis::viridis(6),
                       limits = c(min_log_aod, max_log_aod)) +  
  theme_minimal() +                # Adjust this for your image color scale
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + labs(title = "Testing") +
  guides(fill = "none")#,
#legend.position = "none")
aod_test


combined_plot <- grid.arrange(aod_p, aod_mask, aod_test,
                              aod_p_legend, ncol = 4, widths = c(1, 1, 1, 0.4))

#ggsave("./RDA/pics/AOD_plots.png", combined_plot, width = 9, height = 3)

# check the plot of covariates
V_p <- ggplot(combined_data, aes(x = x, y = y, fill = V1)) + 
  geom_tile() + 
  scale_fill_gradientn(colors = viridis::viridis(6)) +
  theme_minimal() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + labs(title = "V1")
V_p

EVI_p <- ggplot(fea3_df, aes(x = x, y = y, fill = EVI)) + 
  geom_tile() + 
  coord_fixed(xlim = x_range, ylim = y_range) +
  scale_fill_gradientn(colors = viridis::viridis(6)) +
  theme_minimal() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + labs(title = "EVI")
EVI_p

IS_p <- ggplot(fea3_df, aes(x = x, y = y, fill = IS)) + 
  geom_tile() + 
  coord_fixed( xlim = x_range, ylim = y_range) +
  scale_fill_gradientn(colors = viridis::viridis(6)) +
  theme_minimal() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Impervious Surface %")
IS_p

RND_p <- ggplot(fea3_df, aes(x = x, y = y, fill = RND)) + 
  geom_tile() + 
  coord_fixed(xlim = x_range, ylim = y_range) +
  scale_fill_gradientn(colors = viridis::viridis(6)) +
  theme_minimal() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Road Network Density")
RND_p

fea3_plot <- grid.arrange(EVI_p, IS_p, RND_p, ncol = 3)

#ggsave("./RDA/pics/fea3_plots.png", fea3_plot, width = 10, height = 3.5)


# check the residual #
res_dat <- data.frame(residual = AOD_residual2_2, x = combined_data_train$x,
                      y = combined_data_train$y)
res_p <- ggplot(res_dat, aes(x = x, y = y, fill = residual)) + 
  geom_tile() + 
  scale_fill_gradientn(colors = viridis::viridis(6)) +
  theme_minimal() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())
res_p




## Variogram: Figure S21 ##
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
bins=60 # 30
vario.AOD.resid <- variog(coords=coords_train, 
                          data=AOD_residual2, 
                          uvec=(seq(0.5, max.dist, length=bins)))
plot(vario.AOD.resid)

vario.AOD.resid2 <- variog(coords=coords_train, 
                           data=AOD_residual2_2, #AOD_residual2
                           uvec=(seq(0.5, max.dist, length=bins)),
                           estimator.type="modulus")
vfit_wls=variofit(vario.AOD.resid2, 
                  ini.cov.pars=c(0.02, 10), 
                  nugget=0.001, 
                  fix.nugget=FALSE, fix.kappa = FALSE, 
                  cov.model='matern', 
                  weights='cressie')
vfit_wls
fitted_nugget <- vfit_wls$nugget
fitted_sill <- vfit_wls$nugget + vfit_wls$cov.pars[1]
#png("./RDA/pics/res_variogram_3fea.png", width=400, height=300)  # Start PNG device and specify file name
plot(vario.AOD.resid2, xlab="Distance", ylab="Semivariance", 
     ylim = c(-0.001, 0.12))
lines(vfit_wls, col="purple", lwd=2)  # Adjust line color and width
abline(h=fitted_nugget, col="blue", lty=2)  # Nugget line
abline(h=fitted_sill, col="red", lty=2)     # Sill line
# Adding text labels
text(x=max.dist/2, y=fitted_sill+0.002, 
     labels=paste("estimated sill:", round(fitted_sill, 3)), pos=3, col="red",
     cex = 1)
text(x=max.dist/2, y=fitted_nugget+0.002, 
     labels=paste("estimated nugget:", round(fitted_nugget, 3)), pos=3, 
     col="blue", cex = 1)
#dev.off()

## decide the grids for hyper-parameters ##

# method 1.
range(c(1 / Matern.cor.to.range(5, 0.5, cor.target=.05),
        1 / Matern.cor.to.range(25, 0.5, cor.target=.05),
        1 / Matern.cor.to.range(25, 1.75, cor.target=.05),
        1 / Matern.cor.to.range(5, 1.75, cor.target=.05)))

#~0.2 to ~0.8, 0.12 to 1.01
source("utils.R")
phi_grid = c(0.1, 0.4, 0.7, 1.0)
nu_grid = c(0.5, 1.0, 1.5, 1.75)
deltasq_grid <- pick_deltasq(E_sigmasq = 0.088, E_tausq = 0.01, b = 0.088,
                             p_ls = c(0.05, 0.35, 0.65, 0.95))
deltasq_grid

p = 3+2+1
priors <- list(mu_beta = rep(0, p),
               inv_V_beta = 1/4 * diag(p),
               a_sigma = 2,
               b_sigma = 0.05)
K_fold = 10
seed = 123

#pars1 <- paste0("V", 1:64)
pars2 <- c("EVI", "IS", "RND")
run_tag <- FALSE
# fit stacking of means #
if(run_tag){
  CV_fit_LSE <- sp_stacking_K_fold(
    X = cbind(1, as.matrix(
      combined_data_train[, c("x", "y", "EVI", "IS", "RND")])), 
    y = combined_data_train$logAOD, 
    coords = coords_train,
    deltasq_grid = deltasq_grid, phi_grid = phi_grid,
    nu_grid = nu_grid, priors = priors, 
    K_fold = K_fold, seed = seed, label = "LSE")
  CV_fit_LSE$time #9.6 hour
  save(CV_fit_LSE, file = "./RDA/result/CV_fit_LSE_fea3.RData")
  
  cbind(CV_fit_LSE$grid_all[CV_fit_LSE$wts>0.00001, ],
        CV_fit_LSE$wts[CV_fit_LSE$wts>0.00001])
  
  pos_sam_LSE <- 
    stacking_pos_sample(
      Stack_fit = CV_fit_LSE, L1 = 600, L2 = 1800, 
      X.mod = cbind(1, as.matrix(
        combined_data_train[, c("x", "y", "EVI", "IS", "RND")])),
      y.mod = combined_data_train$logAOD, 
      coords.mod = coords_train, priors = priors,
      X.ho = cbind(1, as.matrix(
        combined_data_test[, c("x", "y", "EVI", "IS", "RND")])), 
      coords.ho = as.matrix(
        combined_data_test[, c("x", "y")]), seed = 123)
  
  save(pos_sam_LSE, file = "./RDA/result/pos_sam_LSE_fea3.RData")
}else{
  load("./RDA_carc/result/CV_fit_LSE_fea3.RData")
  load("./RDA_carc/result/pos_sam_LSE_fea3.RData")
}
cbind(CV_fit_LSE$grid_all[CV_fit_LSE$wts>0.00001, ], 
      CV_fit_LSE$wts[CV_fit_LSE$wts>0.00001])

# fit stacking of predictive densities #
if(run_tag){
  CV_fit_LP <- sp_stacking_K_fold(
    X = cbind(1, as.matrix(
      combined_data_train[, c("x", "y", "EVI", "IS", "RND")])), 
    y = combined_data_train$logAOD, 
    coords = coords_train,
    deltasq_grid = deltasq_grid, phi_grid = phi_grid,
    nu_grid = nu_grid, priors = priors, K_fold = K_fold,
    seed = seed, label = "LP", MC = FALSE)
  CV_fit_LP$time # 8.45 hours
  save(CV_fit_LP, file = "./RDA/result/CV_fit_LP_fea3.RData")
  cbind(CV_fit_LP$grid_all[CV_fit_LP$wts>0.00001, ], 
        CV_fit_LP$wts[CV_fit_LP$wts>0.00001])
  
  pos_sam_LP <- 
    stacking_pos_sample(
      Stack_fit = CV_fit_LP, L1 = 1000, L2 = 3000, 
      X.mod = cbind(1, as.matrix(
        combined_data_train[, c("x", "y", "EVI", "IS", "RND")])),
      y.mod = combined_data_train$logAOD, 
      coords.mod = coords_train, priors = priors,
      X.ho = cbind(1, as.matrix(
        combined_data_test[, c("x", "y", "EVI", "IS", "RND")])), 
      coords.ho = as.matrix(
        combined_data_test[, c("x", "y")]), seed = 123)
  
  save(pos_sam_LP, file = "./RDA/result/pos_sam_LP_fea3.RData")
}else{
  load("./RDA_carc/result/CV_fit_LP_fea3.RData")
  load("./RDA_carc/result/pos_sam_LP_fea3.RData")
}
cbind(CV_fit_LP$grid_all[CV_fit_LP$wts>0.00001, ], 
      CV_fit_LP$wts[CV_fit_LP$wts>0.00001])


## For LSE ##
cat("LSE R2:")
1 - sum((rowMeans(pos_sam_LSE$pred_y_U_stack_sam) - log(combined_data_test$AOD))^2)/ 
  sum((log(combined_data_test$AOD) - mean(log(combined_data_test$AOD)))^2)
#0.8455106

cor(combined_data_test$AOD, rowMeans(exp(pos_sam_LSE$pred_y_U_stack_sam)))
# 0.9269146

# 95% CI. coverage
y_U_CI_stack_LSE <- 
  apply(pos_sam_LSE$pred_y_U_stack_sam, 1, 
        function(x){exp(quantile(x, probs = c(0.025, 0.975)))})
cat("LSE 95% CI coverage: ")
sum((y_U_CI_stack_LSE[1, ] < combined_data_test$AOD) & 
      (y_U_CI_stack_LSE[2, ] > combined_data_test$AOD))/4146
# 0.7590449

# 99% CI. coverage
y_U_CI_stack_LSE <- 
  apply(pos_sam_LSE$pred_y_U_stack_sam, 1, 
        function(x){exp(quantile(x, probs = c(0.005, 0.9995)))})
cat("LSE 99% CI coverage: ")
sum((y_U_CI_stack_LSE[1, ] < combined_data_test$AOD) & 
      (y_U_CI_stack_LSE[2, ] > combined_data_test$AOD))/4146
# 0.8699952

# root-mean-squared error (RMSE)
sqrt(mean((combined_data_test$AOD - 
             rowMeans(exp(pos_sam_LSE$pred_y_U_stack_sam)))^2))
#0.0100499

# mean absolute error (MAE)
mean(abs(combined_data_test$AOD - 
           rowMeans(exp(pos_sam_LSE$pred_y_U_stack_sam))))
# 0.007049124

## For LP ##
cat("LP R2:")
1 - sum((rowMeans(pos_sam_LP$pred_y_U_stack_sam) - log(combined_data_test$AOD))^2)/ 
  sum((log(combined_data_test$AOD) - mean(log(combined_data_test$AOD)))^2)
#0.8419137

# correlaton
cor(combined_data_test$AOD, rowMeans(exp(pos_sam_LP$pred_y_U_stack_sam)))
#0.9261988

# 95% CI. coverage
y_U_CI_stack_LP <- 
  apply(pos_sam_LP$pred_y_U_stack_sam, 1, 
        function(x){exp(quantile(x, probs = c(0.025, 0.975)))})
cat("LP 95% CI coverage: ")
sum((y_U_CI_stack_LP[1, ] < combined_data_test$AOD) & 
      (y_U_CI_stack_LP[2, ] > combined_data_test$AOD))/4146
# 0.8212735

# 95% CI. coverage
y_U_CI_stack_LP <- 
  apply(pos_sam_LP$pred_y_U_stack_sam, 1, 
        function(x){exp(quantile(x, probs = c(0.005, 0.9995)))})
cat("LP 99% CI coverage: ")
sum((y_U_CI_stack_LP[1, ] < combined_data_test$AOD) & 
      (y_U_CI_stack_LP[2, ] > combined_data_test$AOD))/4146
# 0.9394597

# root-mean-squared error (RMSE)
sqrt(mean((combined_data_test$AOD - 
             rowMeans(exp(pos_sam_LP$pred_y_U_stack_sam)))^2))
#0.0101573


# mean absolute error (MAE)
mean(abs(combined_data_test$AOD - 
           rowMeans(exp(pos_sam_LP$pred_y_U_stack_sam))))
# 0.007143969


# fit Bayesian linear regression model #
library(brms)
library(bayesplot)
options(mc.cores = parallel::detectCores() - 2)
prior1 <- prior(normal(0,2), class = b) +
  prior(cauchy(0, 0.5), class = sigma)
fit_lm <- brm(logAOD ~ x + y + EVI + IS + RND,
              data = combined_data_train, family = gaussian(), 
              prior = prior1, seed = 123)
summary(fit_lm)
predicted_samples_lm <- 
  posterior_predict(fit_lm, newdata = combined_data_test)

# correlaton
cor(combined_data_test$AOD, colMeans(exp(predicted_samples_lm)))
#0.6932424

# 95% CI. coverage
y_U_CI_lm <- 
  apply(predicted_samples_lm, 2, 
        function(x){exp(quantile(x, probs = c(0.025, 0.975)))})
cat("BLM 95% CI coverage: ")
sum((y_U_CI_lm[1, ] < combined_data_test$AOD) & 
      (y_U_CI_lm[2, ] > combined_data_test$AOD))/4146
# 0.9493488

# 95% CI. coverage
y_U_CI_lm <- 
  apply(predicted_samples_lm, 2, 
        function(x){exp(quantile(x, probs = c(0.005, 0.9995)))})
cat("LP 99% CI coverage: ")
sum((y_U_CI_lm[1, ] < combined_data_test$AOD) & 
      (y_U_CI_lm[2, ] > combined_data_test$AOD))/4146
# 0.9937289

# root-mean-squared error (RMSE)
sqrt(mean((combined_data_test$AOD - 
             colMeans(exp(predicted_samples_lm)))^2))
#0.01941467

# mean absolute error (MAE)
mean(abs(combined_data_test$AOD - 
           exp(colMeans(predicted_samples_lm))))
# 0.01467183

#################################
## compare the predicted plots ##
#################################
# Figure 8
# load results from H2O models #
H2Opred <- read.csv("./RDA/result/H2Opred.csv")
names(H2Opred) <- c("DL", "RF", "GB", "Ensemble")

pick_show <- which(combined_data_test$x>55 & combined_data_test$x<140 & 
        combined_data_test$y>20 & combined_data_test$y<80)
combined_data_test$logAOD = log(combined_data_test$AOD)
x_range <- c(55, 140); y_range <- c(20, 80)

aod_test <- ggplot() + 
  coord_fixed(#ratio = aspect_ratio, 
    xlim = x_range, ylim = y_range) +
  geom_tile(data = combined_data_test[pick_show, ],         # Add the image layer
            aes(x = x, y = y, fill = log(AOD))) +
  scale_fill_gradientn(colors = viridis::viridis(6),
                       limits = c(min_log_aod, max_log_aod),
                       name = "Log AOD") +  
  theme_minimal() +                # Adjust this for your image color scale
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Testing (raw data)") 
aod_p_legend <- gtable_filter(ggplot_gtable(ggplot_build(aod_test)), "guide-box")
aod_test <- aod_test + guides(fill = "none")
aod_test

LSE_data <- data.frame(logAOD = rowMeans(pos_sam_LSE$pred_y_U_stack_sam),
                       x = combined_data_test$x,
                       y = combined_data_test$y)
aod_LSE <- ggplot() + 
  coord_fixed(#ratio = aspect_ratio, 
    xlim = x_range, ylim = y_range) +
  geom_tile(data = LSE_data[pick_show,],         # Add the image layer
            aes(x = x, y = y, fill = logAOD)) +
  geom_tile(data = mask_df,        # Overlay the mask layer
            aes(x = x, y = y), fill = "white", 
            alpha = ifelse(mask_df$mask, 0, 0.95)) +
  scale_fill_gradientn(colors = viridis::viridis(6),
                       limits = c(min_log_aod, max_log_aod)) +  
  theme_minimal() +                # Adjust this for your image color scale
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Stacking of means") +
  guides(fill = "none")#,
aod_LSE

LP_data <- data.frame(logAOD = rowMeans(pos_sam_LP$pred_y_U_stack_sam),
                       x = combined_data_test$x,
                       y = combined_data_test$y)
aod_LP <- ggplot() + 
  coord_fixed(#ratio = aspect_ratio, 
    xlim = x_range, ylim = y_range) +
  geom_tile(data = LP_data[pick_show, ],         # Add the image layer
            aes(x = x, y = y, fill = logAOD)) +
  geom_tile(data = mask_df,        # Overlay the mask layer
            aes(x = x, y = y), fill = "white", 
            alpha = ifelse(mask_df$mask, 0, 0.95)) +
  scale_fill_gradientn(colors = viridis::viridis(6),
                       limits = c(min_log_aod, max_log_aod)) +  
  theme_minimal() +                # Adjust this for your image color scale
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Stacking of pds") +
  guides(fill = "none")#,
aod_LP

BLM_data <- data.frame(logAOD = colMeans(predicted_samples_lm),
                      x = combined_data_test$x,
                      y = combined_data_test$y)
aod_BLM <- ggplot() + 
  coord_fixed(#ratio = aspect_ratio, 
    xlim = x_range, ylim = y_range) +
  geom_tile(data = BLM_data[pick_show, ],         # Add the image layer
            aes(x = x, y = y, fill = logAOD)) +
  geom_tile(data = mask_df,        # Overlay the mask layer
            aes(x = x, y = y), fill = "white", 
            alpha = ifelse(mask_df$mask, 0, 0.95)) +
  scale_fill_gradientn(colors = viridis::viridis(6),
                       limits = c(min_log_aod, max_log_aod)) +  
  theme_minimal() +                # Adjust this for your image color scale
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Bayesian linear model") +
  guides(fill = "none")#,
aod_BLM

RF_data <- data.frame(logAOD = log(H2Opred$RF),
                       x = combined_data_test$x,
                       y = combined_data_test$y)
aod_RF <- ggplot() + 
  coord_fixed(#ratio = aspect_ratio, 
    xlim = x_range, ylim = y_range) +
  geom_tile(data = RF_data[pick_show, ],         # Add the image layer
            aes(x = x, y = y, fill = logAOD)) +
  geom_tile(data = mask_df,        # Overlay the mask layer
            aes(x = x, y = y), fill = "white", 
            alpha = ifelse(mask_df$mask, 0, 0.95)) +
  scale_fill_gradientn(colors = viridis::viridis(6),
                       limits = c(min_log_aod, max_log_aod)) +  
  theme_minimal() +                # Adjust this for your image color scale
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Random Forest") +
  guides(fill = "none")#,
aod_RF


combined_plot_pre <- grid.arrange(aod_test, aod_RF,  aod_p_legend, 
                                  aod_LSE, aod_LP, aod_p_legend, 
                              ncol = 3, nrow = 2, widths = c(1, 1, 0.3))

ggsave("./RDA/pics/AOD_pred.png", combined_plot_pre, width = 8, height = 5.5)

## Figure S22: 95% CI coverage ##
library(ggplot2)
library(dplyr)
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
lp_data$x <- log(combined_data_test$AOD)
cover_lp <- round(sum(lp_data$lower<lp_data$x & lp_data$upper>lp_data$x) / 
                    length(lp_data$x) *100, 1) 
label_lp <- paste0('stacking of pds: ', cover_lp,"% coverage")
lp_data$source <- label_lp
lse_data <- calculate_means_and_cis(pos_sam_LSE$pred_y_U_stack_sam)
lse_data$x <- log(combined_data_test$AOD)
cover_lse <- round(sum(lse_data$lower<lse_data$x & lse_data$upper>lse_data$x) / 
                     length(lse_data$x) *100, 1) 
label_lse <- paste0('stacking of means: ', cover_lse, "% coverage")
lse_data$source <- label_lse
Blm_data <- calculate_means_and_cis(t(predicted_samples_lm))
Blm_data$x <- log(combined_data_test$AOD)
cover_Blm <- round(sum(Blm_data$lower<Blm_data$x & 
                          Blm_data$upper>Blm_data$x) / 
                      length(Blm_data$x) *100, 1) 
label_Blm <- paste0('Bayesian linear: ', cover_Blm, "% coverage")
Blm_data$source <- label_Blm


# Combine data
combined_data <- rbind(lse_data, lp_data, Blm_data)
combined_data$source <- factor(combined_data$source, 
                               levels = c(label_lse, label_lp, label_Blm))

# Determine the common range for x and y if not already known
x_range <- range(combined_data$x, na.rm = TRUE)
y_range <- range(combined_data$mean, combined_data$lower, 
                 combined_data$upper, na.rm = TRUE)
# Plot
pts_y <- ggplot(combined_data, aes(x = x, y = mean)) +
  geom_point(size = 0.2, color = "blue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.005, alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1, color = "green") +
  facet_wrap(~ source) +  # Use 'free_x' to unify y range
  scale_x_continuous(limits = x_range) +
  scale_y_continuous(limits = y_range) +
  theme_bw() +
  labs(x = "log AOD", y = "mean and 95%CI")
pts_y
ggsave(paste0("./RDA/pics/y_U_95CI.png"),
       plot = pts_y,
       width = 8, height = 3, units = "in", dpi = 600)

