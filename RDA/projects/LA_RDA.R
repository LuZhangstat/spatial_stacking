rm(list = ls())
library(reticulate)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)

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

## Some simple data analysis ##
# 1. linear regression:
lm_fit <- lm(AOD~x+y+., data = combined_data)
summary(lm_fit)   # Adjusted R-squared = 0.7326
lm_fit2 <- lm(log(AOD)~x+y+., data = combined_data)
summary(lm_fit2)   # R-squared = 0.737
# check residual
AOD_residual1 <- residuals(lm_fit); hist(AOD_residual1)
qqnorm(scale(AOD_residual1)); abline(a=0,b=1)
AOD_residual2 <- residuals(lm_fit2); hist(AOD_residual2)
qqnorm(scale(AOD_residual2)); abline(a=0,b=1)
#' Obviously the residual of log AOD behaves similar to Gaussian
#' We use Log_transform here in the model fitting.
combined_data$logAOD <- log(combined_data$AOD)

# 2. Variogram



## Visualization ##
# The map of the raw AOD data in central LA ##
aod_p <- ggplot(combined_data, aes(x = x, y = y, fill = logAOD)) + 
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

# Plot the masked image #
# Convert the mask to a data frame
mask_df <- data.frame(mask = c(LA_testmask), x = x_ind, y = y_ind)
aod_mask <- ggplot() + coord_fixed() + # Ensure the aspect ratio remains fixed
  geom_tile(data = combined_data,         # Add the image layer
            aes(x = x, y = y, fill = logAOD)) +
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









