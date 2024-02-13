rm(list = ls())
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)

load("./RDA/LAdemoData/RDALA.RData")

# transform the data from matrix into a dataframe
n_row <- nrow(LA_AOD)
n_col <- ncol(LA_AOD)
x_ind <- rep(1:n_col, each = n_row)
y_ind <- rep(n_row:1, time = n_col)    # y from high to low

df_aod <- data.frame(AOD = c(LA_AOD), x = x_ind, y = y_ind)
df_aod %>% glimpse()
df_aod <- df_aod %>% filter(AOD > 0)   # Remove 0 from the data
df_aod %>% glimpse()
# our data has 16,003 pixles

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

# try h2o for machine learning algorithm
#install.packages("h2o")
library(h2o)
h2o.init()

#h2o.init(nthreads = -1)

names(combined_data)[2:3] <- c("x_ind", "y_ind")
names(combined_data_train)[2:3] <- c("x_ind", "y_ind")
names(combined_data_test)[2:3] <- c("x_ind", "y_ind")

y <- "AOD"  
x <- setdiff(names(combined_data), y)
train <- as.h2o(combined_data_train)
test <- as.h2o(combined_data_test)


## 1. Deep Learning ##
m1 <- h2o.deeplearning(x, y, train, nfolds=10, seed = 1,
                       fold_assignment = "Modulo", 
                       keep_cross_validation_predictions = T,
                       reproducible = TRUE)  
p1 <- h2o.predict(m1, test)  
#as.data.frame(p)
pred_perf1 <- h2o.performance(m1, test)
pred_perf1@metrics$mae # 0.01266137
pred_perf1@metrics$RMSE # 0.01658932
h2o.cor(p1$predict, test$AOD)
# 0.792391

## 2. Distributed Random Forest ##

g <- h2o.grid("randomForest",
              hyper_params = list(
                ntrees = c(50, 100, 120),
                max_depth = c(40, 60),
                min_rows = c(1, 2)
              ), seed = 1,
              x = x, y = y, training_frame = train, nfolds = 10)
g

m2 <- h2o.randomForest(x, y, train, nfolds = 10, model_id = "RF_defaults",
                       fold_assignment = "Modulo", seed = 1,
                       keep_cross_validation_predictions = T,
                       ntrees = 120, max_depth = 60, min_rows = 1)
p2 <- h2o.predict(m2, test) 
pred_perf2 <- h2o.performance(m2, test)
pred_perf2@metrics$mae # 0.01153211
pred_perf2@metrics$RMSE # 0.01526932
h2o.cor(p2$predict, test$AOD)
# 0.8208408

## 3. Gradient Boosting ##
m3 <- h2o.gbm(x, y, train, nfolds = 10, model_id = "GBM_defaults",
              fold_assignment = "Modulo", seed = 1,
              keep_cross_validation_predictions = T)
p3 <- h2o.predict(m3, test) 
pred_perf3 <- h2o.performance(m3, test)
pred_perf3@metrics$mae #0.01221575
pred_perf3@metrics$RMSE #0.01609126
h2o.cor(p3$predict, test$AOD)
# 0.7969876

## 4. Ensemble model ##
library(h2oEnsemble)
RFd <- h2o.randomForest(x, y, train, model_id="RF_defaults", nfolds=10, 
                        seed = 1, fold_assignment = "Modulo", 
                        keep_cross_validation_predictions = T,
                        ntrees = 120, max_depth = 60, min_rows = 1)
GBMd <- h2o.gbm(x, y, train, model_id="GBM_defaults", nfolds=10, seed = 1,
                fold_assignment = "Modulo", 
                keep_cross_validation_predictions = T)
#GLMd <- h2o.glm(x, y, train, model_id="GLM_defaults", nfolds=10, seed = 1,
#                fold_assignment = "Modulo", 
#                keep_cross_validation_predictions = T)
DLd <- h2o.deeplearning(x, y, train, model_id="DL_defaults", nfolds=10, 
                        seed = 1, fold_assignment = "Modulo", 
                        keep_cross_validation_predictions = T, 
                        reproducible = TRUE)

models <- c(RFd, GBMd, DLd)

# Combine base models into a stacked ensemble
ensemble <- h2o.stackedEnsemble(x, y, training_frame = train, seed = 1,
                                base_models = c(RFd, GBMd, GLMd, DLd),
                                metalearner_algorithm = "glm") # Using GLM as metalearner for regression
pred_ensem <- h2o.predict(ensemble, newdata = test)
pred_perf_ensem <- h2o.performance(ensemble, test)
pred_perf_ensem@metrics$mae # 0.01151554
pred_perf_ensem@metrics$RMSE # 0.01532827
h2o.cor(pred_ensem$predict, test$AOD)
# 0.8176291

#h2o.shutdown()


