rm(list=ls())
library(MBA)
library(fields)
library(classInt)
library(geoR)
library(spBayes)
library(rstanarm) # fit Bayesian linear regression
data("WEF.dat")

###################################################
### WEF data
###################################################

set.seed(1)
# simple data cleaning 
indexNA <- (is.na(WEF.dat$East_m) | is.na(WEF.dat$North_m) | 
              is.na(WEF.dat$DBH_cm) | is.na(WEF.dat$Species) | 
              WEF.dat$Live_dead == "D" | WEF.dat$Species == "NF") 
WEF.dat <- droplevels(
  WEF.dat[!indexNA, c("East_m", "North_m", "DBH_cm", "Species")])

table(as.numeric(WEF.dat$Species))
levels(WEF.dat$Species)
ind=sample(1:nrow(WEF.dat), 500, replace=FALSE) # pick 500 observation to check prediction 

### holdout data to assess RMSPE ###
WEF.out=WEF.dat[ind,]   # holdout
WEF.in=WEF.dat[-ind,]   # train 
rm("WEF.dat")

### diameter at breast height for the trees
DBH <- WEF.in$DBH_cm

coords <- as.matrix(WEF.in[,c("East_m","North_m")])

col.br <- colorRampPalette(c("midnightblue", "cyan", "yellow",  "orange", "red"))
col.pal <- col.br(5)

### DBH quantile based color coding of the locations
quant <- classIntervals(DBH, n=5, style="quantile")

quant.col <- findColours(quant, col.pal)

#dev.new()
plot(coords, col=quant.col, pch=19, cex=1.5, main="", xlab="Easting (m)", ylab="Northing (m)")
legend("topleft", fill=attr(quant.col, "palette"), 
       legend=names(attr(quant.col, "table")), bty="n",cex=1.3)

### plot of interpolated surface using mba package ###
surf <- mba.surf(cbind(coords, DBH), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))

### plot of Species type ###
spnum=as.factor(WEF.in$Species)
col.pal2 <- col.br(length(unique(spnum)))

plot(coords, col=col.pal2[spnum], pch=19, cex=0.5, main="", xlab="Easting (m)", ylab="Northing (m)")
legend("topleft", fill=col.pal2, 
       legend=levels(spnum), bty="n")

### linear regression ###
lm.DBH <- lm(DBH~Species, data=WEF.in)
summary(lm.DBH)
DBH.resid <- resid(lm.DBH)

surf <- mba.surf(cbind(coords,DBH.resid), no.X=100, no.Y=100, h=5, m=2, 
                 extend=FALSE)$xyz.est
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", 
           ylab="Northing (m)", col=col.br(25))


### variogram of raw data and residuals ###
max.dist=0.5*max(rdist(coords))
bins=20

vario.DBH <- variog(coords=coords, data=DBH, uvec=(seq(5, max.dist, length=bins)))
plot(vario.DBH)

vario.DBH.resid <- variog(coords=coords, data=DBH.resid, uvec=(seq(0, max.dist, length=bins)))
plot(vario.DBH.resid)


#########################
### Linear regression ###
#########################

model_bayes_lm<- stan_glm(DBH~Species, data=WEF.in, seed=123, 
                       prior = normal(0, 2.5, autoscale = TRUE),
                       #prior_intercept = normal(50, 2.5, autoscale = TRUE),
                       prior_aux = exponential(1, autoscale = TRUE),
                       chains = 4, iter = 2000)


####################
###   stacking   ###
####################
source("utils.R")
# set the prior of stacking to be close to the Bayesian linear regression
prior_summary(model_bayes_lm)
X <- model.matrix(~Species, data = WEF.in) # design matrix 
p = ncol(X)
priors <- list(mu_beta = c(37, 0, 0, 0),
               inv_V_beta = diag(1/c(81, 677.58, 163.85, 185.60)^2),
               a_sigma = 2,
               b_sigma = 1000)

# pick candidate values of phi #
range(1/Matern.cor.to.range(0.6 * max(rdist(coords)), 0.5, cor.target=.05),  # 0.01417578
      1/Matern.cor.to.range(0.1 * max(rdist(coords)), 0.5, cor.target=.05), # 0.08505468
      1/Matern.cor.to.range(0.6 * max(rdist(coords)), 1.75, cor.target=.05), # 0.02397907
      1/Matern.cor.to.range(0.1 * max(rdist(coords)), 1.75, cor.target=.05)) # 0.1438744

# candidate values for the hyperparameters #
deltasq_grid <- c(0.1, 0.5, 1, 2)
phi_grid = seq(0.014, 0.14, length.out = 4)
nu_grid = c(0.5, 1, 1.5, 1.75)

K_fold = 10

# compute the stacking weights #
t <- proc.time()
CV_fit_LSE <- sp_stacking_K_fold(
  X = X, y = WEF.in$DBH_cm, coords = coords,
  deltasq_grid = deltasq_grid, phi_grid = phi_grid,
  nu_grid = nu_grid, priors = priors, K_fold = K_fold,
  seed = 123, label = "LSE")
weights_M_LSE <- round(CV_fit_LSE$wts, digits = 4)
proc.time() - t #150s

t <- proc.time()
CV_fit_LP <- sp_stacking_K_fold(
  X = X, y = WEF.in$DBH_cm, coords = coords,
  deltasq_grid = deltasq_grid, phi_grid = phi_grid,
  nu_grid = nu_grid, priors = priors, K_fold = K_fold,
  seed = 123, label = "LP")
weights_M_LP <- round(CV_fit_LP$wts, digits = 4)
proc.time() - t #150s

#############
### RMSPE ###
#############

# linear regression #
pred_lm=as.vector(as.matrix(trend.spatial(~Species,WEF.out))%*%lm.DBH$coefficients)
rmspe_lm=sqrt(mean((pred_lm-WEF.out$DBH_cm)^2))
rmspe_lm

# Bayesian linear regression #
ytilde <- posterior_predict(model_bayes_lm, WEF.out, draws = 1000)
print(dim(ytilde))
yhat <- apply(ytilde, 2, mean) # posterior prediction for held out locations
rmspe_Blm=sqrt(mean((yhat-WEF.out$DBH_cm)^2))


# stacking #
# design matrix and locations of the data for prediction
X_ho <- model.matrix(~Species, data = WEF.out)
coords_ho <- as.matrix(WEF.out[,c("East_m","North_m")])
N_ho = nrow(X_ho)
## stacking mean squared prediction error ##
y_pred_grid <- matrix(0, nrow = N_ho, ncol = nrow(CV_fit_LSE$grid_all))
w_expect_grid <- matrix(0, nrow = N_ho + nrow(X), 
                        ncol = nrow(CV_fit_LSE$grid_all))

t <- proc.time()
for (i in 1:nrow(CV_fit_LSE$grid_all)){
  if( (CV_fit_LSE$wts[i]>0) | (CV_fit_LP$wts[i]>0)){
    pred_grid <- Conj_predict(X.mod = X, 
                              y.mod = WEF.in$DBH_cm,
                              coords.mod = coords,
                              deltasq_pick = CV_fit_LSE$grid_all$deltasq[i],
                              phi_pick = CV_fit_LSE$grid_all$phi[i], 
                              nu_pick = CV_fit_LSE$grid_all$nu[i],
                              priors,
                              X.ho = X_ho, 
                              coords.ho = coords_ho)
    y_pred_grid[, i] <- pred_grid$y_expect
    w_expect_grid[, i] <- pred_grid$w_expect
  }
}
proc.time() - t

y_pred_stack_LSE = y_pred_grid %*% CV_fit_LSE$wts
RMSPE_stack_LSE <- sqrt(mean((y_pred_stack_LSE - WEF.out$DBH_cm)^2))
y_pred_stack_LP = y_pred_grid %*% CV_fit_LP$wts
RMSPE_stack_LP <- sqrt(mean((y_pred_stack_LP - WEF.out$DBH_cm)^2))

print(round(c(rmspe_lm, rmspe_Blm, RMSPE_stack_LSE, RMSPE_stack_LP), 
            digits = 2))
# 22.84 22.86 20.70 20.79

############
### ELPD ###
############

# Bayesian linear regression #
draws <- as.matrix(model_bayes_lm, pars = NULL, regex_pars = NULL)
M_r <- X_ho%*%t(draws[, 1:4])  - WEF.out$DBH_cm
M_r <- M_r %*% Diagonal(n = nrow(draws), 1 / draws[, 5])
M_r <- -0.5 * (M_r^2 + log(2 * pi) +
                 tcrossprod(rep(1, N_ho), 2*log(draws[, 5])))
lp_expect <- log(rowMeans(exp(M_r)))
ELPD_Blm <- mean(lp_expect)

# stacking #
# stacking Expected log pointwise predictive density 
lp_pred_grid <- matrix(0, nrow = N_ho, ncol = nrow(CV_fit_LSE$grid_all))
for (i in 1:nrow(CV_fit_LSE$grid_all)){
  if((CV_fit_LSE$wts[i] > 0) | (CV_fit_LP$wts[i] > 0)){
    lp_pred_grid[, i] <- Conj_lpd(X.mod = X, y.mod = WEF.in$DBH_cm,
                                  coords.mod = coords, 
                                  deltasq_pick = CV_fit_LSE$grid_all$deltasq[i],
                                  phi_pick = CV_fit_LSE$grid_all$phi[i], 
                                  nu_pick = CV_fit_LSE$grid_all$nu[i],
                                  priors, X.ho = X_ho, 
                                  y.ho = WEF.out$DBH_cm, 
                                  coords.ho = coords_ho)
  }
}

ELPD_stack_LSE = mean(log(exp(lp_pred_grid) %*% CV_fit_LSE$wts))
ELPD_stack_LP = mean(log(exp(lp_pred_grid) %*% CV_fit_LP$wts))

print(round(c(ELPD_Blm, ELPD_stack_LSE, ELPD_stack_LP), digits = 2))
# -4.55 -4.45 -4.44

#######################
### 95% CI coverage ###
#######################

# Bayesian linear regression #
CI_Blm = predictive_interval(model_bayes_lm, newdata = WEF.out, prob = 0.95)
CP_Blm = mean(CI_Blm[ , 1]<WEF.out$DBH_cm & CI_Blm[, 2]>WEF.out$DBH_cm)
CIW_Blm = mean(CI_Blm[, 2]-CI_Blm[, 1])

# stacking #
# generate posterior samples 
set.seed(123)
pos_sam_list <- list()
L = 1000
for (i in 1:nrow(CV_fit_LSE$grid_all)){
  if((CV_fit_LSE$wts[i] > 0.0001) | (CV_fit_LP$wts[i] > 0.0001)){
    cat(i, "\t")
    pos_sam_list[[i]] <- Conj_pos_sam(X.mod = X, y.mod = WEF.in$DBH_cm,
                                      coords.mod = coords, 
                                      deltasq_pick = CV_fit_LSE$grid_all$deltasq[i],
                                      phi_pick = CV_fit_LSE$grid_all$phi[i], 
                                      nu_pick = CV_fit_LSE$grid_all$nu[i],
                                      priors, X.ho = X_ho, 
                                      y.ho = WEF.out$DBH_cm, 
                                      coords.ho = coords_ho, L = L)
  }
}

sample_int_LSE <- sample(1:nrow(CV_fit_LSE$grid_all), L, 
                         replace = TRUE, prob = weights_M_LSE)
sample_int_LP <- sample(1:nrow(CV_fit_LSE$grid_all), L, 
                        replace = TRUE, prob = weights_M_LP)
table(sample_int_LSE)
table(sample_int_LP)

y_pred_sam_LSE <- matrix(0, nrow = nrow(X_ho), ncol = L)
y_pred_sam_LP <- matrix(0, nrow = nrow(X_ho), ncol = L)
for(l in 1:L){
  y_pred_sam_LSE[, l] = pos_sam_list[[sample_int_LSE[l]]]$y_U_expect[, l]
  y_pred_sam_LP[, l] = pos_sam_list[[sample_int_LP[l]]]$y_U_expect[, l]
}

CI_LSE <- apply(y_pred_sam_LSE, 1, 
                f <- function(x){quantile(x, c(0.025, 0.975))})
CI_LP <- apply(y_pred_sam_LP, 1, 
               f <- function(x){quantile(x, c(0.025, 0.975))})


CP_LSE=mean(CI_LSE[1,]<WEF.out$DBH_cm & CI_LSE[2, ]>WEF.out$DBH_cm)
CIW_LSE=mean(CI_LSE[2, ]-CI_LSE[1, ])

CP_LP=mean(CI_LP[1,]<WEF.out$DBH_cm & CI_LP[2, ]>WEF.out$DBH_cm)
CIW_LP=mean(CI_LP[2, ]-CI_LP[1, ])

print(c(CP_Blm, CP_LSE, CP_LP))
# 0.928 0.934 0.928
print(c(CIW_Blm, CIW_LSE, CIW_LP))
# 86.14428 77.91968 75.40655


library(dplyr)
library(kableExtra)
table <- matrix(c(rmspe_Blm, RMSPE_stack_LSE, RMSPE_stack_LP,
                  ELPD_Blm, ELPD_stack_LSE, ELPD_stack_LP), 
                nrow = 2, ncol = 3, byrow = TRUE)
rownames(table) <- c("RMSPE", "MLPD")
table %>%
  kbl(caption="Summary Statistics of Western experimental forest inventory data analysis",
      format="latex",
      digits = 2, 
      col.names = c("Bayesian linear regression","stacking of means",
                    "stacking of predictive densities"),
      align="c") %>%
  kable_minimal(full_width = F,  html_font = "Source Sans Pro")























