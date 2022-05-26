rm(list=ls())
library(MBA)
library(fields)
library(classInt)
library(geoR)
library(spBayes)
data("WEF.dat")
#WEFraw.dat = WEF.dat

###################################################
### WEF data
###################################################
#WEF.dat=read.csv("./RDA/data/WEF_data/WEFsmall.csv")

set.seed(1)
indexNA <- (is.na(WEF.dat$East_m) | is.na(WEF.dat$North_m) | 
              is.na(WEF.dat$DBH_cm) | is.na(WEF.dat$Species) | 
              WEF.dat$Live_dead == "D" | WEF.dat$Species == "NF") 
WEF.dat <- droplevels(
  WEF.dat[!indexNA, c("East_m", "North_m", "DBH_cm", "Species")])

table(as.numeric(WEF.dat$Species))
levels(WEF.dat$Species)
ind=sample(1:nrow(WEF.dat), 500, replace=FALSE)

### holdout data to assess RMSPE ###
WEF.out=WEF.dat[ind,]
WEF.in=WEF.dat[-ind,]
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


### Linear regression ###
lm.DBH <- lm(DBH~Species, data=WEF.in)
summary(lm.DBH)
DBH.resid <- resid(lm.DBH)

surf <- mba.surf(cbind(coords,DBH.resid), no.X=100, no.Y=100, h=5, m=2, extend=FALSE)$xyz.est
#dev.new()
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", col=col.br(25))


### variogram of raw data and residuals ###
max.dist=0.5*max(rdist(coords))
bins=20

#dev.new()
vario.DBH <- variog(coords=coords, data=DBH, uvec=(seq(5, max.dist, length=bins)))
plot(vario.DBH)


#dev.new()
vario.DBH.resid <- variog(coords=coords, data=DBH.resid, uvec=(seq(0, max.dist, length=bins)))
plot(vario.DBH.resid)


### spatial mle ###
t <- proc.time()
mle <- likfit(coords = coords, data = DBH, 
              trend = trend.spatial(~Species, WEF.in), 
              ini.cov.pars=c(100,100),
              nugget = 300, kappa = 0.5, fix.kappa = FALSE,
              cov.model="matern", nospatial=TRUE)
proc.time() - t

mle

mle$AIC
mle$BIC

mle$nospatial$AIC
mle$nospatial$BIC

### RMSPE ###
krig_mlefit=krige.conv(coords=coords, data=DBH,
                       locations=WEF.out[,c("East_m","North_m")], 
                       krige=krige.control(type.krige="OK", obj.model=mle, 
                                           trend.d=trend.spatial(~Species,WEF.in),
                                           trend.l=trend.spatial(~Species,WEF.out)))

pred_spatial=krig_mlefit$predict
rmspe_spatial=sqrt(mean((pred_spatial-WEF.out$DBH_cm)^2))

pred_lm=as.vector(as.matrix(trend.spatial(~Species,WEF.out))%*%lm.DBH$coefficients)
rmspe_lm=sqrt(mean((pred_lm-WEF.out$DBH_cm)^2))

rmspe_spatial
rmspe_lm


### CP ###
CI_spatial=pred_spatial+1.96*sqrt(krig_mlefit$krige.var)%*%t(c(-1,1))  ## confidence interval ##
CP_spatial=mean(CI_spatial[,1]<WEF.out$DBH_cm & CI_spatial[,2]>WEF.out$DBH_cm) ## coverage probability ##
CIW_spatial=mean(CI_spatial[,2]-CI_spatial[,1]) ## confidence interval width ##

CP_spatial
CIW_spatial

m_r <- pred_spatial - WEF.out$DBH_cm 
m_r <- m_r / sqrt(krig_mlefit$krige.var) 
m_r <- -0.5 * (m_r^2 + log(2 * pi) + log(krig_mlefit$krige.var))
ELPD_mle <- mean(m_r)

N_ho=nrow(WEF.out)
CI_lm=pred_lm+1.96*summary(lm.DBH)$sigma*cbind(-rep(1,N_ho),rep(1,N_ho))
CP_lm=mean(CI_lm[,1]<WEF.out$DBH_cm & CI_lm[,2]>WEF.out$DBH_cm)
CIW_lm=mean(CI_lm[,2]-CI_lm[,1])

CP_lm
CIW_lm

## stacking ##
source("utils.R")
# pick candidate values of phi #
range(1/Matern.cor.to.range(0.6 * max(rdist(coords)), 0.5, cor.target=.05),  # 0.01417578
      1/Matern.cor.to.range(0.1 * max(rdist(coords)), 0.5, cor.target=.05), # 0.08505468
      1/Matern.cor.to.range(0.6 * max(rdist(coords)), 1.75, cor.target=.05), # 0.02397907
      1/Matern.cor.to.range(0.1 * max(rdist(coords)), 1.75, cor.target=.05)) # 0.1438744

# candidate values for the hyperparameters #
deltasq_grid <- c(0.1, 0.5, 1, 2)
phi_grid = seq(0.014, 0.14, length.out = 4)
nu_grid = c(0.5, 1, 1.5, 1.75)

# precomputation and priors #
X <- model.matrix(~Species, data = WEF.in) # design matrix 
p = ncol(X)
priors <- list(mu_beta = rep(0, p),
               inv_V_beta = 1/400 * diag(p),
               a_sigma = 2,
               b_sigma = 300)
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

# prediction #
# design matrix and locations of the data for prediction
X_ho <- model.matrix(~Species, data = WEF.out)
coords_ho <- as.matrix(WEF.out[,c("East_m","North_m")])
## stacking mean squared prediction error ##
y_pred_grid <- matrix(0, nrow = nrow(X_ho), ncol = nrow(CV_fit_LSE$grid_all))
w_expect_grid <- matrix(0, nrow = nrow(X_ho) + nrow(X), 
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


## stacking Expected log pointwise predictive density ##
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

## generate posterior samples ##
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

# generate posterior samples 
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

CP_LSE
CIW_LSE

CP_LP=mean(CI_LP[1,]<WEF.out$DBH_cm & CI_LP[2, ]>WEF.out$DBH_cm)
CIW_LP=mean(CI_LP[2, ]-CI_LP[1, ])

CP_LP
CIW_LP






