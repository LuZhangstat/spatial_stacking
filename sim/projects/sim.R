rm(list = ls())
source("utils.R")
load("./sim/data/simdata.RData")


## fit conjugate spatial model ##
alpha_grid = c(0.25, 0.5, 0.75)
phi_grid_0 = c(0.1, 0.3, 0.8)
D_max = sqrt(2)
phi_grid = 3 / (D_max * phi_grid_0)
priors <- list(mu_beta = rep(0, length(beta)),
               inv_V_beta = 1/4 * diag(length(beta)),
               a_sigma = 2,
               b_sigma = 2)
L = 200

set.seed(1234)



results <- sp_stacking(X, y, coords, alpha_grid, phi_grid, priors, L,
                       label = "LSE")
results$wts
results$grid_all

set.seed(1234)
results <- sp_stacking(X, y, coords, alpha_grid, phi_grid, priors, L, 
                       label = "LP")
results$wts
results$grid_all

# check the approximation 
sapply(1:9, function(ind){
  round(mean(results$post_sam_ls[[ind]]$sigma.sq.sam) * 
             results$grid_all$deltasq[ind], 2)
})
sapply(1:9, function(ind){
  round(mean(results$post_sam_ls[[ind]]$sigma.sq.sam), 2)
})
sapply(1:9, function(ind){
  round(abs(mean(results$post_sam_ls[[ind]]$sigma.sq.sam) * 
              results$grid_all[ind, 2]), 2)
})
phi*sigma.sq
results$wts
sapply(lp_pred_ls, function(x){summary(x$pareto_k)})
sapply(results$y_pred_ls, function(x){summary(x$pareto_k)})
summary(results$y_pred_ls[[3]]$pareto_k)
summary(results$y_pred_ls[[3]]$n_eff)

deltasq2alpha(results$grid_all[3, 1]) # 0.25
((3 / results$grid_all[3, 2]) / D_max) # 0.75

## check the loo prediction ##
#--------------------- fitted w compare ---------------------#
library(coda)
library(spBayes)
library(MBA)
library(fields)
library(classInt)
library(RColorBrewer)
library(sp)

h <- 12
surf.raw <- mba.surf(cbind(coords, y), no.X=300, no.Y=300, 
                     exten=TRUE, sp=TRUE, h=h)$xyz.est


surf.brks <- classIntervals(surf.raw$z,500, 'pretty')$brks
col.pal <- colorRampPalette(brewer.pal(11,'RdBu')[11:1])
xlim <- c(0,1.13)
zlim <- range(c(y, c(Y_hat)))

# size for the mapping of w               
par(mfrow = c(3, 4))

##Obs
i <- as.image.SpatialGridDataFrame(surf.raw)
plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", xlab="x", 
     main = "raw") 
#main = "true")
axis(2, las=1)
axis(1)
image.plot(i, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)

for(ind in 1:9){
  
  surf.fit <- mba.surf(cbind(coords, Y_hat[, ind]), no.X=300, no.Y=300,
                        exten=TRUE, sp=TRUE, h=h)$xyz.est
  
  
  i <- as.image.SpatialGridDataFrame(surf.fit)
  plot(coords, typ="n", cex=0.5, xlim=xlim, axes=FALSE, ylab="y", xlab="x",
       main = paste(ind))
  axis(2, las=1)
  axis(1)
  image.plot(i, add=TRUE, col=rev(col.pal(length(surf.brks)-1)), zlim=zlim)
}

colSums((Y_hat - y)^2)/N




## second round ##
alpha_grid = c(0.25 / 2, 0.25, (0.5 + 0.25) / 2)
phi_grid_0 = c((0.5 + 0.75) / 2, 0.75, (0.75 + 1) / 2)
phi_grid = 3 / (D_max * phi_grid_0)
set.seed(1234)
results <- sp_stacking(X, y, coords, alpha_grid, phi_grid, priors, L)
results$wts
results$grid_all






