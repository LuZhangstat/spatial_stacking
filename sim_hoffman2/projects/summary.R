rm(list = ls())
library(ggplot2)
library(MASS)
library(spBayes)
library(ggplot2)
library("gridExtra")

# colorblind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sim_ind = 1 # simulation index 1 or 2


load(paste0("./sim_hoffman2/results/sim", sim_ind, "_1.Rdata"))
K_fold = 10
N_list = length(samplesize_ls)
N_sim = 60
SPE_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_M0 = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_M0 = matrix(NA, nrow = N_sim, ncol = N_list)
ELPD_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_stack_LP = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_stack_LSE = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_M0 = matrix(NA, nrow = N_sim, ncol = N_list)
SPE_w_MCMC = matrix(NA, nrow = N_sim, ncol = N_list)

## recover the estimated hyper-parameters ##
deltasq_grid <- c(0.1, 0.5, 1, 2)
phi_grid = c(3, 9, 15, 21)   #3/(0.6*sqrt(2)) to 3/(0.1*sqrt(2)) phi_grid = c(3, 13, 23, 33)
nu_grid = c(0.5, 1, 1.5, 1.75)
grid_all <- expand.grid(deltasq_grid, phi_grid, nu_grid)
colnames(grid_all) <- c("deltasq", "phi", "nu")
weights_M_LSE_all = array(0, dim = c(length(deltasq_grid) * length(phi_grid) * 
                                       length(nu_grid), N_list, N_sim))
weights_M_LP_all = array(0, dim = c(length(deltasq_grid) * length(phi_grid) * 
                                      length(nu_grid), N_list, N_sim))


for(i in 1:N_sim){
  filename <- paste0("./sim_hoffman2/results/sim", sim_ind, "_", i, ".Rdata")
  load(filename)
  
  SPE_stack_LP[i, ] = DIV_matrix[, "SPE_stack_LP"]
  SPE_stack_LSE[i, ] = DIV_matrix[, "SPE_stack_LSE"]
  SPE_M0[i, ] = DIV_matrix[, "SPE_M0"]
  SPE_MCMC[i, ] = DIV_matrix[, "SPE_MCMC"]
  ELPD_stack_LSE[i, ] = DIV_matrix[, "ELPD_stack_LSE"]
  ELPD_stack_LP[i, ] = DIV_matrix[, "ELPD_stack_LP"]
  ELPD_M0[i, ] = DIV_matrix[, "ELPD_M0"]
  ELPD_MCMC[i, ] = DIV_matrix[, "ELPD_MCMC"]
  SPE_w_stack_LP[i, ] = DIV_matrix[, "SPE_w_stack_LP"]
  SPE_w_stack_LSE[i, ] = DIV_matrix[, "SPE_w_stack_LSE"]
  SPE_w_M0[i, ] = DIV_matrix[, "SPE_w_M0"]
  SPE_w_MCMC[i, ] = DIV_matrix[, "SPE_w_MCMC"]
  
  
  weights_M_LSE_all[, , i] <- weights_M_LSE
  weights_M_LP_all[, , i] <- weights_M_LP
  
}


####

type = c("M0", "stacking of means", "stacking of predictive distributions", 
         "MCMC")

test = c("MSPE", "MSEZ", "MLPD")


dat_check <- data.frame(N_sample = rep(rep(rep(paste(samplesize_ls), 
                                               each = N_sim), 
                                       length(type)), length(test)),
                        value = c(c(c(SPE_M0), c(SPE_stack_LSE), 
                                    c(SPE_stack_LP), c(SPE_MCMC)),
                                c(c(SPE_w_M0), c(SPE_w_stack_LSE), 
                                  c(SPE_w_stack_LP), c(SPE_w_MCMC)),
                                c(c(ELPD_M0), c(ELPD_stack_LSE), 
                                 c(ELPD_stack_LP), c(ELPD_MCMC))),
                        label = rep(rep(c(1:4), each = N_list * N_sim), 
                                    length(test)),
                        test = rep(1:3, each = N_list * N_sim * length(type)))

dat_check$label <- factor(dat_check$label, levels = 1:4,
                          labels = type)

dat_check$test <- factor(dat_check$test, levels = 1:3,
                          labels = test)

p_summary <- ggplot(dat_check, aes(x = N_sample, y = value, color = label)) +
  geom_violin(draw_quantiles = c(0.5)) +theme_bw() + xlab("sample size") + 
  facet_wrap(~ test, ncol = 1, scales = "free_y", strip.position="right") +
  theme(legend.position="top", legend.title = element_blank()) + ylab(" ") +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"))
  
ggsave(paste0("./sim_hoffman2/pics/CVexperiment_sim", sim_ind, ".png"), 
       plot = p_summary, 
       width = 6.5, height = 4.5, units = "in", dpi = 600)


# On average, only 3.5 out of 64 models have no-zero weights
weights_nonzero_LSE = (weights_M_LSE_all > 0.001)
sum(weights_nonzero_LSE) / (64 * 8) # 3.6
weights_nonzero_LP = (weights_M_LP_all > 0.001)
sum(weights_nonzero_LP) / (64 * 8) # 25.24

weight_data <- data.frame(
  nonzero_count = c(c(apply(weights_nonzero_LSE, 3:2, sum)), 
                    c(apply(weights_nonzero_LP, 3:2, sum))),
  N_sample = rep(rep(paste(samplesize_ls), each = N_sim), 2),
  label = rep(1:2, each = N_list * N_sim)
)

weight_data$label <- factor(weight_data$label, levels = 1:2,
                          labels = c("stacking of means", 
                                     "stacking of predictive densities"))

p_nonzero_counts <- 
  ggplot(weight_data, aes(x = N_sample, y = nonzero_count, color = label)) +
  geom_violin(draw_quantiles = c(0.5))  + theme_bw() + ylim(c(0, 40)) +
  theme(legend.position = c(0.8, 0.4), legend.title = element_blank()) +
  xlab("sample size") + ylab("No. of nonzero weights") +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"))
p_nonzero_counts 

## Obviously, stacking based on LP has more nonzero weights than stacking based on LSE
ggsave(paste0("./sim_hoffman2/pics/nonzero_check_sim", sim_ind, ".png"), 
       plot = p_nonzero_counts, 
       width = 6.5, height = 3, units = "in", dpi = 600)


## check the inference of hyper-parameters ##
expect_w_phi <- matrix(0, nrow = length(phi_grid), ncol = N_sim * N_list)
for(i in 1:length(phi_grid)){
  expect_w_phi[i, ] <- c(apply(weights_M_LP_all[which(grid_all$phi == phi_grid[i]), , ], 2:3, sum))
}
summary(c(t(expect_w_phi) %*% phi_grid))
hist(t(expect_w_phi) %*% phi_grid )
plot(phi_grid, rowMeans(expect_w_phi))


expect_w_phi <- matrix(0, nrow = length(phi_grid), ncol = N_sim)
for(i in 1:length(phi_grid)){
  expect_w_phi[i, ] <- c(colSums(weights_M_LP_all[which(grid_all$phi == phi_grid[i]), 2, ]))
}
summary(c(t(expect_w_phi) %*% phi_grid))
hist(t(expect_w_phi) %*% phi_grid )
#plot(phi_grid, rowMeans(expect_w_phi))


apply(weights_M_LSE_all[which(grid_all$phi == phi_grid[i]), , ], 2:3, sum)


## plot the interpolated map of the latent process ##
load(paste0("./sim_hoffman2/results/sim", sim_ind, "_10.Rdata"))
## check the plots of latent process ##
library(coda)
library(spBayes)
library(MBA)
library(fields)
library(classInt)
library(RColorBrewer)
library(sp)

r = 8
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

png(paste0("./sim_hoffman2/pics/w_all_sim", sim_ind, ".png"), 
    width = 600, height = 400, units = "px", pointsize = 16)
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
     xlab="x", main = "stacking-EP") 
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
dev.off()

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

png(paste0("./sim_hoffman2/pics/w_pred_sim", sim_ind, ".png"), 
    width = 600, height = 400, units = "px", pointsize = 16)
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
dev.off()

