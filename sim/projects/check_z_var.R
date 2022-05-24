library(geoR)
library(rdist)
library(ggplot2)
rm(list = ls())

## experimental analysis of the posterior variance E2 ##
library(fields)
1 / Matern.cor.to.range(0.5, 0.5, cor.target=.05) # 6.0
1 / Matern.cor.to.range(0.5, 1.0, cor.target=.05) # 8.0
1 / Matern.cor.to.range(0.5, 1.5, cor.target=.05) # 9.5

samplesize_ls  = c(100, 200, 400, 800, 1600, 3200)
N_ls <- length(samplesize_ls)
pos_var_z_d1 = list()
phi_ls = c(6.0, 8.0, 9.5)
nu_ls = c(0.5, 1.0, 1.5)
N_psi = length(phi_ls)

# d = 1 #
set.seed(2)
coords <- runif(samplesize_ls[N_ls])
for(i in 1:N_ls){
  cat("\n", i, "\t")
  N = samplesize_ls[i]
  D <- as.matrix(dist(coords[1:N])) 
  pos_var_z_d1[[i]] = matrix(NA, nrow = N, ncol = N_psi)
  
  for (j in 1:N_psi){
    cat(j, "\t")
    phi = phi_ls[j]
    nu = nu_ls[j]
    #phi = 6.0; nu = 0.5
    deltasq = 1.0
    R <- matern(D, phi = 1/phi, kappa = nu)
    cholR = chol(R)
    # diagonal elements of the covariance of posterior distribution of z
    pos_var_z_d1[[i]][, j] <- diag(solve(deltasq^{-1}*diag(N) + chol2inv(cholR)))
  }
}

E2_d1_med_dat <- 
  data.frame(median = c(sapply(pos_var_z_d1, 
                               f <- function(x)(apply(x, 2, median)))),
             U = c(sapply(
               pos_var_z_d1, f <- function(x){
                 apply(x, 2, g <- function(x){quantile(x, 0.975)})})),
             L = c(sapply(
               pos_var_z_d1, f <- function(x){
                 apply(x, 2, g <- function(x){quantile(x, 0.025)})})),
             N_sample = rep(1:N_ls, each = N_psi), 
             label = rep(1:N_psi, N_ls))

E2_d1_med_dat$N_sample <- factor(E2_d1_med_dat$N_sample, levels = 1:N_ls,
                                 labels = paste(samplesize_ls))
E2_d1_med_dat$label <- factor(E2_d1_med_dat$label, levels = 1:N_psi,
                              labels = c("phi=6.0, nu=0.5, deltasq=1.0",
                                         "phi=8.0, nu=1.0, deltasq=1.0",
                                         "phi=9.5, nu=1.5, deltasq=1.0"))

p1_1 <- ggplot(data = E2_d1_med_dat, 
               aes(x = N_sample, y = median, group = label)) + 
  geom_errorbar(aes(ymin=L, ymax=U, linetype = label),
                position = position_dodge(width=0.3),
                width = 0.3) + 
  geom_line(aes(linetype = label), position = position_dodge(width=0.3)) + 
  geom_point(aes(shape = label), position = position_dodge(width=0.3)) + 
  theme_bw() +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank(), 
        legend.background = element_blank()) +
  xlab("sample size") + ylab("variance of z") + 
  scale_shape_discrete(
    name  ="",
    breaks=c("phi=6.0, nu=0.5, deltasq=1.0", "phi=8.0, nu=1.0, deltasq=1.0",
             "phi=9.5, nu=1.5, deltasq=1.0"),
    labels=c(expression(phi ~ "="~ 6.0 ~ "," ~ nu ~ "=" ~ 0.5 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0), 
             expression(phi ~ "="~ 8.0 ~ "," ~ nu ~ "=" ~ 1.0 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0),
             expression(phi ~ "="~ 9.5 ~ "," ~ nu ~ "=" ~ 1.5 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0))) +
  scale_linetype_discrete(
    name  ="",
    breaks=c("phi=6.0, nu=0.5, deltasq=1.0", "phi=8.0, nu=1.0, deltasq=1.0",
             "phi=9.5, nu=1.5, deltasq=1.0"),
    labels=c(expression(phi ~ "="~ 6.0 ~ "," ~ nu ~ "=" ~ 0.5 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0), 
             expression(phi ~ "="~ 8.0 ~ "," ~ nu ~ "=" ~ 1.0 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0),
             expression(phi ~ "="~ 9.5 ~ "," ~ nu ~ "=" ~ 1.5 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0))) 
p1_1
ggsave(paste0("./sim/pics/check_z_var_d1.png"), 
       plot = p1_1, 
       width = 6.5, height = 3, units = "in", dpi = 600)



# d = 2 #
samplesize_ls  = c(100, 300, 900, 2700, 8100, 24300)
N_ls <- length(samplesize_ls)
pos_var_z_d2 = list()
phi_ls = c(6.0, 8.0, 9.5)
nu_ls = c(0.5, 1.0, 1.5)
N_psi = length(phi_ls)

set.seed(1)
coords <- cbind(runif(samplesize_ls[N_ls]), runif(samplesize_ls[N_ls]))
t <- proc.time()
for(i in 1:N_ls){
  cat("\n", i, "\t")
  N = samplesize_ls[i]
  pos_var_z_d2[[i]] = matrix(NA, nrow = N, ncol = N_psi)
  D <- as.matrix(dist(coords[1:N, ])) 
  
  for (j in 1:N_psi){
    cat(j, "\t")
    phi = phi_ls[j]
    nu = nu_ls[j]
    deltasq = 1.0
    R <- matern(D, phi = 1/phi, kappa = nu)
    cholR = chol(R)
    # diagonal elements of the covariance of posterior distribution of z
    pos_var_z_d2[[i]][, j] <- diag(solve(deltasq^{-1}*diag(N) + chol2inv(cholR)))
  }
}
proc.time() - t


E2_d2_med_dat <- 
  data.frame(median = c(sapply(pos_var_z_d2, 
                               f <- function(x)(apply(x, 2, median)))),
             U = c(sapply(
               pos_var_z_d2, f <- function(x){
                 apply(x, 2, g <- function(x){quantile(x, 0.975)})})),
             L = c(sapply(
               pos_var_z_d2, f <- function(x){
                 apply(x, 2, g <- function(x){quantile(x, 0.025)})})),
             N_sample = rep(1:N_ls, each = N_psi), 
             label = rep(1:N_psi, N_ls))

E2_d2_med_dat$N_sample <- factor(E2_d2_med_dat$N_sample, levels = 1:N_ls,
                             labels = paste(samplesize_ls))
E2_d2_med_dat$label <- factor(E2_d2_med_dat$label, levels = 1:N_psi,
                          labels = c("phi=6.0, nu=0.5, deltasq=1.0",
                                     "phi=8.0, nu=1.0, deltasq=1.0",
                                     "phi=9.5, nu=1.5, deltasq=1.0"))

p1_2 <- ggplot(data = E2_d2_med_dat, 
               aes(x = N_sample, y = median, group = label)) + 
  geom_errorbar(aes(ymin=L, ymax=U, linetype = label),
                position = position_dodge(width=0.3),
                width = 0.3) + 
  geom_line(aes(linetype = label), position = position_dodge(width=0.3)) + 
  geom_point(aes(shape = label), position = position_dodge(width=0.3)) + 
  theme_bw() +
  theme(legend.position = c(0.6, 0.8), legend.title = element_blank(), 
        legend.background = element_blank()) +
  xlab("sample size") + ylab("variance of z") + 
  scale_shape_discrete(
    name  ="",
    breaks=c("phi=6.0, nu=0.5, deltasq=1.0", "phi=8.0, nu=1.0, deltasq=1.0",
             "phi=9.5, nu=1.5, deltasq=1.0"),
    labels=c(expression(phi ~ "="~ 6.0 ~ "," ~ nu ~ "=" ~ 0.5 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0), 
             expression(phi ~ "="~ 8.0 ~ "," ~ nu ~ "=" ~ 1.0 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0),
             expression(phi ~ "="~ 9.5 ~ "," ~ nu ~ "=" ~ 1.5 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0))) +
  scale_linetype_discrete(
    name  ="",
    breaks=c("phi=6.0, nu=0.5, deltasq=1.0", "phi=8.0, nu=1.0, deltasq=1.0",
             "phi=9.5, nu=1.5, deltasq=1.0"),
    labels=c(expression(phi ~ "="~ 6.0 ~ "," ~ nu ~ "=" ~ 0.5 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0), 
             expression(phi ~ "="~ 8.0 ~ "," ~ nu ~ "=" ~ 1.0 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0),
             expression(phi ~ "="~ 9.5 ~ "," ~ nu ~ "=" ~ 1.5 ~ 
                          "," ~ delta^2 ~ "=" ~ 1.0))) 
p1_2
ggsave(paste0("./sim/pics/check_z_var_d2.png"), 
       plot = p1_2, 
       width = 6.5, height = 3, units = "in", dpi = 600)


# E2_d2_dat <- 
#   data.frame(pos_var_z = unlist(pos_var_z_d2),
#              N_sample =  unlist(sapply(1:N_ls, f <- function(x){
#                rep(x, samplesize_ls[x]*N_psi)})),
#              label = unlist(sapply(1:N_ls, f <- function(x){
#                rep(1:N_psi, each = samplesize_ls[x])})))
# E2_d2_dat$N_sample <- factor(E2_d2_dat$N_sample, levels = 1:N_ls,
#                              labels = paste(samplesize_ls))
# E2_d2_dat$label <- factor(E2_d2_dat$label, levels = 1:N_psi,
#                           labels = c("phi=6.0, nu=0.5, deltasq=1.0",
#                                      "phi=9.0, nu=1.0, deltasq=1.0"))
# p2 <- ggplot(E2_d2_dat, aes(x = N_sample, y = pos_var_z, color = label)) +
#   geom_boxplot()  + theme_bw() + ylim(c(0, 0.5)) +
#   # geom_line(data = E2_d2_med_dat, aes(x = N_sample, y = median, color = label), 
#   #           group= E2_d2_med_dat$label, position = position_dodge(width=0.75)) + 
#   # geom_point(data = E2_d2_med_dat, aes(x = N_sample, y = median, color = label), 
#   #            position = position_dodge(width=0.75))+
#   theme(legend.position = c(0.45, 0.9), legend.title = element_blank(), 
#         legend.background = element_blank()) +
#   xlab("sample size") + ylab("variance of z") +
#   scale_colour_manual(values=c("#E69F00", "#56B4E9"))
# p2
# #time 616.429












