library(geoR)
library(rdist)
library(ggplot2)
rm(list = ls())

## experimental analysis of the trace of H22 ##
library(fields)
1 / Matern.cor.to.range(0.5, 0.5, cor.target=.05) # 6.0
1 / Matern.cor.to.range(0.5, 1.0, cor.target=.05) # 8.0
1 / Matern.cor.to.range(0.5, 1.5, cor.target=.05) # 9.5

samplesize_ls  = c(100, 200, 400, 800, 1600, 3200)
N_ls <- length(samplesize_ls)
phi_ls = c(6.0, 8.0, 9.5)
nu_ls = c(0.5, 1.0, 1.5)
deltasq_ls = c(1.0, 1.0, 1.0)
N_psi = length(phi_ls)
trU_d1 = matrix(NA, nrow = N_psi*2, ncol = N_ls)
trB_d1 = matrix(NA, nrow = N_psi*2, ncol = N_ls)
trC_d1 = matrix(NA, nrow = N_psi*2, ncol = N_ls)
trD_d1 = matrix(NA, nrow = N_psi*2, ncol = N_ls)
trH22_d1 = matrix(NA, nrow = N_psi*2, ncol = N_ls)

# design matrix: intercept + standard normal
set.seed(1)
X1 = cbind(1, rnorm(samplesize_ls[N_ls]))
X2 = cbind(rnorm(samplesize_ls[N_ls]), rnorm(samplesize_ls[N_ls]))
p = ncol(X1)

# inv_V_beta #
#inv_V_beta = diag(p)*1/4
V_beta = diag(p)*4
inv_V_beta = solve(V_beta)
inv_L_V_beta = solve(t(chol(V_beta)))


# d = 1 #
set.seed(2)
coords <- runif(samplesize_ls[N_ls])
for(i in 1:N_ls){
  cat("\n", i, "\t")
  N = samplesize_ls[i]
  D <- as.matrix(dist(coords[1:N])) 
  
  for (j in 1:N_psi){
    cat(j, "\t")
    phi = phi_ls[j]
    nu = nu_ls[j]
    #phi = 6.0; nu = 0.5
    deltasq = deltasq_ls[j]
    delta = sqrt(deltasq)
    R <- matern(D, phi = 1/phi, kappa = nu)
    Lphi = solve(t(chol(R)))
    
    # check assumptions in thm 2 with intercept
    xtilde = rbind(
      cbind(1/delta * X1[1:N, ], 1/delta * diag(N)), 
      cbind(inv_L_V_beta , matrix(0, nrow = p, ncol = N)),
      cbind(matrix(0, nrow = N, ncol = p), Lphi))
    U = chol2inv(chol(crossprod(xtilde)))
    B = crossprod(cbind(X1[1:N, ], diag(N)) %*% U)
    C = U[1:p, 1:p] %*% inv_V_beta %*% U[1:p, 1:p]
    DD = crossprod(Lphi %*% U[(p+1):(N+p), 1:p])
    H22 = Lphi %*% U[(p+1):(N+p), (p+1):(N+p)] %*% t(Lphi)
    
    trU_d1[j, i] = sum(diag(U[1:p, 1:p]))
    trB_d1[j, i] = sum(diag(B[1:p, 1:p]))
    trC_d1[j, i] = sum(diag(C))
    trD_d1[j, i] = sum(diag(DD))
    trH22_d1[j, i] = sum(diag(H22)) / N
    
    # check assumptions in thm 2 without intercept
    xtilde = rbind(
      cbind(1/delta * X2[1:N, ], 1/delta * diag(N)), 
      cbind(inv_L_V_beta , matrix(0, nrow = p, ncol = N)),
      cbind(matrix(0, nrow = N, ncol = p), Lphi))
    U = chol2inv(chol(crossprod(xtilde)))
    B = crossprod(cbind(X2[1:N, ], diag(N)) %*% U)
    C = U[1:p, 1:p] %*% inv_V_beta %*% U[1:p, 1:p]
    DD = crossprod(Lphi %*% U[(p+1):(N+p), 1:p])
    H22 = Lphi %*% U[(p+1):(N+p), (p+1):(N+p)] %*% t(Lphi)
    
    # check assumptions in thm 2 with intercept
    trU_d1[j+N_psi, i] = sum(diag(U[1:p, 1:p]))
    trB_d1[j+N_psi, i] = sum(diag(B[1:p, 1:p]))
    trC_d1[j+N_psi, i] = sum(diag(C))
    trD_d1[j+N_psi, i] = sum(diag(DD))
    trH22_d1[j+N_psi, i] = sum(diag(H22)) / N
    
  }
}

thm2_cond_d1_dat <- 
  data.frame(trUBCD = c(c(trU_d1), c(trB_d1), c(trC_d1), c(trD_d1)),
             N_sample = rep(rep(1:N_ls, each = N_psi*2), 4), 
             label = rep(rep(1:N_psi, N_ls*2), 4),
             type = rep(1:4, each = N_psi*N_ls*2), 
             group = rep(rep(1:2, each = N_psi), N_ls*4))

thm2_cond_d1_dat$group <- factor(thm2_cond_d1_dat$group, levels = 1:2,
                                 labels = c("with intercept", 
                                            "without intercept"))

thm2_cond_d1_dat$N_sample <- factor(thm2_cond_d1_dat$N_sample, levels = 1:N_ls,
                                    labels = paste(samplesize_ls))
thm2_cond_d1_dat$label <- factor(thm2_cond_d1_dat$label, levels = 1:N_psi,
                                 labels = c("phi=6.0, nu=0.5, deltasq=1.0",
                                            "phi=8.0, nu=1.0, deltasq=1.0",
                                            "phi=9.5, nu=1.5, deltasq=1.0"))

thm2_cond_d1_dat$type <- factor(thm2_cond_d1_dat$type, levels = 1:4,
<<<<<<< HEAD
                                labels = c("U", "B", "C", "D"))
=======
                                labels = c("tr(U[1:p, 1:p])", "tr(B[1:p, 1:p])", 
                                           "tr(C)", "tr(D)"))
>>>>>>> dd64eb2be66346625452c64fd685f1d1b13dd9c8

p1_1 <- ggplot(data = thm2_cond_d1_dat, 
               aes(x = N_sample, y = trUBCD, group = interaction(label, type),
                   color = type)) + 
  geom_line(aes(linetype = label, color = type)) +#, position = position_dodge(width=0.3)) + 
  geom_point(aes(shape = label, color = type)) + #, position = position_dodge(width=0.3)) + 
  theme_bw() + facet_grid(cols = vars(group)) +
  theme(legend.position = c(0.8, 0.58), legend.title = element_blank(), 
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key.size = unit(0.5, 'cm')) +
  xlab("sample size") + ylab("") + #ylim(0.6, 1) + 
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
ggsave(paste0("./sim/pics/check_thm1cond_d1.png"), 
       plot = p1_1, 
       width = 6.5, height = 3, units = "in", dpi = 600)


trH22_d1_dat <- 
  data.frame(trH22n = c(trH22_d1),
             N_sample = rep(1:N_ls, each = N_psi * 2), 
             label = rep(1:N_psi, N_ls * 2),
             group = rep(rep(1:2, each = N_psi), N_ls))

trH22_d1_dat$group <- factor(trH22_d1_dat$group, levels = 1:2,
                                 labels = c("with intercept", 
                                            "without intercept"))

trH22_d1_dat$N_sample <- factor(trH22_d1_dat$N_sample, levels = 1:N_ls,
                                    labels = paste(samplesize_ls))
trH22_d1_dat$label <- factor(trH22_d1_dat$label, levels = 1:N_psi,
                                 labels = c("phi=6.0, nu=0.5, deltasq=1.0",
                                            "phi=8.0, nu=1.0, deltasq=1.0",
                                            "phi=9.5, nu=1.5, deltasq=1.0"))


p1_2 <- ggplot(data = trH22_d1_dat, 
               aes(x = N_sample, y = trH22n, group = label)) + 
  # geom_errorbar(aes(ymin=L, ymax=U, linetype = label),
  #               position = position_dodge(width=0.3),
  #               width = 0.3) + 
  geom_line(aes(linetype = label), position = position_dodge(width=0.3)) + 
  geom_point(aes(shape = label), position = position_dodge(width=0.3)) + 
  theme_bw() + facet_grid(cols = vars(group)) +
  theme(legend.position = c(0.8, 0.3), legend.title = element_blank(), 
        legend.background = element_blank()) +
  xlab("sample size") + ylab("tr(H22) / n") + ylim(0.6, 1) + 
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
ggsave(paste0("./sim/pics/check_trH22_d1.png"), 
       plot = p1_2, 
       width = 6.5, height = 3, units = "in", dpi = 600)



# d = 2 #
samplesize_ls  = c(100, 250, 625, 1563, 3906, 9765) #c(100, 400, 800, 1600, 3200, 6400) #
N_ls <- length(samplesize_ls)
phi_ls = c(6.0, 8.0, 9.5)
nu_ls = c(0.5, 1.0, 1.5)
deltasq_ls = c(1.0, 1.0, 1.0)
N_psi = length(phi_ls)
trU_d2 = matrix(NA, nrow = N_psi*2, ncol = N_ls)
trB_d2 = matrix(NA, nrow = N_psi*2, ncol = N_ls)
trC_d2 = matrix(NA, nrow = N_psi*2, ncol = N_ls)
trD_d2 = matrix(NA, nrow = N_psi*2, ncol = N_ls)
trH22_d2 = matrix(NA, nrow = N_psi*2, ncol = N_ls)

# design matrix: intercept + standard normal
set.seed(1)
X1 = cbind(1, rnorm(samplesize_ls[N_ls]))
X2 = cbind(rnorm(samplesize_ls[N_ls]), rnorm(samplesize_ls[N_ls]))
p = ncol(X1)

# inv_V_beta #
#inv_V_beta = diag(p)*1/4
V_beta = diag(p)*4
inv_V_beta = solve(V_beta)
inv_L_V_beta = solve(t(chol(V_beta)))


set.seed(1)
coords <- cbind(runif(samplesize_ls[N_ls]), runif(samplesize_ls[N_ls]))

t <- proc.time()
for(i in 1:N_ls){
  cat("\n", i, "\t")
  N = samplesize_ls[i]
  D <- as.matrix(dist(coords[1:N, ])) 
  
  for (j in 1:N_psi){
    cat(j, "\t")
    phi = phi_ls[j]
    nu = nu_ls[j]
    deltasq = deltasq_ls[j]
    delta = sqrt(deltasq)
    R <- matern(D, phi = 1/phi, kappa = nu)
    Lphi = solve(t(chol(R)))
    
    # check assumptions in thm 2 with intercept
    xtilde = rbind(
      cbind(1/delta * X1[1:N, ], 1/delta * diag(N)), 
      cbind(inv_L_V_beta , matrix(0, nrow = p, ncol = N)),
      cbind(matrix(0, nrow = N, ncol = p), Lphi))
    U = chol2inv(chol(crossprod(xtilde)))
    B = crossprod(cbind(X1[1:N, ], diag(N)) %*% U)
    C = U[1:p, 1:p] %*% inv_V_beta %*% U[1:p, 1:p]
    DD = crossprod(Lphi %*% U[(p+1):(N+p), 1:p])
    H22 = Lphi %*% U[(p+1):(N+p), (p+1):(N+p)] %*% t(Lphi)
    
    trU_d2[j, i] = sum(diag(U[1:p, 1:p]))
    trB_d2[j, i] = sum(diag(B[1:p, 1:p]))
    trC_d2[j, i] = sum(diag(C))
    trD_d2[j, i] = sum(diag(DD))
    trH22_d2[j, i] = sum(diag(H22)) / N
    
    # check assumptions in thm 2 without intercept
    xtilde = rbind(
      cbind(1/delta * X2[1:N, ], 1/delta * diag(N)), 
      cbind(inv_L_V_beta , matrix(0, nrow = p, ncol = N)),
      cbind(matrix(0, nrow = N, ncol = p), Lphi))
    U = chol2inv(chol(crossprod(xtilde)))
    B = crossprod(cbind(X2[1:N, ], diag(N)) %*% U)
    C = U[1:p, 1:p] %*% inv_V_beta %*% U[1:p, 1:p]
    DD = crossprod(Lphi %*% U[(p+1):(N+p), 1:p])
    H22 = Lphi %*% U[(p+1):(N+p), (p+1):(N+p)] %*% t(Lphi)
    
    # check assumptions in thm 2 with intercept
    trU_d2[j+N_psi, i] = sum(diag(U[1:p, 1:p]))
    trB_d2[j+N_psi, i] = sum(diag(B[1:p, 1:p]))
    trC_d2[j+N_psi, i] = sum(diag(C))
    trD_d2[j+N_psi, i] = sum(diag(DD))
    trH22_d2[j+N_psi, i] = sum(diag(H22)) / N
  }
}
proc.time() - t


thm2_cond_d2_dat <- 
  data.frame(trUBCD = c(c(trU_d2), c(trB_d2), c(trC_d2), c(trD_d2)),
             N_sample = rep(rep(1:N_ls, each = N_psi*2), 4), 
             label = rep(rep(1:N_psi, N_ls*2), 4),
             type = rep(1:4, each = N_psi*N_ls*2), 
             group = rep(rep(1:2, each = N_psi), N_ls*4))

thm2_cond_d2_dat$group <- factor(thm2_cond_d2_dat$group, levels = 1:2,
                                 labels = c("with intercept", 
                                            "without intercept"))

thm2_cond_d2_dat$N_sample <- factor(thm2_cond_d2_dat$N_sample, levels = 1:N_ls,
                                    labels = paste(samplesize_ls))
thm2_cond_d2_dat$label <- factor(thm2_cond_d2_dat$label, levels = 1:N_psi,
                                 labels = c("phi=6.0, nu=0.5, deltasq=1.0",
                                            "phi=8.0, nu=1.0, deltasq=1.0",
                                            "phi=9.5, nu=1.5, deltasq=1.0"))

thm2_cond_d2_dat$type <- factor(thm2_cond_d2_dat$type, levels = 1:4,
                                labels = c("tr(U[1:p, 1:p])", "tr(B[1:p, 1:p])", 
                                           "tr(C)", "tr(D)"))

p2_1 <- ggplot(data = thm2_cond_d2_dat, 
               aes(x = N_sample, y = trUBCD, group = interaction(label, type),
                   color = type)) + 
  geom_line(aes(linetype = label, color = type)) +#, position = position_dodge(width=0.3)) + 
  geom_point(aes(shape = label, color = type)) + #, position = position_dodge(width=0.3)) + 
  theme_bw() + facet_grid(cols = vars(group)) +
  theme(legend.position = c(0.8, 0.58), legend.title = element_blank(), 
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.spacing.y = unit(0.1, 'cm'),
        legend.key.size = unit(0.5, 'cm')) +
  xlab("sample size") + ylab("") + #ylim(0.6, 1) + 
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

p2_1
ggsave(paste0("./sim/pics/check_thm1cond_d2.png"), 
       plot = p2_1, 
       width = 6.5, height = 3, units = "in", dpi = 600)


trH22_d2_dat <- 
  data.frame(trH22n = c(trH22_d2),
             N_sample = rep(1:N_ls, each = N_psi * 2), 
             label = rep(1:N_psi, N_ls * 2),
             group = rep(rep(1:2, each = N_psi), N_ls))

trH22_d2_dat$group <- factor(trH22_d2_dat$group, levels = 1:2,
                             labels = c("with intercept", 
                                        "without intercept"))

trH22_d2_dat$N_sample <- factor(trH22_d2_dat$N_sample, levels = 1:N_ls,
                                labels = paste(samplesize_ls))
trH22_d2_dat$label <- factor(trH22_d2_dat$label, levels = 1:N_psi,
                             labels = c("phi=6.0, nu=0.5, deltasq=1.0",
                                        "phi=8.0, nu=1.0, deltasq=1.0",
                                        "phi=9.5, nu=1.5, deltasq=1.0"))


p2_2 <- ggplot(data = trH22_d2_dat, 
               aes(x = N_sample, y = trH22n, group = label)) + 
  # geom_errorbar(aes(ymin=L, ymax=U, linetype = label),
  #               position = position_dodge(width=0.3),
  #               width = 0.3) + 
  geom_line(aes(linetype = label), position = position_dodge(width=0.3)) + 
  geom_point(aes(shape = label), position = position_dodge(width=0.3)) + 
  theme_bw() + facet_grid(cols = vars(group)) +
<<<<<<< HEAD
  theme(legend.position = c(0.8, 0.3), legend.title = element_blank(), 
=======
  theme(legend.position = c(0.85, 0.3), legend.title = element_blank(), 
>>>>>>> dd64eb2be66346625452c64fd685f1d1b13dd9c8
        legend.background = element_blank()) +
  xlab("sample size") + ylab("tr(H22) / n") + ylim(0.6, 1) + 
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
p2_2
ggsave(paste0("./sim/pics/check_trH22_d2.png"), 
       plot = p2_2, 
       width = 6.5, height = 3, units = "in", dpi = 600)









