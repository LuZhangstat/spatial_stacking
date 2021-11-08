## util functions ##
library(invgamma)
library(Matrix)
library(loo)
library(rbenchmark) # test the efficiency

alpha2deltasq <- function(alpha){
  # change alpha \in (0, 1) into deltasq \in (0, \infty) #
  if(alpha <= 0 | alpha >1){
    stop("alpha should be in (0, 1]")
  }
  deltasq <- (1/alpha - 1)
  return(deltasq)
}

sp_stacking <- function(X, y, coords, alpha_grid, phi_grid, priors, L){
  
  # L: number of samples for each model #
  
  # generate all candidates of the fixed parameter
  deltasq_grid <- sapply(alpha_grid, alpha2deltasq)
  grid_all <- expand.grid(deltasq_grid, phi_grid)
  colnames(grid_all) <- c("deltasq", "phi")
  
  # # pre-computation
  D <- as.matrix(dist(coords)) # distance matrix
  N = length(y)
  p = ncol(X) # number of predictors, including intercept
  
  post_sam_ls <- list()
  y_pred_ls <- list()
  for (i in 1:nrow(grid_all)){
    cat(i, "\t")
    post_sam_ls[[i]] <- conj_sp_sample(X, y, D, N, p, 
                               deltasq_pick = grid_all$deltasq[i], 
                               phi_pick = grid_all$phi[i], priors, L)
    y_pred_ls[[i]] <- IS_y_expect(X, y, post_sam_ls[[i]], N, p, 
                          deltasq_pick = grid_all$deltasq[i], L)
  }
  Y_hat <- sapply(y_pred_ls, function(x){x$y_expect})
  
  
}

# deltasq_pick <- grid_all$deltasq[1]
# phi_pick <- grid_all$phi[1]
# L = 100

conj_sp_sample <- function(X, y, D, N, p, deltasq_pick, phi_pick, priors, L){
  ## generate posterior samples of the conjugate spatial model ##
  
  # pre-computation #
  chol_inv_V = chol(rbind(cbind(crossprod(X) / deltasq_pick + priors$inv_V_beta, 
                                t(X) / deltasq_pick),
                          cbind(X / deltasq_pick, 
                                chol2inv(chol(exp(- phi_pick * D))) +
                                  diag(rep(1/deltasq_pick, N)))))
  
  inv_V_mu_beta = priors$inv_V_beta %*% priors$mu_beta
  
  mu_star = c(inv_V_mu_beta + crossprod(X, y) / deltasq_pick, y / deltasq_pick)
  
  mu_star = forwardsolve(chol_inv_V, mu_star, transpose = TRUE, upper.tri = TRUE)
  
  b_star = priors$b_sigma + 0.5 * (sum(y^2) / deltasq_pick + 
                                     sum(priors$mu_beta * inv_V_mu_beta) - 
                                     sum(mu_star^2))
  
  a_star = priors$a_sigma + N / 2
  
  
  ## generate posterior samples ##
  sigma.sq.sam =  rinvgamma(L, shape = a_star, rate = b_star) 
  # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
  gamma.sam = matrix(rnorm((N + p) * L), nrow = N + p, ncol = L) %*% 
    Diagonal(n = L, sqrt(sigma.sq.sam)) + mu_star
  gamma.sam = backsolve(chol_inv_V, gamma.sam)
  
  return(list(sigma.sq.sam = sigma.sq.sam, gamma.sam = gamma.sam))
}

IS_y_expect <- function(X, y, post_sam, N, p, deltasq_pick, L){
  ## generate the LOO expectation ##
  
  # expectation
  M_expect <-  X %*% post_sam$gamma.sam[1:p, ] + 
    post_sam$gamma.sam[(p+1):(N+p), ]
  # the matrix of log ratios (-lp)
  M_r <- M_expect - y
  M_r <- M_r %*% Diagonal(n = L, sqrt(1 / (deltasq_pick * post_sam$sigma.sq.sam)))
  M_r <- M_r^2 + tcrossprod(rep(1, N), log(deltasq_pick * post_sam$sigma.sq.sam))
  # weight matrix
  M_w <- apply(M_r, 1, function(x){
    psis_result <- psis(x, r_eff = NA)
    w <- c(weights(psis_result, log = FALSE))
    return(c(psis_result$diagnostics$pareto_k, psis_result$diagnostics$n_eff,
             w))
  })
  y_expect <- rowSums(M_expect * t(M_w[-(1:2), ]))
  
  return(list(y_expect = y_expect, pareto_k = M_w[1, ], n_eff = M_w[2, ]))
}

lagrange_stacking_weight <- function(Y_hat, y){
  K = ncol(Y_hat)
  YhatTYhat = crossprod(Y_hat)
  YhatTy = crossprod(Y_hat, y)
  CholYTY = chol(YhatTYhat)
  u1 = forwardsolve(CholYTY, rep(1, K), upper.tri = TRUE, transpose = TRUE)
  u2 =  forwardsolve(CholYTY, YhatTy, upper.tri = TRUE, transpose = TRUE)
  lambda = 2 * (1 - sum(u2 * u1))/sum(u1^2)
  u3 = forwardsolve(CholYTY, rep(lambda, K) + 2 * YhatTy, upper.tri = TRUE, 
                    transpose = TRUE) 
  w = 0.5 * backsolve(CholYTY, u3)
  return(w)
}

