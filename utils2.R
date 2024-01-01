## util functions ##
library(invgamma)
library(Matrix)
library(quadprog)   # Quadratic Programming Problems
library(cvTools) #run the above line if you don't have this library
library(rdist)
library(CVXR)  # stacking weights computation
library(Rmosek) # need to install mosek 
library(GPBayes)

QP_stacking_weight <- function(Y_hat, y){
  
  ## compute stacking weights based on expectation
  
  stopifnot(is.matrix(Y_hat))
  N <- nrow(Y_hat)
  K <- ncol(Y_hat)
  if (K < 2) {
    stop("At least two models are required for stacking weights.")
  }
  
  y_tilde = y - Y_hat[, K];
  Y_tilde = Y_hat[, -K] - Y_hat[, K]
  
  ## minimize: y_tilde^T y_tilde - 2 y_tilde^T Y_tilde w  + 1/2 w^T 2Y_tilde^\top Y_tilde w
  ## under the constraints: A^T b >= b0
  ## with b0 = (0,...0, -1)^T
  ## and (1  0)
  ## A = (0 1 )
  ##     (-1 -1)
  ## we can use solve.QP.compact as follows:
  ##
  Dmat <- 2 * crossprod(Y_tilde)
  diag(Dmat) <- diag(Dmat) + 0.0001 * min(diag(Dmat))
  #Dmat = Dmat + diag(K-1) * 3  + tcrossprod(rep(1, K-1)) * 3
  dvec <- 2 * crossprod(Y_tilde, y_tilde) 
  #dvec = dvec + rep(1, K-1) * 3
  Amat <- cbind(diag(K-1), rep(-1, K - 1))
  bvec <- c(rep(0, K-1), -1)
  w <- solve.QP(Dmat,dvec,Amat, bvec=bvec)$solution
  wts <- structure(
    c(w, 1 - sum(w)),
    names = paste0("model", 1:K),
    class = c("stacking_weights")
  )
  
  return(wts)
}


stacking_weights <- function(lpd_point){
  
  lpd_m = mean(lpd_point)
  lpd_point = lpd_point - lpd_m ## rescale the log-density for numerical stability
  exp_lpd_point <- exp(lpd_point)
  G <- ncol(lpd_point)
  
  w <- Variable(G)
  obj <- sum(log(exp_lpd_point %*% w))
  constr <- list(sum(w) == 1, w >= 0)
  prob <- Problem(Maximize(obj), constr)
  result <- psolve(prob, solver = "MOSEK")
  
  wts <- structure(
    result$getValue(w)[,1],
    names = paste0("model", 1:G),
    class = c("stacking_weights")
  )
  
  return(wts)
}

# modified based on https://github.com/stan-dev/loo/blob/e9f93fab6d408d8cb6831b237424e33374c4ca0c/R/loo_model_weights.R#L240
stacking_weights_old <-
  function(lpd_point,
           optim_method = "BFGS",
           optim_control = list()){
    
    ## Compute stacking weights based on log-pointwise predictive density ##
    
    stopifnot(is.matrix(lpd_point))
    N <- nrow(lpd_point)
    K <- ncol(lpd_point)
    if (K < 2) {
      stop("At least two models are required for stacking weights.")
    }
    
    lpd_m = mean(lpd_point)
    lpd_point = lpd_point - lpd_m ## rescale the log-density for numerical stability
    
    exp_lpd_point <- exp(lpd_point)
    negative_log_score_loo <- function(w) {
      # objective function: log score
      stopifnot(length(w) == K - 1)
      w_full <- c(w, 1 - sum(w))
      sum <- 0
      for (i in 1:N) {
        sum <- sum + log(exp(lpd_point[i, ]) %*% w_full)
      }
      return(-as.numeric(sum))
    }
    
    gradient <- function(w) {
      # gradient of the objective function
      stopifnot(length(w) == K - 1)
      w_full <- c(w, 1 - sum(w))
      grad <- rep(0, K - 1)
      for (k in 1:(K - 1)) {
        for (i in 1:N) {
          grad[k] <- grad[k] +
            (exp_lpd_point[i, k] - exp_lpd_point[i, K]) / (exp_lpd_point[i,]  %*% w_full)
        }
      }
      return(-grad)
    }
    
    ui <- rbind(rep(-1, K - 1), diag(K - 1))  # K-1 simplex constraint matrix
    ci <- c(-1, rep(0, K - 1))
    w <- constrOptim(
      theta = rep(1 / K, K - 1),
      f = negative_log_score_loo,
      grad = gradient,
      ui = ui,
      ci = ci,
      method = optim_method,
      control = optim_control
    )$par
    
    wts <- structure(
      c(w, 1 - sum(w)),
      names = paste0("model", 1:K),
      class = c("stacking_weights")
    )
    
    return(wts)
  }

Conj_predict <- function(X.mod, y.mod, coords.mod, deltasq_pick, phi_pick, 
                         nu_pick, priors, X.ho, coords.ho){
  
  chol_nk = chol(geoR::matern(as.matrix(dist(coords.mod)), phi = 1/phi_pick,
                              kappa = nu_pick))
  #chol_nk = chol(exp(- phi_pick * as.matrix(dist(coords.mod))))
  N.mod = nrow(coords.mod)
  w_expect = rep(NA, nrow(coords.mod) + nrow(coords.ho))
  # compute the mu_star #
  if(!is.null(X.mod)){
    p = ncol(X.mod)
    chol_inv_V = chol(rbind(cbind(crossprod(X.mod) / deltasq_pick + 
                                    priors$inv_V_beta, 
                                  t(X.mod) / deltasq_pick),
                            cbind(X.mod / deltasq_pick, 
                                  chol2inv(chol_nk) +
                                    diag(rep(1 / deltasq_pick, N.mod)))))
    
    inv_V_mu_beta = priors$inv_V_beta %*% priors$mu_beta
    
    mu_star = c(inv_V_mu_beta + crossprod(X.mod, y.mod) / 
                  deltasq_pick, y.mod / deltasq_pick)
    mu_star = forwardsolve(chol_inv_V, mu_star, transpose = TRUE, 
                           upper.tri = TRUE)
    mu_star = backsolve(chol_inv_V, mu_star)
    
    
    ## compute the expected w on unobserved locations ##
    C_k_nk = geoR::matern(cdist(coords.ho, coords.mod), phi = 1/phi_pick,
                          kappa = nu_pick)
    
    w_U_expect =  forwardsolve(chol_nk, mu_star[-(1:p)], transpose = TRUE, 
                               upper.tri = TRUE)
    w_U_expect =  backsolve(chol_nk, w_U_expect)
    w_U_expect = (C_k_nk %*% w_U_expect)
    y_expect = X.ho %*% mu_star[1:p] + w_U_expect
    w_expect = c(mu_star[-(1:p)], w_U_expect)
  }else{ # p = 0
    chol_inv_V = chol(chol2inv(chol_nk) + diag(rep(1 / deltasq_pick, N.mod)))
    
    mu_star = y.mod / deltasq_pick
    mu_star = forwardsolve(chol_inv_V, mu_star, transpose = TRUE, 
                           upper.tri = TRUE)
    mu_star = backsolve(chol_inv_V, mu_star)
    
    
    ## compute the expected w on unobserved locations ##
    C_k_nk = geoR::matern(cdist(coords.ho, coords.mod), phi = 1/phi_pick,
                          kappa = nu_pick)
    
    w_U_expect =  forwardsolve(chol_nk, mu_star, transpose = TRUE, 
                               upper.tri = TRUE)
    w_U_expect =  backsolve(chol_nk, w_U_expect)
    w_U_expect = (C_k_nk %*% w_U_expect)
    y_expect = w_U_expect
    w_expect = c(mu_star, w_U_expect)
  }
  return(list(y_expect = y_expect, w_expect = w_expect))
}


Conj_lpd_old <- function(X.mod, y.mod, coords.mod, deltasq_pick, phi_pick, nu_pick,
                         priors, X.ho, y.ho, coords.ho, L = 300, seed = 123,
                         MC = TRUE){
  
  set.seed(seed)
  chol_nk = chol(geoR::matern(as.matrix(dist(coords.mod)), phi = 1/phi_pick,
                              kappa = nu_pick))
  #chol_nk = chol(exp(- phi_pick * as.matrix(dist(coords.mod))))
  N.mod = nrow(coords.mod)
  N.ho = nrow(coords.ho)
  lp_expect <- c()
  if(!is.null(X.mod)){
    p = ncol(X.mod)
    # compute the mu_star #
    chol_inv_V = chol(rbind(cbind(crossprod(X.mod) / deltasq_pick + 
                                    priors$inv_V_beta, 
                                  t(X.mod) / deltasq_pick),
                            cbind(X.mod / deltasq_pick, 
                                  chol2inv(chol_nk) +
                                    diag(rep(1 / deltasq_pick, N.mod)))))
    
    inv_V_mu_beta = priors$inv_V_beta %*% priors$mu_beta
    
    mu_star = c(inv_V_mu_beta + crossprod(X.mod, y.mod) / 
                  deltasq_pick, y.mod / deltasq_pick)
    
    mu_star = forwardsolve(chol_inv_V, mu_star, transpose = TRUE, 
                           upper.tri = TRUE)
    
    #mu_star = backsolve(chol_inv_V, mu_star)
    
    
    b_star = priors$b_sigma + 0.5 * (sum(y.mod^2) / deltasq_pick + 
                                       sum(priors$mu_beta * inv_V_mu_beta) - 
                                       sum(mu_star^2))
    
    a_star = priors$a_sigma + N.mod / 2
    
    
    ## generate posterior samples ##
    sigma.sq.sam =  rinvgamma(L, shape = a_star, rate = b_star) 
    # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
    gamma.sam = matrix(rnorm((N.mod + p) * L), nrow = N.mod + p, ncol = L) %*% 
      Diagonal(n = L, sqrt(sigma.sq.sam)) + mu_star
    gamma.sam = backsolve(chol_inv_V, gamma.sam)
    
    
    ## compute the expected w on unobserved locations ##
    C_k_nk = geoR::matern(cdist(coords.ho, coords.mod), phi = 1/phi_pick,
                          kappa = nu_pick)
    
    w_U_expect =  forwardsolve(chol_nk, gamma.sam[-(1:p), ], transpose = TRUE, 
                               upper.tri = TRUE)
    
    w_U_expect =  backsolve(chol_nk, w_U_expect)
    w_U_expect = (C_k_nk %*% w_U_expect)
    y_U_expect = (X.ho %*% gamma.sam[(1:p), ] + w_U_expect)
    
    # the matrix of log ratios (lp)
    M_r <- y_U_expect - y.ho 
    M_r <- M_r %*% Diagonal(n = L, sqrt(1 / (deltasq_pick * sigma.sq.sam)))
    M_r <- -0.5 * (M_r^2 + log(2 * pi) +
                     tcrossprod(rep(1, N.ho), 
                                log(deltasq_pick * sigma.sq.sam)))
    
    #lp_expect <- rowMeans(M_r)
    lp_expect <- log(rowMeans(exp(M_r)))
  }else{ # p = 0
    # compute the mu_star #
    chol_inv_V = chol(chol2inv(chol_nk) + diag(rep(1 / deltasq_pick, N.mod)))
    
    mu_star = c(y.mod / deltasq_pick)
    
    mu_star = forwardsolve(chol_inv_V, mu_star, transpose = TRUE, 
                           upper.tri = TRUE)
    
    b_star = priors$b_sigma + 0.5 * (sum(y.mod^2) / deltasq_pick -
                                       sum(mu_star^2))
    
    a_star = priors$a_sigma + (N.mod) / 2
    
    
    ## generate posterior samples ##
    sigma.sq.sam =  rinvgamma(L, shape = a_star, rate = b_star) 
    # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
    gamma.sam = matrix(rnorm(N.mod * L), nrow = N.mod, ncol = L) %*% 
      Diagonal(n = L, sqrt(sigma.sq.sam)) + mu_star
    gamma.sam = backsolve(chol_inv_V, gamma.sam)
    
    
    ## compute the expected w on unobserved locations ##
    C_k_nk = geoR::matern(cdist(coords.ho, coords.mod), phi = 1/phi_pick,
                          kappa = nu_pick)
    
    w_U_expect =  forwardsolve(chol_nk, gamma.sam, transpose = TRUE, 
                               upper.tri = TRUE)
    w_U_expect =  backsolve(chol_nk, w_U_expect)
    w_U_expect = (C_k_nk %*% w_U_expect)
    y_U_expect = w_U_expect
    
    # the matrix of log ratios (lp)
    M_r <- y_U_expect - y.ho 
    M_r <- M_r %*% Diagonal(n = L, sqrt(1 / (deltasq_pick * sigma.sq.sam)))
    M_r <- -0.5 * (M_r^2 + log(2 * pi) +
                     tcrossprod(rep(1, N.ho), 
                                log(deltasq_pick * sigma.sq.sam)))
    
    lp_expect <- log(rowMeans(exp(M_r)))
  }
  return(lp_expect = lp_expect)
}

Conj_pos_sam <- function(X.mod, y.mod, coords.mod, deltasq_pick, phi_pick, nu_pick,
                         priors, X.ho, y.ho, coords.ho, L = 300, seed = 123,
                         MC = TRUE){
  
  set.seed(seed)
  chol_nk = chol(geoR::matern(as.matrix(dist(coords.mod)), phi = 1/phi_pick,
                              kappa = nu_pick))
  #chol_nk = chol(exp(- phi_pick * as.matrix(dist(coords.mod))))
  N.mod = nrow(coords.mod)
  N.ho = nrow(coords.ho)
  lp_expect <- c()
  if(!is.null(X.mod)){
    p = ncol(X.mod)
    # compute the mu_star #
    chol_inv_V = chol(rbind(cbind(crossprod(X.mod) / deltasq_pick + 
                                    priors$inv_V_beta, 
                                  t(X.mod) / deltasq_pick),
                            cbind(X.mod / deltasq_pick, 
                                  chol2inv(chol_nk) +
                                    diag(rep(1 / deltasq_pick, N.mod)))))
    
    inv_V_mu_beta = priors$inv_V_beta %*% priors$mu_beta
    
    mu_star = c(inv_V_mu_beta + crossprod(X.mod, y.mod) / 
                  deltasq_pick, y.mod / deltasq_pick)
    
    mu_star = forwardsolve(chol_inv_V, mu_star, transpose = TRUE, 
                           upper.tri = TRUE)
    
    #mu_star = backsolve(chol_inv_V, mu_star)
    
    
    b_star = priors$b_sigma + 0.5 * (sum(y.mod^2) / deltasq_pick + 
                                       sum(priors$mu_beta * inv_V_mu_beta) - 
                                       sum(mu_star^2))
    
    a_star = priors$a_sigma + N.mod / 2
    
    
    ## generate posterior samples ##
    sigma.sq.sam =  rinvgamma(L, shape = a_star, rate = b_star) 
    # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
    gamma.sam = matrix(rnorm((N.mod + p) * L), nrow = N.mod + p, ncol = L) %*% 
      Diagonal(n = L, sqrt(sigma.sq.sam)) + mu_star
    gamma.sam = backsolve(chol_inv_V, gamma.sam)
    
    
    ## compute the expected w on unobserved locations ##
    C_k_nk = geoR::matern(cdist(coords.ho, coords.mod), phi = 1/phi_pick,
                          kappa = nu_pick)
    
    w_U_expect =  forwardsolve(chol_nk, gamma.sam[-(1:p), ], transpose = TRUE, 
                               upper.tri = TRUE)
    
    w_U_expect =  backsolve(chol_nk, w_U_expect)
    w_U_expect = (C_k_nk %*% w_U_expect)
    y_U_expect = (X.ho %*% gamma.sam[(1:p), ] + w_U_expect) + 
      matrix(rnorm(N.ho * L), nrow = N.ho, ncol = L) %*% 
      diag(sqrt(sigma.sq.sam * deltasq_pick))
    
  }else{ # p = 0
    # compute the mu_star #
    chol_inv_V = chol(chol2inv(chol_nk) + diag(rep(1 / deltasq_pick, N.mod)))
    
    mu_star = c(y.mod / deltasq_pick)
    
    mu_star = forwardsolve(chol_inv_V, mu_star, transpose = TRUE, 
                           upper.tri = TRUE)
    
    b_star = priors$b_sigma + 0.5 * (sum(y.mod^2) / deltasq_pick -
                                       sum(mu_star^2))
    
    a_star = priors$a_sigma + (N.mod) / 2
    
    
    ## generate posterior samples ##
    sigma.sq.sam =  rinvgamma(L, shape = a_star, rate = b_star) 
    # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
    gamma.sam = matrix(rnorm(N.mod * L), nrow = N.mod, ncol = L) %*% 
      Diagonal(n = L, sqrt(sigma.sq.sam)) + mu_star
    gamma.sam = backsolve(chol_inv_V, gamma.sam)
    
    
    ## compute the expected w on unobserved locations ##
    C_k_nk = geoR::matern(cdist(coords.ho, coords.mod), phi = 1/phi_pick,
                          kappa = nu_pick)
    
    w_U_expect =  forwardsolve(chol_nk, gamma.sam, transpose = TRUE, 
                               upper.tri = TRUE)
    w_U_expect =  backsolve(chol_nk, w_U_expect)
    w_U_expect = (C_k_nk %*% w_U_expect)
    y_U_expect = w_U_expect + 
      matrix(rnorm(N.ho * L), nrow = N.ho, ncol = L) %*% 
      diag(sqrt(sigma.sq.sam * deltasq_pick))
  }
  return(list(sigma.sq.sam = sigma.sq.sam,
              gamma.sam = gamma.sam,
              w_U_expect = w_U_expect,
              y_U_expect = y_U_expect))
}


Conj_lpd <- function(X.mod, y.mod, coords.mod, deltasq_pick, phi_pick, nu_pick,
                     priors, X.ho, y.ho, coords.ho, L = 300, seed = 123,
                     MC = FALSE){
  set.seed(seed)
  chol_nk = chol(geoR::matern(as.matrix(dist(coords.mod)), phi = 1/phi_pick,
                              kappa = nu_pick))
  invR_nk <- chol2inv(chol_nk)
  R_k_nk = geoR::matern(cdist(coords.ho, coords.mod), phi = 1/phi_pick,
                        kappa = nu_pick)
  #chol_nk = chol(exp(- phi_pick * as.matrix(dist(coords.mod))))
  N.mod = nrow(coords.mod)
  N.ho = nrow(coords.ho)
  lp_expect <- c()
  
  if(is.null(X)){p = 0}else{
    p = ncol(X) # number of predictors, including intercept
    inv_V_mu_beta = priors$inv_V_beta %*% priors$mu_beta
  }
  
  if(p == 0){
    chol_inv_V = chol(invR_nk + diag(rep(1 / deltasq_pick, N.mod)))
    mu_star = c(y.mod / deltasq_pick)
  }else{ 
    # compute the mu_star #
    chol_inv_V = chol(rbind(cbind(crossprod(X.mod) / deltasq_pick + 
                                    priors$inv_V_beta, 
                                  t(X.mod) / deltasq_pick),
                            cbind(X.mod / deltasq_pick, 
                                  invR_nk +
                                    diag(rep(1 / deltasq_pick, N.mod)))))
    mu_star = c(inv_V_mu_beta + crossprod(X.mod, y.mod) / 
                  deltasq_pick, y.mod / deltasq_pick)
  }
  
  mu_star = forwardsolve(chol_inv_V, mu_star, transpose = TRUE, 
                         upper.tri = TRUE)
  if(p == 0){
    b_star = priors$b_sigma + 0.5 * (sum(y.mod^2) / deltasq_pick -
                                       sum(mu_star^2))
  }else{
    b_star = priors$b_sigma + 0.5 * (sum(y.mod^2) / deltasq_pick + 
                                       sum(priors$mu_beta * inv_V_mu_beta) - 
                                       sum(mu_star^2))
  }
  a_star = priors$a_sigma + N.mod / 2
  
  if(MC){
    ## use Monte Carlo method to estimate lpd
    ## generate posterior samples ##
    sigma.sq.sam =  rinvgamma(L, shape = a_star, rate = b_star) 
    # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
    gamma.sam = matrix(rnorm((N.mod + p) * L), nrow = N.mod + p, ncol = L) %*% 
      Diagonal(n = L, sqrt(sigma.sq.sam)) + mu_star
    gamma.sam = backsolve(chol_inv_V, gamma.sam)
    
    ## compute the expected w on unobserved locations ##
    if(p == 0){
      w_U_expect =  (R_k_nk %*% invR_nk) %*% gamma.sam
      y_U_expect = w_U_expect
    }else{
      w_U_expect =  (R_k_nk %*% invR_nk) %*% gamma.sam[-(1:p), ]
      y_U_expect = (X.ho %*% gamma.sam[(1:p), ] + w_U_expect)
    }
    
    # the matrix of log ratios (lp)
    M_r <- y_U_expect - y.ho 
    M_r <- M_r %*% Diagonal(n = L, sqrt(1 / (deltasq_pick * sigma.sq.sam)))
    M_r <- -0.5 * (M_r^2 + log(2 * pi) +
                     tcrossprod(rep(1, N.ho), 
                                log(deltasq_pick * sigma.sq.sam)))
    lp_expect <- log(rowMeans(exp(M_r)))
  }else{
    u <- backsolve(chol_inv_V, mu_star) # expected beta and z
    
    lp_c <- -0.5 * log(2 * pi) + lgamma(a_star + 0.5) - 
      lgamma(a_star) + a_star * log(b_star)
    
    if(p == 0){
      H <- tcrossprod(invR_nk, R_k_nk)
      y_U_expect <- c(u %*% H)
      H <- forwardsolve(chol_inv_V, H, transpose = TRUE, 
                        upper.tri = TRUE)
    }else{
      H <- rbind(t(X.ho), tcrossprod(invR_nk, R_k_nk))
      y_U_expect <- c(u %*% H)
      H <- forwardsolve(chol_inv_V, H, transpose = TRUE, 
                        upper.tri = TRUE)
    }
    
    Vs_ls <- apply(H, 2, f <- function(s)(sum(s^2))) + deltasq_pick
    lp_expect <- 
      lp_c - 0.5 * log(Vs_ls) - (a_star + 0.5) * 
      log(b_star + (y.ho - y_U_expect)^2 / (2 * Vs_ls))
  }
  return(lp_expect = lp_expect)
}

expects_MCMC <- function(theta.recover, beta.recover, y.mod, X.mod, coords.mod, 
                         X.ho, y.ho, coords.ho){
  t0 <- proc.time()
  N.mod = nrow(coords.mod)
  N.ho = nrow(coords.ho)
  N.all = N.mod + N.ho
  n.sam = nrow(theta.recover)
  dist.M = as.matrix(dist(rbind(coords.mod, coords.ho)))
  w.recover.sample <- sapply(1:n.sam, f <- function(ind){
    diag_ele = c(rep(1 / theta.recover[ind, "tau.sq"], N.mod), rep(0, N.ho))
    Chol_Cov_w <- chol2inv(chol(
      geoR::matern(dist.M, phi = 1 / theta.recover[ind, "phi"], 
                   kappa = theta.recover[ind, "nu"]))) / 
      theta.recover[ind, "sigma.sq"] 
    
    diag(Chol_Cov_w)  =  diag(Chol_Cov_w) + diag_ele
    Chol_Cov_w <- chol(Chol_Cov_w)
    u = c((y.mod - X.mod %*% beta.recover[ind, ])/theta.recover[ind, "tau.sq"], 
          rep(0, N.ho))
    u <- forwardsolve(Chol_Cov_w, u, upper.tri = TRUE, transpose = TRUE)
    v = rnorm(N.all)
    u <- backsolve(Chol_Cov_w, u + v, upper.tri = TRUE, transpose = FALSE)
    return(u)
  })
  
  w_expect_MCMC <- rowMeans(w.recover.sample)
  
  y.ho.pred.sample <- tcrossprod(X.ho, beta.recover) + 
    w.recover.sample[(N.mod+1):N.all, ] 
  
  lp.ho.recover <- -0.5 * ((y.ho - y.ho.pred.sample)^2 %*% 
                             diag(1 / theta.recover[, "tau.sq"]) + log(2 * pi) + 
                             tcrossprod(rep(1, N.ho), log(theta.recover[, "tau.sq"])))
  
  y.ho.sample <- y.ho.pred.sample + 
    matrix(rnorm(N.ho * n.sam), nrow = N.ho) %*% 
    diag(sqrt(theta.recover[, "tau.sq"]))
  
  
  w_expect_MCMC <- rowMeans(w.recover.sample)
  y_expect_MCMC <- rowMeans(y.ho.sample)
  #lp_expect_MCMC <- rowMeans(lp.ho.recover)
  lp_expect_MCMC <- log(rowMeans(exp(lp.ho.recover)))
  t1 <- proc.time() - t0
  return(list(w_expect_MCMC = w_expect_MCMC,
              y_expect_MCMC = y_expect_MCMC,
              lp_expect_MCMC = lp_expect_MCMC,
              time = t1))
  
}

sp_stacking_K_fold <- function(X, y, coords, deltasq_grid, phi_grid, nu_grid,
                               priors, K_fold = 10, seed = 1, J = 200,
                               label = "LSE", MC = FALSE){
  
  # J: number of samples for computing log point-wise predictive density label LP #
  # K: number of folds
  
  # generate all candidates of the fixed parameter
  t <- proc.time()
  grid_all <- expand.grid(deltasq_grid, phi_grid, nu_grid)
  colnames(grid_all) <- c("deltasq", "phi", "nu")
  rownames(grid_all) <- paste0("model", 1:nrow(grid_all))
  grid_phi_nv = expand.grid(phi_grid, nu_grid)
  colnames(grid_phi_nv) <- c("phi", "nu")
  L_grid_deltasq =length(deltasq_grid)
  
  # # pre-computation
  N = length(y)
  if(is.null(X)){p = 0}else{
    p = ncol(X) # number of predictors, including intercept
    inv_V_mu_beta = priors$inv_V_beta %*% priors$mu_beta
  }
  
  ## CV-folds ##
  set.seed(seed)
  folds <- cvFolds(N, K_fold, type = "random")
  # pre-computation and pre-allocation#
  ind_k_list <- lapply(1:K_fold, function(x){c(folds$subsets)[(folds$which == x)]})
  nk_list <- sapply(ind_k_list, length)
  if(p > 0){
    XTX_list <- lapply(1:K_fold, function(x){crossprod(X[-ind_k_list[[x]], ])})
    XTy_list <- lapply(1:K_fold, function(x){crossprod(X[-ind_k_list[[x]], ],
                                                       y[-ind_k_list[[x]]]) })
  }
  if(label == "LSE"){
    y_expect <- matrix(NA, N, nrow(grid_all))
  }else if(label == "LP"){
    lp_expect <- matrix(NA, N, nrow(grid_all))
    y_sq_sum_list <- sapply(1:K_fold, function(x){sum(y[-ind_k_list[[x]]]^2)})
  }else{print("label has to be LSE or LP")}
  
  for (i1 in 1:nrow(grid_phi_nv)){
    phi_pick <- grid_phi_nv$phi[i1]
    nu_pick <- grid_phi_nv$nu[i1]
    for (k in 1:K_fold){
      invR_nk = chol2inv(chol(geoR::matern(as.matrix(dist(coords[-ind_k_list[[k]], ])),
                                           phi = 1 / phi_pick, kappa = nu_pick)))
      R_k_nk = geoR::matern(cdist(coords[ind_k_list[[k]], ],
                                  coords[-ind_k_list[[k]], ]),
                            phi = 1 / phi_pick, kappa = nu_pick)
      for (i2 in 1:L_grid_deltasq){ # should use parallel computing
        deltasq_pick <- deltasq_grid[i2]
        # Compute Cholesky decomposition of M_*^{-1}
        if(p == 0){
          chol_inv_M = chol(invR_nk + diag(rep(1 / deltasq_pick, N - nk_list[k])))
          u = y[-ind_k_list[[k]]] / deltasq_pick
        }else{
          chol_inv_M = chol(
            rbind(cbind(XTX_list[[k]] / deltasq_pick + priors$inv_V_beta,
                        t(X[-ind_k_list[[k]], ]) / deltasq_pick),
                  cbind(X[-ind_k_list[[k]], ] / deltasq_pick,
                        invR_nk + diag(rep(1 / deltasq_pick, N - nk_list[k])))))
          u = c(inv_V_mu_beta + XTy_list[[k]] / deltasq_pick,
                y[-ind_k_list[[k]]] / deltasq_pick)
          
        }
        u <- forwardsolve(chol_inv_M, u, transpose = TRUE, upper.tri = TRUE)
        
        if (label == "LSE"){
          ## Stacking based on expectation
          # compute expectation of response in fold k
          u <- backsolve(chol_inv_M, u)
          if(p == 0){
            w_U_expect <- R_k_nk %*% (invR_nk %*% u)
            y_expect[ind_k_list[[k]], (i1-1)*L_grid_deltasq + i2] <- w_U_expect
          }else{
            w_U_expect <- R_k_nk %*% (invR_nk %*% u[-(1:p)])
            y_expect[ind_k_list[[k]], (i1-1)*L_grid_deltasq + i2] <-
              X[ind_k_list[[k]], ] %*% u[(1:p)] + w_U_expect
          }
        }else{
          ## Stacking based on log point-wise predictive density
          
          if(p == 0){
            b_star = priors$b_sigma +
              0.5 * (y_sq_sum_list[k] / deltasq_pick - sum(u^2))
          }else{
            b_star = priors$b_sigma +
              0.5 * (y_sq_sum_list[k] / deltasq_pick +
                       sum(priors$mu_beta * inv_V_mu_beta) - sum(u^2))
          }
          
          a_star = priors$a_sigma + (N - nk_list[k]) / 2
          
          if(MC){
            ## use Monte Carlo method to estimate lpd
            ## generate posterior samples ##
            sigma.sq.sam =  rinvgamma(J, shape = a_star, rate = b_star)
            # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
            
            ## compute the expected response on unobserved locations ##
            gamma.sam = matrix(rnorm((N - nk_list[k] + p) * J),
                               nrow = N - nk_list[k] + p, ncol = J) %*%
              Diagonal(n = J, sqrt(sigma.sq.sam)) + u
            gamma.sam = backsolve(chol_inv_M, gamma.sam)
            
            if(p == 0){
              w_U_expect = (R_k_nk %*% invR_nk) %*% gamma.sam
              y_U_expect = w_U_expect
            }else{
              w_U_expect = (R_k_nk %*% invR_nk) %*% gamma.sam[-(1:p), ]
              y_U_expect = (X[ind_k_list[[k]], ] %*% gamma.sam[(1:p), ] +
                              w_U_expect)
            }
            # the matrix of log ratios (lp)
            M_r <- y_U_expect - y[ind_k_list[[k]]]
            M_r <- M_r %*% Diagonal(n = J, sqrt(1 / (deltasq_pick * sigma.sq.sam)))
            M_r <- -0.5 * (M_r^2 + log(2 * pi) +
                             tcrossprod(rep(1, nk_list[k]),
                                        log(deltasq_pick * sigma.sq.sam)))
            lp_expect[ind_k_list[[k]], (i1-1)*L_grid_deltasq + i2] <-
              log(rowMeans(exp(M_r))) 
          }else{
            u <- backsolve(chol_inv_M, u) # expected beta and z
            
            lp_c <- -0.5 * log(2 * pi) + lgamma(a_star + 0.5) - 
              lgamma(a_star) + a_star * log(b_star)
            
            if(p == 0){
              H <- tcrossprod(invR_nk, R_k_nk)
              y_U_expect <- c(u %*% H)
              H <- forwardsolve(chol_inv_M, H, transpose = TRUE, 
                                upper.tri = TRUE)
            }else{
              H <- rbind(t(X[ind_k_list[[k]], ]), tcrossprod(invR_nk, R_k_nk))
              y_U_expect <- c(u %*% H)
              H <- forwardsolve(chol_inv_M, H, transpose = TRUE, 
                                upper.tri = TRUE)
            }
            
            Vs_ls <- apply(H, 2, f <- function(s)(sum(s^2))) + deltasq_pick
            lp_expect[ind_k_list[[k]], (i1-1)*L_grid_deltasq + i2] <- 
              lp_c - 0.5 * log(Vs_ls) - (a_star + 0.5) * 
              log(b_star + (y[ind_k_list[[k]]] - y_U_expect)^2 / (2 * Vs_ls))
          }
        }
      }
    }
  }
  
  ## Compute stacking weights ##
  if(label == "LSE"){
    wts <- QP_stacking_weight(y_expect, y)
    
    time <- proc.time()-t
    return(list(wts = wts,
                grid_all = grid_all,
                time = time
    ))
  }else{
    time <- proc.time()-t
    wts <- stacking_weights(lp_expect)
    return(list(wts = wts,
                grid_all = grid_all,
                time = time
    ))
  }
}


#### test for different algorithms for prefixs ####
sp_stacking_K_fold2 <- function(X, y, coords, deltasq_grid, phi_nu_ls,
                                priors, K_fold = 10, seed = 1, J = 200,
                                label = "LSE", MC = FALSE){
  
  # J: number of samples for computing log point-wise predictive density label LP #
  # K: number of folds
  # phi_nu_ls: a matrix of two columns recording candidate values for phi & nu.
  
  # generate all candidates of the fixed parameter
  t <- proc.time()
  grid_all <- expand.grid(deltasq_grid, phi_nu_ls[, "phi"])
  grid_all[, 3] <- rep(phi_nu_ls[, "nu"], each = length(deltasq_grid))
  colnames(grid_all) <- c("deltasq", "phi", "nu")
  rownames(grid_all) <- paste0("model", 1:nrow(grid_all))
  grid_phi_nv = as.data.frame(phi_nu_ls)
  colnames(grid_phi_nv) <- c("phi", "nu")
  L_grid_deltasq =length(deltasq_grid)
  
  # # pre-computation
  N = length(y)
  if(is.null(X)){p = 0}else{
    p = ncol(X) # number of predictors, including intercept
    inv_V_mu_beta = priors$inv_V_beta %*% priors$mu_beta
  }
  
  ## CV-folds ##
  set.seed(seed)
  folds <- cvFolds(N, K_fold, type = "random")
  # pre-computation and pre-allocation#
  ind_k_list <- lapply(1:K_fold, function(x){c(folds$subsets)[(folds$which == x)]})
  nk_list <- sapply(ind_k_list, length)
  if(p > 0){
    XTX_list <- lapply(1:K_fold, function(x){crossprod(X[-ind_k_list[[x]], ])})
    XTy_list <- lapply(1:K_fold, function(x){crossprod(X[-ind_k_list[[x]], ],
                                                       y[-ind_k_list[[x]]]) })
  }
  if(label == "LSE"){
    y_expect <- matrix(NA, N, nrow(grid_all))
  }else if(label == "LP"){
    lp_expect <- matrix(NA, N, nrow(grid_all))
    y_sq_sum_list <- sapply(1:K_fold, function(x){sum(y[-ind_k_list[[x]]]^2)})
  }else{print("label has to be LSE or LP")}
  
  for (i1 in 1:nrow(grid_phi_nv)){
    phi_pick <- grid_phi_nv$phi[i1]
    nu_pick <- grid_phi_nv$nu[i1]
    for (k in 1:K_fold){
      t <- proc.time()
      invR_nk = chol2inv(chol(geoR::matern(as.matrix(dist(coords[-ind_k_list[[k]], ])),
                                           phi = 1 / phi_pick, kappa = nu_pick)))
      proc.time() - t
      R_k_nk = geoR::matern(cdist(coords[ind_k_list[[k]], ],
                                  coords[-ind_k_list[[k]], ]),
                            phi = 1 / phi_pick, kappa = nu_pick)
      for (i2 in 1:L_grid_deltasq){ # should use parallel computing
        deltasq_pick <- deltasq_grid[i2]
        # Compute Cholesky decomposition of M_*^{-1}
        if(p == 0){
          chol_inv_M = chol(invR_nk + diag(rep(1 / deltasq_pick, N - nk_list[k])))
          u = y[-ind_k_list[[k]]] / deltasq_pick
        }else{
          chol_inv_M = chol(
            rbind(cbind(XTX_list[[k]] / deltasq_pick + priors$inv_V_beta,
                        t(X[-ind_k_list[[k]], ]) / deltasq_pick),
                  cbind(X[-ind_k_list[[k]], ] / deltasq_pick,
                        invR_nk + diag(rep(1 / deltasq_pick, N - nk_list[k])))))
          u = c(inv_V_mu_beta + XTy_list[[k]] / deltasq_pick,
                y[-ind_k_list[[k]]] / deltasq_pick)
          
        }
        u <- forwardsolve(chol_inv_M, u, transpose = TRUE, upper.tri = TRUE)
        
        if (label == "LSE"){
          ## Stacking based on expectation
          # compute expectation of response in fold k
          u <- backsolve(chol_inv_M, u)
          if(p == 0){
            w_U_expect <- R_k_nk %*% (invR_nk %*% u)
            y_expect[ind_k_list[[k]], (i1-1)*L_grid_deltasq + i2] <- w_U_expect
          }else{
            w_U_expect <- R_k_nk %*% (invR_nk %*% u[-(1:p)])
            y_expect[ind_k_list[[k]], (i1-1)*L_grid_deltasq + i2] <-
              X[ind_k_list[[k]], ] %*% u[(1:p)] + w_U_expect
          }
        }else{
          ## Stacking based on log point-wise predictive density
          
          if(p == 0){
            b_star = priors$b_sigma +
              0.5 * (y_sq_sum_list[k] / deltasq_pick - sum(u^2))
          }else{
            b_star = priors$b_sigma +
              0.5 * (y_sq_sum_list[k] / deltasq_pick +
                       sum(priors$mu_beta * inv_V_mu_beta) - sum(u^2))
          }
          
          a_star = priors$a_sigma + (N - nk_list[k]) / 2
          
          if(MC){
            ## use Monte Carlo method to estimate lpd
            ## generate posterior samples ##
            sigma.sq.sam =  rinvgamma(J, shape = a_star, rate = b_star)
            # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
            
            ## compute the expected response on unobserved locations ##
            gamma.sam = matrix(rnorm((N - nk_list[k] + p) * J),
                               nrow = N - nk_list[k] + p, ncol = J) %*%
              Diagonal(n = J, sqrt(sigma.sq.sam)) + u
            gamma.sam = backsolve(chol_inv_M, gamma.sam)
            
            if(p == 0){
              w_U_expect = (R_k_nk %*% invR_nk) %*% gamma.sam
              y_U_expect = w_U_expect
            }else{
              w_U_expect = (R_k_nk %*% invR_nk) %*% gamma.sam[-(1:p), ]
              y_U_expect = (X[ind_k_list[[k]], ] %*% gamma.sam[(1:p), ] +
                              w_U_expect)
            }
            # the matrix of log ratios (lp)
            M_r <- y_U_expect - y[ind_k_list[[k]]]
            M_r <- M_r %*% Diagonal(n = J, sqrt(1 / (deltasq_pick * sigma.sq.sam)))
            M_r <- -0.5 * (M_r^2 + log(2 * pi) +
                             tcrossprod(rep(1, nk_list[k]),
                                        log(deltasq_pick * sigma.sq.sam)))
            lp_expect[ind_k_list[[k]], (i1-1)*L_grid_deltasq + i2] <-
              log(rowMeans(exp(M_r))) 
          }else{
            u <- backsolve(chol_inv_M, u) # expected beta and z
            
            lp_c <- -0.5 * log(2 * pi) + lgamma(a_star + 0.5) - 
              lgamma(a_star) + a_star * log(b_star)
            
            if(p == 0){
              H <- tcrossprod(invR_nk, R_k_nk)
              y_U_expect <- c(u %*% H)
              H <- forwardsolve(chol_inv_M, H, transpose = TRUE, 
                                upper.tri = TRUE)
            }else{
              H <- rbind(t(X[ind_k_list[[k]], ]), tcrossprod(invR_nk, R_k_nk))
              y_U_expect <- c(u %*% H)
              H <- forwardsolve(chol_inv_M, H, transpose = TRUE, 
                                upper.tri = TRUE)
            }
            
            Vs_ls <- apply(H, 2, f <- function(s)(sum(s^2))) + deltasq_pick
            lp_expect[ind_k_list[[k]], (i1-1)*L_grid_deltasq + i2] <- 
              lp_c - 0.5 * log(Vs_ls) - (a_star + 0.5) * 
              log(b_star + (y[ind_k_list[[k]]] - y_U_expect)^2 / (2 * Vs_ls))
          }
        }
      }
    }
  }
  
  ## Compute stacking weights ##
  if(label == "LSE"){
    wts <- QP_stacking_weight(y_expect, y)
    
    time <- proc.time()-t
    return(list(wts = wts,
                grid_all = grid_all,
                time = time
    ))
  }else{
    time <- proc.time()-t
    wts <- stacking_weights(lp_expect)
    return(list(wts = wts,
                grid_all = grid_all,
                time = time
    ))
  }
}

#### test for marginal posterior for prefixs ####
sp_stacking_K_fold3 <- function(X, y, coords, all_prefix_ls,
                                priors, K_fold = 10, seed = 1, J = 200,
                                label = "LSE", MC = FALSE){
  
  # J: number of samples for computing log point-wise predictive density label LP #
  # K: number of folds
  # all_prefix_ls: a data.frame of three variables recording candidate values for deltasq, phi & nu.
  
  # generate all candidates of the fixed parameter
  t <- proc.time()
  grid_all <- as.data.frame(all_prefix_ls)
  colnames(grid_all) <- c("deltasq", "phi", "nu")
  rownames(grid_all) <- paste0("model", 1:nrow(grid_all))
  
  # # pre-computation
  N = length(y)
  if(is.null(X)){p = 0}else{
    p = ncol(X) # number of predictors, including intercept
    inv_V_mu_beta = priors$inv_V_beta %*% priors$mu_beta
  }
  
  ## CV-folds ##
  set.seed(seed)
  folds <- cvFolds(N, K_fold, type = "random")
  # pre-computation and pre-allocation#
  ind_k_list <- lapply(1:K_fold, function(x){c(folds$subsets)[(folds$which == x)]})
  nk_list <- sapply(ind_k_list, length)
  if(p > 0){
    XTX_list <- lapply(1:K_fold, function(x){crossprod(X[-ind_k_list[[x]], ])})
    XTy_list <- lapply(1:K_fold, function(x){crossprod(X[-ind_k_list[[x]], ],
                                                       y[-ind_k_list[[x]]]) })
  }
  if(label == "LSE"){
    y_expect <- matrix(NA, N, nrow(grid_all))
  }else if(label == "LP"){
    lp_expect <- matrix(NA, N, nrow(grid_all))
    y_sq_sum_list <- sapply(1:K_fold, function(x){sum(y[-ind_k_list[[x]]]^2)})
  }else{print("label has to be LSE or LP")}
  
  for (i1 in 1:nrow(grid_all)){ # should use parallel computing
    phi_pick <- all_prefix_ls$phi[i1]
    nu_pick <- all_prefix_ls$nu[i1]
    deltasq_pick <- all_prefix_ls$delatsq[i1]
    for (k in 1:K_fold){
      invR_nk = chol2inv(chol(geoR::matern(as.matrix(dist(coords[-ind_k_list[[k]], ])),
                                           phi = 1 / phi_pick, kappa = nu_pick)))
      R_k_nk = geoR::matern(cdist(coords[ind_k_list[[k]], ],
                                  coords[-ind_k_list[[k]], ]),
                            phi = 1 / phi_pick, kappa = nu_pick)
      # Compute Cholesky decomposition of M_*^{-1}
      if(p == 0){
        chol_inv_M = chol(invR_nk + diag(rep(1 / deltasq_pick, N - nk_list[k])))
        u = y[-ind_k_list[[k]]] / deltasq_pick
      }else{
        chol_inv_M = chol(
          rbind(cbind(XTX_list[[k]] / deltasq_pick + priors$inv_V_beta,
                      t(X[-ind_k_list[[k]], ]) / deltasq_pick),
                cbind(X[-ind_k_list[[k]], ] / deltasq_pick,
                      invR_nk + diag(rep(1 / deltasq_pick, N - nk_list[k])))))
        u = c(inv_V_mu_beta + XTy_list[[k]] / deltasq_pick,
              y[-ind_k_list[[k]]] / deltasq_pick)
        
      }
      u <- forwardsolve(chol_inv_M, u, transpose = TRUE, upper.tri = TRUE)
      
      if (label == "LSE"){
        ## Stacking based on expectation
        # compute expectation of response in fold k
        u <- backsolve(chol_inv_M, u)
        if(p == 0){
          w_U_expect <- R_k_nk %*% (invR_nk %*% u)
          y_expect[ind_k_list[[k]], i1] <- w_U_expect
        }else{
          w_U_expect <- R_k_nk %*% (invR_nk %*% u[-(1:p)])
          y_expect[ind_k_list[[k]], i1] <-
            X[ind_k_list[[k]], ] %*% u[(1:p)] + w_U_expect
        }
      }else{
        ## Stacking based on log point-wise predictive density
        
        if(p == 0){
          b_star = priors$b_sigma +
            0.5 * (y_sq_sum_list[k] / deltasq_pick - sum(u^2))
        }else{
          b_star = priors$b_sigma +
            0.5 * (y_sq_sum_list[k] / deltasq_pick +
                     sum(priors$mu_beta * inv_V_mu_beta) - sum(u^2))
        }
        
        a_star = priors$a_sigma + (N - nk_list[k]) / 2
        
        if(MC){
          ## use Monte Carlo method to estimate lpd
          ## generate posterior samples ##
          sigma.sq.sam =  rinvgamma(J, shape = a_star, rate = b_star)
          # the rate in this function is the scale in wikipedia, mean = b_star/(a_star-1)
          
          ## compute the expected response on unobserved locations ##
          gamma.sam = matrix(rnorm((N - nk_list[k] + p) * J),
                             nrow = N - nk_list[k] + p, ncol = J) %*%
            Diagonal(n = J, sqrt(sigma.sq.sam)) + u
          gamma.sam = backsolve(chol_inv_M, gamma.sam)
          
          if(p == 0){
            w_U_expect = (R_k_nk %*% invR_nk) %*% gamma.sam
            y_U_expect = w_U_expect
          }else{
            w_U_expect = (R_k_nk %*% invR_nk) %*% gamma.sam[-(1:p), ]
            y_U_expect = (X[ind_k_list[[k]], ] %*% gamma.sam[(1:p), ] +
                            w_U_expect)
          }
          # the matrix of log ratios (lp)
          M_r <- y_U_expect - y[ind_k_list[[k]]]
          M_r <- M_r %*% Diagonal(n = J, sqrt(1 / (deltasq_pick * sigma.sq.sam)))
          M_r <- -0.5 * (M_r^2 + log(2 * pi) +
                           tcrossprod(rep(1, nk_list[k]),
                                      log(deltasq_pick * sigma.sq.sam)))
          lp_expect[ind_k_list[[k]], i1] <-
            log(rowMeans(exp(M_r))) 
        }else{
          u <- backsolve(chol_inv_M, u) # expected beta and z
          
          lp_c <- -0.5 * log(2 * pi) + lgamma(a_star + 0.5) - 
            lgamma(a_star) + a_star * log(b_star)
          
          if(p == 0){
            H <- tcrossprod(invR_nk, R_k_nk)
            y_U_expect <- c(u %*% H)
            H <- forwardsolve(chol_inv_M, H, transpose = TRUE, 
                              upper.tri = TRUE)
          }else{
            H <- rbind(t(X[ind_k_list[[k]], ]), tcrossprod(invR_nk, R_k_nk))
            y_U_expect <- c(u %*% H)
            H <- forwardsolve(chol_inv_M, H, transpose = TRUE, 
                              upper.tri = TRUE)
          }
          
          Vs_ls <- apply(H, 2, f <- function(s)(sum(s^2))) + deltasq_pick
          lp_expect[ind_k_list[[k]], i1] <- 
            lp_c - 0.5 * log(Vs_ls) - (a_star + 0.5) * 
            log(b_star + (y[ind_k_list[[k]]] - y_U_expect)^2 / (2 * Vs_ls))
        }
      }
    }
  }
  
  ## Compute stacking weights ##
  if(label == "LSE"){
    wts <- QP_stacking_weight(y_expect, y)
    
    time <- proc.time()-t
    return(list(wts = wts,
                grid_all = grid_all,
                time = time
    ))
  }else{
    time <- proc.time()-t
    wts <- stacking_weights(lp_expect)
    return(list(wts = wts,
                grid_all = grid_all,
                time = time
    ))
  }
}


# test function for calculating the candidate value of prefixed parameters
decay_est <-  function(eff_r_ls, nu_ls){
  # eff_r_ls: list of effective range
  # nu_ls: list of candidate smoothness in matern
  decay_ls <- matrix(NA, nrow = length(eff_r_ls), ncol = length(nu_ls))
  # decay_ls records candidate decay (row) for each nu (col)
  for(i in 1:length(eff_r_ls)){
    for(j in 1:length(nu_ls)){
      decay_ls[i, j] <- 1/cor.to.par(
        d = eff_r_ls[i], #effective range
        param=list(nu = nu_ls[j]),
        family = "matern",
        cor.target = 0.05
      )
    }
  }
  return(decay_ls)
}



stacking_pos_sample <- function(Stack_fit, L1 = 300, L2 = 900, 
                                X.mod, y.mod, coords.mod, priors,
                                X.ho, coords.ho, seed = 123){
  # L1 sample for each candidate model
  # L2 sample size of the returned posterior
  set.seed(seed)
  t <- proc.time()
  j = 1
  pick_mods <- Stack_fit$grid_all[(Stack_fit$wts>0.00001), ]
  pos_y_U <- c()
  pos_w_U <- c()
  pos_sigmasq <- c()
  cat("No. of models:", length(pick_mods$deltasq), "\n")
  
  for (j in 1:length(pick_mods$deltasq)){
    cat(j, "\t")
    t1 <- proc.time()
    pred_pos_sam <- 
      Conj_pos_sam(X.mod = X.mod, y.mod = y.mod,
                   coords.mod = coords.mod,
                   deltasq_pick = pick_mods$deltasq[j],
                   phi_pick = pick_mods$phi[j], 
                   nu_pick = pick_mods$nu[j], priors = priors,
                   X.ho = X.ho, 
                   coords.ho = coords.ho,
                   L = L1)
    pos_y_U[[j]] <- pred_pos_sam$y_U_expect
    pos_sigmasq[[j]] <- pred_pos_sam$sigma.sq.sam
    cat("use time: ", (proc.time() - t1)[3], "\n")
  }
  proc.time() - t
  stack_prob <- Stack_fit$wts[(Stack_fit$wts>0.00001)]
  num_counts <- c(rmultinom(n = 1, size = L2, prob= stack_prob))
  pick_ind <- lapply(1:length(stack_prob), 
                     function(x){sort(sample.int(L, num_counts[x], 
                                                 replace = TRUE))})
  #save the predictive samples for y
  if(length(pick_mods$deltasq) == 1){
    pred_y_U_stack_sam <- pos_y_U[[1]][, pick_ind[[1]]]
  }else{
    pred_y_U_stack_sam <- do.call(cbind,
                                  sapply(1:length(stack_prob), function(x){
                                    pos_y_U[[x]][, pick_ind[[x]]]
                                  }))
  }
  sigmasq_sam <- unlist(sapply(1:length(stack_prob), function(x){
    pos_sigmasq[[x]][pick_ind[[x]]]
  }))
  return(list(sigmasq_sam = sigmasq_sam, 
              pred_y_U_stack_sam = pred_y_U_stack_sam))
}


