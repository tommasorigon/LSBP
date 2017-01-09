library(BayesLogit)
library(Matrix)
library(mvtnorm)


sb <- function(nu) {
  nu    <- c(nu,1)
  H     <- length(nu)
  prob  <- nu * c(1,cumprod(1 - nu[-H]))
  return(prob)
}

#' Bayesian nonparametric density estimation
#'
#' The dependent logistic stick-breaking process (DLSBP) model posterior mod, estimated trhough the EM algorithm
#' @param y A vector containing the response vector
#' @param X A n \times p design matrix containing the covariates
#' @param H The number of mixtures
#' @param maxiter Logical: should be putted a prior on alpha?
#' @param prior Logical: should be putted a prior on alpha?
#' @param verbose Logical: should be putted a prior on alpha?
#' @param tol Logical: should be putted a prior on alpha?
#' @param random Logical: should be putted a prior on alpha?
#' @param pstep Logical: should be putted a prior on alpha?
#' @export
#' 
DLSBP_EM <- function(y, X, H = 2, maxiter = 1000, prior = NULL, verbose = TRUE, tol=1e-4, random=FALSE, pstep=10) {
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  X <- Matrix(X)
  
  # In prior hyperparameters not defined, use the MLE, defined as follow
  if (is.null(prior)) {
    prior <- list(b = rep(0, p), 
                  B = Diagonal(p, Inf), 
                  mu_mu =  0, 
                  tau_mu = 0, 
                  a_tau = 1, 
                  b_tau = 0)
  }
  
  # Hyperparameters
  b <- prior$b
  mu_mu  <- prior$mu_mu
  B <- prior$B; P <- solve(B); Pb <- P %*% b
  tau_mu <- prior$tau_mu
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Initialization
  tau     <- rep(1/diff(quantile(y,c(.25,0.75))), H)
  mu      <- rep(0,H)
  beta    <- Matrix(1e-2, H - 1, p) # A little jitter is added
  logpost <- -Inf
  
  if(random){
    #tau   <- rgamma(H,a_tau,b_tau)    # Initialized a robust estimate of the global sample variance
    mu    <- rnorm(H,mean(y),sd(y))
    beta  <- Matrix(rnorm(p*(H-1),0,.1), H-1,p)
  }
  
  # EM Algorithm
  for (r in 1:maxiter) {
    # Step 2,3,4: performed within the cluster.
    for (h in 1:H) {
      
      # (Expectation) Step 1 - Cluster allocation
      Expectation <- Expectation_step(y, as.matrix(X), as.matrix(beta), mu, tau)
      z <- Expectation$z
      
      # (Expectation) Step 2 - PolyaGamma
      if (h < H) {
        linpred <- as.numeric(X%*%beta[h, ])
        
        # Fixed quantities
        z_sum <- rowSums(z[, h:H])
        omega <- z_sum/(2 * linpred) * tanh(linpred/2)
        is_omega_nan <- is.nan(omega); omega[is_omega_nan] <- z_sum[is_omega_nan]/4
        #diag(Omega) <- omega
        
        # (Maximization) Step 3
        Sigma_beta1 <- crossprod(X*sqrt(omega)) + P#t(X) %*% Omega %*% X + P
        beta[h, ] <- solve(Sigma_beta1, crossprod(X ,z[, h] - z_sum/2) + Pb)
      }

      # (Maximization) Step 4
      tau_tilde <- tau[h]*sum(z[,h]) + tau_mu
      mu[h]     <- (tau[h]*sum(z[,h]* y) + tau_mu*mu_mu)/tau_tilde
      
      residual <- as.numeric(y - mu[h])
      tau[h] <- max(0,(a_tau - 1 + sum(z[, h])/2)/(b_tau + sum(z[, h] * residual^2)/2))
    }
    
    # Convergence checks
    loglik         <- Expectation$loglik
    logpriors      <- sum(dmvnorm(as.matrix(beta), mean = b,sigma = as.matrix(B),log=TRUE)) + sum(dnorm(mu,  mu_mu,1/sqrt(tau_mu),log=TRUE)) + sum(dgamma(tau, a_tau,b_tau,log=TRUE))
    logpost_new    <- loglik + logpriors
    if(!is.finite(logpriors)) logpost_new <- loglik
    
    # Break the loop at convergence
    if((logpost_new - logpost)<tol)break
    
    # Otherwise continue!
    logpost <- logpost_new
    
    # Display status
    if (verbose) {
      if(r%%pstep==0) cat(paste("log-posterior:",round(logpost,15),", iteration:", r, "\n",sep=""))
    }
  }
  
  # Output
  cluster <- apply(z,1,which.max)
  pred    <- prediction(y, as.matrix(X), as.matrix(beta), mu, tau)
  list(beta = beta, mu = mu, tau=tau, pred = pred, cluster=cluster,z=z, logposterior=logpost)
}

#' Bayesian nonparametric density estimation
#'
#' The dependent logistic stick-breaking process (DLSBP) model posterior mod, estimated trhough the EM algorithm
#' @param y A vector containing the response vector
#' @param X A n \times p design matrix containing the covariates
#' @param H The number of mixtures
#' @param maxiter Logical: should be putted a prior on alpha?
#' @param prior Logical: should be putted a prior on alpha?
#' @param verbose Logical: should be putted a prior on alpha?
#' @param tol Logical: should be putted a prior on alpha?
#' @param random Logical: should be putted a prior on alpha?
#' @param pstep Logical: should be putted a prior on alpha?
#' @export
#' 
DLSBP_Gibbs <- function(y, X, H = 2, R = 5000, burn.in = 2000, prior = NULL, verbose = TRUE, parameters=TRUE, useEM=TRUE, pstep=100){
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  X <- Matrix(X)
  
  # In prior hyperparameters not defined, use the MLE, defined as follow
  if (is.null(prior)) {
    prior <- list(b = rep(0, p), 
                  B = Diagonal(p, Inf), 
                  mu_mu =  0, 
                  tau_mu = 0, 
                  a_tau = 1, 
                  b_tau = 0)
  }
  
  # Hyperparameters
  b <- prior$b
  mu_mu  <- prior$mu_mu
  B <- prior$B; P <- solve(B); Pb <- P %*% b
  tau_mu <- prior$tau_mu
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Output 
  if(parameters) {
    beta_out   <- array(0, c(R, H - 1, p)) 
    mu_out     <- matrix(0, R, H)
    tau_out    <- matrix(0, R, H)
  }
  prediction <- matrix(0, R, n)
  
  # Initialization
  
  tau  <- rep(1/diff(quantile(y,c(.25,0.75))), H)
  mu   <- rep(0, H)
  beta <- Matrix(0, H - 1, p)
  
  if(useEM) {
    EM    <- DLSBP_EM(y, X, H=H, prior=prior, maxiter=500,tol=1e-4)
    tau   <- EM$tau
    mu    <- EM$mu
    beta  <- EM$beta
  }
  # Gibbs sampling
  
  cat("Starting the Gibbs sampler...","\n")
  for (r in 1:(R + burn.in)) {
    
    # Step 1 - Cluster allocation
    Update <- G_update(y, as.matrix(X), as.matrix(beta), mu, tau)
    G      <- c(Update$G)
    G_mixt <- c(Update$G_mixt)
    
    # Step 2 - 3: performed within the cluster.
    for (h in 1:H) {
      if (h < H) {
        # Subsetting observations
        index <-G > h - 1
        nh <- sum(index)
        zh <- G[index] == h  
        Xh <- Matrix(matrix(X[index, ], nh, p))
        
        # Step 2 - Polya-gamma weights, beta coefficients
        
        linh <- as.numeric(Xh %*% beta[h, ])
        omega <- rpg.devroye(num = nh, n = 1, z = linh)
        #Omega <- Diagonal(nh)
        #diag(Omega) <- omega
        
        Sigma_beta <- solve(crossprod(Xh*sqrt(omega)) + P) #solve(t(Xh) %*% Omega %*% Xh + P)
        mu_beta    <- Sigma_beta %*% (crossprod(Xh, zh - 1/2) + Pb)
        beta[h, ]  <- c(rmvnorm(1, mean = mu_beta, sigma = as.matrix(Sigma_beta)))
      }
      
      # Step 3 - Mixture components
      indexG <- G == h
      nG <- sum(indexG)
      yG <- y[indexG]
      
      tau_tilde <- tau_mu + tau[h]*nG
      mu_tilde  <- (tau_mu*mu_mu + tau[h]*sum(yG))/tau_tilde
      mu[h]     <- rnorm(1,mu_tilde,1/sqrt(tau_tilde))
      
      residual  <- as.numeric(yG - mu[h])
      tau[h]    <- rgamma(1, a_tau + nG/2, b_tau + sum(residual^2)/2)
    }
    
    # Output
    if (r > burn.in) {
      if(parameters){ 
        beta_out[r - burn.in, , ] <- as.matrix(beta)
        mu_out[r - burn.in, ]     <- mu
        tau_out[r - burn.in, ]    <- tau
      }
      # Prediction
      prediction[r - burn.in, ] <- rnorm(n, mu[G_mixt], as.numeric(1/sqrt(tau[G_mixt])))
    }
    
    if (verbose) {
      if(r%%pstep==0) cat(paste("Sampling iteration:", r, "\n",sep=""))
    }
  }
  if(parameters) return(list(pred = prediction,beta=beta_out,mu=mu_out,tau=tau_out))
  list(pred=prediction)
}

#' Bayesian nonparametric density estimation
#'
#' The dependent logistic stick-breaking process (DLSBP) model posterior mod, estimated trhough the EM algorithm
#' @param y A vector containing the response vector
#' @param X A n x p design matrix containing the covariates
#' @param H The number of mixtures
#' @param maxiter Logical: should be putted a prior on alpha?
#' @param prior Logical: should be putted a prior on alpha?
#' @param verbose Logical: should be putted a prior on alpha?
#' @param tol Logical: should be putted a prior on alpha?
#' @param random Logical: should be putted a prior on alpha?
#' @param pstep Logical: should be putted a prior on alpha?
#' @export
#' 
DLSBP_VB <- function(y, X, H = 4, maxiter = 1000, prior = NULL, verbose = TRUE, tol=1e-4, useEM=FALSE, pstep=10) {
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  X <- Matrix(X)
  log2pi <- log(2*pi)
  
  # In prior hyperparameters not defined, use the MLE, defined as follow
  if (is.null(prior)) {
    prior <- list(b = rep(0, p), 
                  B = Diagonal(p, Inf), 
                  mu_mu =  0, 
                  tau_mu = 0, 
                  a_tau = 1, 
                  b_tau = 0)
  }
  
  # Hyperparameters
  b <- prior$b
  mu_mu  <- prior$mu_mu
  B <- prior$B; P <- solve(B); Pb <- P%*%b
  tau_mu <- prior$tau_mu
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Initialization
  tau_tilde <- numeric(H)
  mu_tilde  <- numeric(H)
  tau       <- as.numeric(rep(1/diff(quantile(y,c(.25,0.75))),H))
  ltau      <- a_tilde <- b_tilde <- numeric(H) 
  xi        <- matrix(1,n,H-1)
  
  mu_beta    <- matrix(rnorm((H-1)*p,0,0.01), H - 1, p) # A little jitter is added
  Sigma_beta <- list(); for(h in 1:(H-1)) Sigma_beta[[h]] <- Matrix(0,p,p)
  
  rho        <- matrix(runif(n*(H-1)),n,H-1)
  z          <- t(apply(rho,1,sb))
  
  if(useEM) {
    EM       <- DLSBP_EM(y, X, H=H, prior=prior, maxiter=500,tol=1e-4)
    tau      <- as.numeric(EM$tau)
    ltau     <- digamma(tau)
    mu_tilde <- EM$mu
    z        <- EM$z
    mu_beta  <- as.matrix(EM$beta)
    for(h in 1:(H-1)) {
      rho[,h] <- plogis(as.numeric(X%*%mu_beta[h,]))
    }
  }
  
  lower1 <- lower4 <- lower6 <- numeric(H)
  lower2 <- lower3 <- lower5 <- lower7 <- numeric(H-1)
  
  lowerbound <- -Inf
  
  # VB Algorithm
  for (r in 1:maxiter) {
  
    # Step 2,3,4: performed within the cluster.
    for (h in 1:H) {
      # Step 2 - Logistic regressions
      if (h < H) {
        
        # Variability
        omega <- tanh(xi[,h]/2)/(2 * xi[,h])
        omega[is.nan(omega)] <- 1/4
        #diag(Omega) <- omega
        
        # Parameters
        Sigma_beta[[h]] <- solve(crossprod(X*sqrt(omega)) + P)
        mu_beta[h, ]    <- as.numeric(Sigma_beta[[h]] %*% (crossprod(X,rho[, h] - 1/2) + Pb))
        
        # Step 4 - Maximization of variational parameters
        xi[,h] <- sqrt(as.numeric((X %*% mu_beta[h, ])^2) + mahalanobis(X,center=0,cov=as.matrix(Sigma_beta[[h]]),inverted=TRUE))
      }
      
      # Step 3 - Parameter updating
      tau_tilde[h] <- tau[h]*sum(z[,h]) + tau_mu
      mu_tilde[h]  <- (tau[h]*sum(z[,h]* y) + tau_mu*mu_mu)/tau_tilde[h]
      
      residual  <- y^2 - 2*y*mu_tilde[h] + mu_tilde[h]^2 + 1/tau_tilde[h]
      a_tilde[h]<- a_tau  + sum(z[, h])/2
      b_tilde[h]   <- b_tau + sum(z[, h] * residual)/2
      tau[h]    <- a_tilde[h] / b_tilde[h]
      ltau[h]   <- digamma(a_tilde[h]) - log(b_tilde[h])
    }
    
    # Step 1 - Cluster allocation
    Variational <- Variational_step(rho, y, as.matrix(X), mu_beta, mu_tilde, tau_tilde, tau, ltau)
    rho <- pmin(pmax(Variational$rho,1e-16),1-1e-16)
    z   <- Variational$z
    
    
    # ---------Lowerbound computatations
    for (h in 1:H) {
      residual   <-  y^2 - 2*y*mu_tilde[h] + mu_tilde[h]^2 + 1/tau_tilde[h]
      lower1[h]  <-  sum(z[,h]*(ltau[h]/2 - log2pi/2 - tau[h]/2*residual))
      lower4[h]  <- -log2pi/2 + log(tau_mu)/2 - tau_mu/2*(mu_tilde[h]^2 + 1/tau_tilde[h] + mu_mu^2 -2*mu_mu*mu_tilde[h]) + a_tau*log(b_tau) - lgamma(a_tau) - b_tau*tau[h] + (a_tau-1)*ltau[h]
      lower6[h]  <-  -log2pi/2 + log(tau_tilde[h])/2 - 0.5 + a_tilde[h]*log(b_tilde[h]) - lgamma(a_tilde[h]) - b_tilde[h]*tau[h] + (a_tilde[h]-1)*ltau[h]
      
      if (h < H) {
        eta         <- X%*%mu_beta[h,]
        lower2[h]   <- sum(rho[,h]*eta - (eta + xi[,h])/2 + log(plogis(xi[,h])))
        lower3[h]   <- -p*log2pi/2 - determinant(B,logarithm=TRUE)$modulus/2 - 0.5*(sum(diag(P%*%Sigma_beta[[h]])) + t(mu_beta[h,] - b)%*%P%*%(mu_beta[h,] - b))
        lower5[h]  <-  -p*log2pi/2 - determinant(Sigma_beta[[h]],logarithm=TRUE)$modulus/2 - p/2
        lower7[h]  <-  sum(rho[,h]*log(rho[,h]) + (1 - rho[,h])*log(1-rho[,h]))
      }
    }
    # Convergence checks
    lowerbound_new <- sum(lower1) + sum(lower2) + sum(lower3) + sum(lower4) - sum(lower5) - sum(lower6) - sum(lower7)
    
    # Break the loop at convergence
    if(lowerbound_new - lowerbound <tol)break
    
    # Otherwise continue!
    lowerbound <- lowerbound_new
    
    # Display status
    if (verbose) {
      if(r%%pstep==0) cat(paste("Lower-bound:",round(lowerbound,15),", iteration:", r, "\n",sep=""))
    }
  }
  
  # Output
  cluster <- apply(z,1,which.max)
  pred    <- prediction(y,as.matrix(X),mu_beta,mu_tilde,pmax(0,(a_tilde-1)/b_tilde))
  list(mu_beta = mu_beta, Sigma_beta=Sigma_beta, mu_tilde = mu_tilde, tau_tilde=tau_tilde, a_tilde=a_tilde, b_tilde=b_tilde, cluster=cluster,z=z, lowerbound=lowerbound,pred=pred)
}