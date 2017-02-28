DLSBP_ECM_multi <- function(y, X1, X2, H, prior, maxiter, tol, method_init, verbose, verbose_step) {
  
  # Fixed quantities
  n  <- length(y)
  p_kernel <- NCOL(X1); p_mixing <- NCOL(X2)
  
  # Hyperparameters
  b_mixing <- prior$b_mixing
  b_kernel <- prior$b_kernel
  B_mixing <- prior$B_mixing; P_mixing <- solve(B_mixing); Pb_mixing <- P_mixing %*% b_mixing;eig_Bmixing <- eigen(B_mixing)
  B_kernel<- prior$B_kernel; P_kernel <- solve(B_kernel); Pb_kernel <- P_kernel %*% b_kernel;eig_Bmixingkernel <- eigen(B_kernel)
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Initializing quantities
  tau   <- as.numeric(rep(1/diff(quantile(y,c(.25,0.75))), H))
  beta_kernel <- cbind(rep(median(y),H),matrix(0, H, p_kernel-1))
  beta_mixing  <- matrix(.1, H - 1, p_mixing) 
  logpost <- -Inf
  
  if(method_init=="cluster") {
    if(verbose) cat("Clustering observation before starting ECM algorithm...\n")
    cluster <- clara(X2,H)$clustering
    z       <- model.matrix(y~as.factor(cluster)-1)
  } else if (method_init=="random"){
    if(verbose) cat("Random initialization of the ECM algorithm...\n")
    beta_kernel <- matrix(rnorm(p_kernel*H,0,.5), H,p_kernel)
    beta_mixing  <- matrix(rnorm(p_mixing*(H-1),0,.1), H-1,p_mixing)
    Expectation <- Expectation_step_multi(y, X1, X2, beta_mixing, beta_kernel, tau)
    z <- Expectation$z
  } else {
    if(verbose) cat("Deterministic initialization of the ECM algorithm...\n")
    Expectation <- Expectation_step_multi(y, X1, X2, beta_mixing, beta_kernel, tau)
    z <- Expectation$z
  }
  
  # ECM Algorithm
  if(verbose) cat("Starting the ECM optimization. \n")
  for (r in 1:maxiter) {
    for (h in 1:H) {
      if (h < H) {
        # (Expectation) Step  - PolyaGamma
        linpred <- as.numeric(X2%*%beta_mixing[h, ])
        
        # Fixed quantities
        z_sum <- rowSums(z[, h:H])
        omega <- z_sum/(2 * linpred) * tanh(linpred/2)
        is_omega_nan <- is.nan(omega); omega[is_omega_nan] <- z_sum[is_omega_nan]/4
        
        # (Maximization) Step
        Sigma_mixing1 <- crossprod(X2*sqrt(omega)) + P_mixing 
        beta_mixing[h, ] <- solve(Sigma_mixing1, crossprod(X2 ,z[, h] - z_sum/2) + Pb_mixing)
      }
      
      # (Maximization) Step
      Sigma_kernel1 <- tau[h]*crossprod(X1*sqrt(z[, h])) + P_kernel
      beta_kernel[h, ] <- solve(Sigma_kernel1, crossprod(X1*z[, h]*tau[h] ,y) + Pb_kernel)
      residual  <- as.numeric(y - X1 %*% beta_kernel[h, ])
      tau[h]    <- max(0,(a_tau - 1 + sum(z[, h])/2)/(b_tau + sum(z[, h] * residual^2)/2))
    }
    
    # (Expectation) Step - Cluster allocation
    Expectation <- Expectation_step_multi(y, X1, X2, beta_mixing, beta_kernel, tau)
    z <- Expectation$z
    
    # Computation of the log-posterior
    loglik         <- Expectation$loglik
    logpriors      <- sum(ldmvnorm(beta_mixing, b_mixing, eig_Bmixing)) + sum(ldmvnorm(beta_kernel, b_kernel,eig_Bmixingkernel)) + sum(dgamma(tau, a_tau,b_tau,log=TRUE))
    logpost_new    <- loglik + logpriors
    
    if(!is.finite(logpriors)) logpost_new <- loglik
    
    # Break the loop at convergence
    if((logpost_new - logpost)<tol) {
      if(verbose) cat(paste("Convergence reached after",r,"iterations."))
      break
    }
    
    logpost <- logpost_new
    
    # Display status
    if (verbose) {
      if(r%%verbose_step==0) cat(paste("log-posterior:",round(logpost,4),", iteration:", r, "\n",sep=""))
    }
  }
  if(r==maxiter) warning(paste("Convergence has not been reached after",r,"iterations."))
  
  # Output
  cluster <- apply(z,1,which.max)
  out <- list(param=list(beta_mixing = beta_mixing, beta_kernel = beta_kernel, tau=tau), cluster=cluster, z=z, logposterior=logpost)
  return(out)
}

DLSBP_ECM_univ <- function(y, X, H, prior, maxiter, tol, method_init, verbose, verbose_step) {
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)

  # Hyperparameters
  b_mixing <- prior$b_mixing
  B_mixing <- prior$B_mixing; P_mixing <- solve(B_mixing); Pb_mixing <- P_mixing %*% b_mixing;eig_Bmixing <- eigen(B_mixing)
  mu_mu  <- as.numeric(prior$b_kernel)
  tau_mu <- as.numeric(1/prior$B_kernel)
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Initialization
  tau     <- as.numeric(rep(1/diff(quantile(y,c(.25,0.75))), H))
  mu      <- rep(median(y),H)
  beta    <- matrix(1e-2, H - 1, p)    # A little jitter is added
  logpost <- -Inf
  
  if(method_init=="cluster") {
    if(verbose) cat("Clustering observation before starting ECM algorithm...\n")
    cluster <- clara(X,H)$clustering 
    z       <- model.matrix(y~as.factor(cluster)-1)
  } else if (method_init=="random"){
    if(verbose) cat("Random initialization of the ECM algorithm...\n")
    beta  <- matrix(rnorm(p*(H-1),0,.1), H-1,p)
    mu    <- rnorm(H,mean(y),.1)
    Expectation <- Expectation_step(y, X, beta, mu, tau)
    z <- Expectation$z
  } else {
    if(verbose) cat("Deterministic initialization of the ECM algorithm...\n")
    Expectation <- Expectation_step(y, X, beta, mu, tau)
    z <- Expectation$z
  }
  
  # ECM Algorithm
  if(verbose) cat("Starting the ECM optimization. \n")
  for (r in 1:maxiter) {
    for (h in 1:H) {
      if (h < H) {
        # (Expectation) Step- PolyaGamma
        linpred <- as.numeric(X%*%beta[h, ])
        
        # Fixed quantities
        z_sum <- rowSums(z[, h:H])
        omega <- z_sum/(2 * linpred) * tanh(linpred/2)
        is_omega_nan <- is.nan(omega); omega[is_omega_nan] <- z_sum[is_omega_nan]/4
        
        # (Maximization) Step 
        Sigma_beta1 <- crossprod(X*sqrt(omega)) + P_mixing 
        beta[h, ] <- solve(Sigma_beta1, crossprod(X ,z[, h] - z_sum/2) + Pb_mixing)
      }
      
      # (Maximization) Step 
      tau_tilde <- tau[h]*sum(z[,h]) + tau_mu
      mu[h]     <- (tau[h]*sum(z[,h]* y) + tau_mu*mu_mu)/tau_tilde
      
      residual <- as.numeric(y - mu[h])
      tau[h] <- max(0,(a_tau - 1 + sum(z[, h])/2)/(b_tau + sum(z[, h] * residual^2)/2))
    }
    
    # (Expectation) Step - Cluster allocation
    Expectation <- Expectation_step(y, X, beta, mu, tau)
    z <- Expectation$z
    
    # Convergence checks
    loglik         <- Expectation$loglik
    logpriors      <- sum(ldmvnorm(beta, b_mixing, eig_Bmixing)) + sum(dnorm(mu,  mu_mu,1/sqrt(tau_mu),log=TRUE)) + sum(dgamma(tau, a_tau,b_tau,log=TRUE))
    logpost_new    <- loglik + logpriors
    if(!is.finite(logpriors)) logpost_new <- loglik
    
    # Break the loop at convergence
    if((logpost_new - logpost)<tol) {
      if(verbose) cat(paste("Convergence reached after",r,"iterations."))
      break
    }
    
    # Otherwise continue!
    logpost <- logpost_new
    
    # Display status
    if (verbose) {
      if(r%%verbose_step==0) cat(paste("log-posterior:",round(logpost,4),", iteration:", r, "\n",sep=""))
    }
  }
  if(r==maxiter) warning(paste("Convergence has not been reached after",r,"iterations.\n"))
  
  # Output
  cluster <- apply(z,1,which.max)
  out <- list(param=list(beta_mixing = beta, beta_kernel = mu, tau=tau), cluster=cluster,z=z, logposterior=logpost)
  return(out)
}

