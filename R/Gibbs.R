LSBP_Gibbs_univ <- function(y, X, prior, H, R, burn_in, method_init, verbose, verbose_step) {
   
   # Fixed quantities
   n <- length(y)
   p <- NCOL(X)
   
   # Hyperparameters
   b_mixing <- prior$b_mixing
   B_mixing <- prior$B_mixing
   P_mixing <- solve(B_mixing)
   Pb_mixing <- P_mixing %*% b_mixing
   eig_B <- eigen(B_mixing)
   mu_mu <- as.numeric(prior$b_kernel)
   tau_mu <- as.numeric(1/prior$B_kernel)
   a_tau <- prior$a_tau
   b_tau <- prior$b_tau
   
   # Initialization
   tau <- as.numeric(rep(1/diff(quantile(y, c(0.25, 0.75))), H))
   mu <- rep(median(y), H)
   beta_mixing <- matrix(0.1, H - 1, p)  # A little jitter is added
   
   # Output
   beta_mixing_out <- array(0, c(R, H - 1, p))
   mu_out <- matrix(0, R, H)
   tau_out <- matrix(0, R, H)
   logpost <- numeric(R + burn_in)
   
   if (method_init == "cluster") {
      if (verbose) 
         cat("Clustering observation before starting Gibbs sampling...\n")
      G <- clara(X, H)$clustering
   } else {
      G <- G_update(y, X, beta_mixing, mu, tau)$G
   }
   
   # Gibbs sampling
   if (verbose) 
      cat(paste("Starting the Gibbs sampling with R = ", R, " and burn_in = ", burn_in, ". \n", sep = ""))
   for (r in 1:(R + burn_in)) {
      
      # Step 2 - 4: performed within the cluster.
      for (h in 1:H) {
         if (h < H) {
            # Subsetting observations
            index <- G > h - 1
            nh <- sum(index)
            zh <- G[index] == h
            Xh <- matrix(X[index, ], nh, p)
            
            # Polya-gamma weights, beta_mixing coefficients
            linh <- as.numeric(Xh %*% beta_mixing[h, ])
            omega <- rpg.devroye(num = nh, n = 1, z = linh)
            
            eig <- eigen(crossprod(Xh * sqrt(omega)) + P_mixing, symmetric = TRUE)
            Sigma_mixing <- crossprod(t(eig$vectors)/sqrt(eig$values))
            mu_mixing <- Sigma_mixing %*% (crossprod(Xh, zh - 1/2) + Pb_mixing)
            
            A1 <- t(eig$vectors)/sqrt(eig$values)
            beta_mixing[h, ] <- mu_mixing + c(matrix(rnorm(1 * p), 1, p) %*% A1)
         }
         
         # Mixture components
         indexG <- G == h
         nG <- sum(indexG)
         yG <- y[indexG]
         
         tau_tilde <- tau_mu + tau[h] * nG
         mu_tilde <- (tau_mu * mu_mu + tau[h] * sum(yG))/tau_tilde
         mu[h] <- rnorm(1, mu_tilde, 1/sqrt(tau_tilde))
         
         residual <- as.numeric(yG - mu[h])
         tau[h] <- rgamma(1, a_tau + nG/2, b_tau + sum(residual^2)/2)
      }
      
      # Cluster allocation (Step 1)
      Update <- G_update(y, X, beta_mixing, mu, tau)
      G <- c(Update$G)
      G_mixt <- c(Update$G_mixt)
      
      # Convergence checks
      loglik <- Update$loglik
      logpriors <- sum(ldmvnorm(beta_mixing, b_mixing, eig_B)) + sum(dnorm(mu, mu_mu, 1/sqrt(tau_mu), 
         log = TRUE)) + sum(dgamma(tau, a_tau, b_tau, log = TRUE))
      logpost[r] <- loglik + logpriors
      
      # Output
      if (r > burn_in) {
         beta_mixing_out[r - burn_in, , ] <- beta_mixing
         mu_out[r - burn_in, ] <- mu
         tau_out[r - burn_in, ] <- tau
      }
      
      if (verbose) {
         if (r%%verbose_step == 0) 
            cat(paste("Sampling iteration: ", r, ", logposterior: ", round(logpost[r], 4), ".\n", 
              sep = ""))
      }
   }
   list(param = list(beta_mixing = beta_mixing_out, beta_kernel = mu_out, tau = tau_out), logposterior = logpost)
}

LSBP_Gibbs_multi <- function(y, X1, X2,  H, R, prior, burn_in, method_init, verbose, verbose_step) {
   
   # Fixed quantities
   n <- length(y)
   p_kernel <- NCOL(X1)
   p_mixing <- NCOL(X2)
   
   # Hyperparameters
   b_mixing <- prior$b_mixing
   b_kernel <- prior$b_kernel
   B_mixing <- prior$B_mixing
   P_mixing <- solve(B_mixing)
   Pb_mixing <- P_mixing %*% b_mixing
   eig_B <- eigen(B_mixing)
   B_kernel <- prior$B_kernel
   P_kernel <- solve(B_kernel)
   Pb_kernel <- P_kernel %*% b_kernel
   eig_Bkernel <- eigen(B_kernel)
   a_tau <- prior$a_tau
   b_tau <- prior$b_tau
   
   # Initializing quantities
   tau <- as.numeric(rep(1/diff(quantile(y, c(0.25, 0.75))), H))
   beta_kernel <- cbind(rep(median(y), H), matrix(0, H, p_kernel - 1))
   beta_mixing <- matrix(0.1, H - 1, p_mixing)
   
   # Output
   beta_mixing_out <- array(0, c(R, H - 1, p_mixing))
   beta_kernel_out <- array(0, c(R, H, p_kernel))
   tau_out <- matrix(0, R, H)
   logpost <- numeric(R + burn_in)
   
   if (method_init == "cluster") {
      if (verbose) 
         cat("Clustering observation before starting Gibbs sampling...\n")
      G <- clara(X2, H)$clustering
   } else {
      G <- G_update_multi(y, X1, X2, beta_mixing, beta_kernel, tau)$G
   }
   
   # Gibbs sampling
   if (verbose) 
      cat(paste("Starting the Gibbs sampling with R = ", R, " and burn_in = ", burn_in, ". \n", sep = ""))
   for (r in 1:(R + burn_in)) {
      
      # Step 2 - 3: performed within the cluster.
      for (h in 1:H) {
         if (h < H) {
            # Subsetting observations
            index <- G > h - 1
            nh <- sum(index)
            zh <- G[index] == h
            X2h <- matrix(X2[index, ], nh, p_mixing)
            
            # Polya-gamma weights, beta_mixing coefficients
            linh <- as.numeric(X2h %*% beta_mixing[h, ])
            omega <- rpg.devroye(num = nh, n = 1, z = linh)
            
            eig <- eigen(crossprod(X2h * sqrt(omega)) + P_mixing, symmetric = TRUE)
            Sigma_mixing <- crossprod(t(eig$vectors)/sqrt(eig$values))
            mu_mixing <- Sigma_mixing %*% (crossprod(X2h, zh - 1/2) + Pb_mixing)
            
            A1 <- t(eig$vectors)/sqrt(eig$values)
            beta_mixing[h, ] <- mu_mixing + c(matrix(rnorm(1 * p_mixing), 1, p_mixing) %*% A1)
         }
         
         # Mixture components
         indexG <- G == h
         nG <- sum(indexG)
         yG <- y[indexG]
         X1G <- matrix(X1[indexG, ], nG, p_kernel)
         
         eig <- eigen(tau[h] * crossprod(X1G) + P_kernel, symmetric = TRUE)
         Sigma_kernel <- crossprod(t(eig$vectors)/sqrt(eig$values))
         mu_kernel <- Sigma_kernel %*% (tau[h] * crossprod(X1G, yG) + Pb_kernel)
         
         A1 <- t(eig$vectors)/sqrt(eig$values)
         beta_kernel[h, ] <- mu_kernel + c(matrix(rnorm(1 * p_kernel), 1, p_kernel) %*% A1)
         
         residual <- as.numeric(yG - X1G %*% beta_kernel[h, ])
         tau[h] <- rgamma(1, a_tau + nG/2, b_tau + sum(residual^2)/2)
      }
      # Cluster allocation
      Update <- G_update_multi(y, X1, X2, beta_mixing, beta_kernel, tau)
      G <- c(Update$G)
      G_mixt <- c(Update$G_mixt)
      
      # Convergence checks
      loglik <- Update$loglik
      logpriors <- sum(ldmvnorm(beta_mixing, b_mixing, eig_B)) + sum(ldmvnorm(beta_kernel, b_kernel, 
         eig_Bkernel)) + sum(dgamma(tau, a_tau, b_tau, log = TRUE))
      logpost[r] <- loglik + logpriors
      
      # Output
      if (r > burn_in) {
         beta_mixing_out[r - burn_in, , ] <- beta_mixing
         beta_kernel_out[r - burn_in, , ] <- beta_kernel
         tau_out[r - burn_in, ] <- tau
      }
      
      if (verbose) {
         if (r%%verbose_step == 0) 
            cat(paste("Sampling iteration: ", r, ", logposterior: ", round(logpost[r], 4), ".\n", 
              sep = ""))
      }
   }
   list(param = list(beta_mixing = beta_mixing_out, beta_kernel = beta_kernel_out, tau = tau_out), logposterior = logpost)
}
