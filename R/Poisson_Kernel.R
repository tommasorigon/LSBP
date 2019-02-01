Poisson_VB <- function(y, X, H, prior, maxiter, tol, method_init, verbose, verbose_step) {
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  log2pi <- log(2 * pi)  # Compute the logarithm of pi.
  
  # Hyperparameters
  b_mixing <- prior$b_mixing
  B_mixing <- prior$B_mixing
  P_mixing <- solve(B_mixing)
  Pb_mixing <- P_mixing %*% b_mixing
  eig_B <- eigen(B_mixing)
  logdetB <- sum(log(eig_B$values))
  
  # Prior hypeparameters for the mean parameter of the Poisson kernel
  # For notational consistency, the mean parameters are called "tau"
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # The posterior mean (and its logarithm) of the "tau" parameters are stored in "a_tilde","b_tilde"
  tau  <- rep(mean(y), H)
  ltau <- a_tilde <- b_tilde <- numeric(H)
  
  # All these quantities are, instead, related to the mixing.
  xi        <- matrix(1, n, H - 1)
  rho       <- matrix(0, n, H - 1)
  linpred   <- matrix(0, n, H - 1)
  mu_mixing <- matrix(0, H - 1, p)
  Sigma_mixing <- array(0, c(H - 1, p, p))
  
  
  if (method_init == "cluster") {
    if (verbose) 
      cat("Clustering observation before starting VB algorithm...\n")
    cluster <- clara(X, H)$clustering
    z <- model.matrix(y ~ as.factor(cluster) - 1)
    z <- z + 1e-08
    z <- z/rowSums(z)
    for (h in 1:(H - 1)) {
      rho[, h] <- z[, h]/rowSums(z[, h:H])
    }
  } else {
    if (verbose) 
      cat("Random initialization of the VB algorithm...\n")
    # Random initialization
    mu_mixing <- matrix(rnorm((H - 1) * p, -0.1, 0.5), H - 1, p)
    Sigma_mixing <- array(0, c(H - 1, p, p))
    
    # Random and uniform starting probabilities
    rho <- matrix(runif(n * (H - 1)), n, H - 1)
    z <- t(apply(rho, 1, sb))
  }
  
  
  # Initialization of different pieces of the lowerbound
  lower1 <- lower4 <- lower6 <- numeric(H)
  lower2 <- lower3 <- lower5 <- lower7 <- numeric(H - 1)
  lowerbound <- -Inf
  
  # VB Algorithm
  for (r in 1:maxiter) {
    for (h in 1:H) {
      
      # Stick-breaking component
      if (h < H) {
        omega <- tanh(xi[, h]/2)/(2 * xi[, h])
        omega[is.nan(omega)] <- 1/4
        
        # Parameters for beta_mixing
        Sigma_mixing[h, , ] <- solve(crossprod(X * sqrt(omega)) + P_mixing)
        mu_mixing[h, ]      <- as.numeric(Sigma_mixing[h, , ] %*% (crossprod(X, rho[, h] - 1/2) + Pb_mixing))
        linpred[, h]        <- X %*% mu_mixing[h, ]
        
        # Variational step
        xi[, h] <- sqrt(linpred[, h]^2 + rowSums(X %*% Sigma_mixing[h, , ] * X))
      }
      
      # Kernel parameters
      a_tilde[h] <- a_tau + sum(z[, h] * y)
      b_tilde[h] <- b_tau + sum(z[, h])
      tau[h] <- a_tilde[h]/b_tilde[h]
      ltau[h] <- digamma(a_tilde[h]) - log(b_tilde[h])
    }
    
    # Mixture probabilities
    Variational <- Variational_step_pois(rho, linpred, y, tau, ltau)
    rho <- pmin(pmax(Variational$rho, 1e-16), 1 - 1e-16)  # Added for numerical stability
    z <- Variational$z
    
    # Lowerbound computations
    for (h in 1:H) {
      lower1[h] <-  sum(z[, h] * (y*ltau[h] - tau[h]))
      lower4[h] <-  a_tau * log(b_tau) - lgamma(a_tau) - b_tau * tau[h] + (a_tau - 1) * ltau[h]
      lower6[h] <-  a_tilde[h] * log(b_tilde[h]) - lgamma(a_tilde[h]) - b_tilde[h] * tau[h] + (a_tilde[h] - 1) * ltau[h]
      
      if (h < H) {
        lower2[h] <- sum(rho[, h] * linpred[, h] - (linpred[, h] + xi[, h])/2 + log(plogis(xi[, h])))
        lower3[h] <- -p * log2pi/2 - logdetB/2 - 0.5 * (sum(diag(P_mixing %*% Sigma_mixing[h,, ])) + t(mu_mixing[h, ] - b_mixing) %*% P_mixing %*% (mu_mixing[h, ] - b_mixing))
        lower5[h] <- -p * log2pi/2 - ldet(Sigma_mixing[h, , ])/2 -  p/2
        lower7[h] <- sum(rho[, h] * log(rho[, h]) + (1 - rho[, h]) * log(1 - rho[, h]))
      }
    }
    # Convergence checks
    lowerbound_new <- sum(lower1) + sum(lower2) + sum(lower3) + sum(lower4) - sum(lower5) - sum(lower6) - 
      sum(lower7)
    # print(lowerbound_new) 
    
    # Break the loop at convergence
    if (lowerbound_new - lowerbound < tol) {
      if (verbose) 
        cat(paste("Convergence reached after", r, "iterations."))
      break
    }
    
    # Otherwise continue
    lowerbound <- lowerbound_new
    
    # Display status
    if (verbose) {
      if (r%%verbose_step == 0) 
        cat(paste("Lower-bound: ", round(lowerbound, 4), ", iteration: ", r, "\n", sep = ""))
    }
  }
  
  if (r == maxiter) 
    warning(paste("Convergence has not been reached after", r, "iterations."))
  
  # Output
  cluster <- apply(z, 1, which.max)
  list(param = list(mu_mixing = mu_mixing, Sigma_mixing = Sigma_mixing, 
                    a_tilde = a_tilde, b_tilde = b_tilde), cluster = cluster, z = z, lowerbound = lowerbound)
}

Poisson_Gibbs <- function(y, X, prior, H, R, burn_in, method_init, verbose, verbose_step) {
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  
  # Hyperparameters
  b_mixing <- prior$b_mixing
  B_mixing <- prior$B_mixing
  P_mixing <- solve(B_mixing)
  Pb_mixing <- P_mixing %*% b_mixing
  eig_B <- eigen(B_mixing)
  
  # Prior hyperparameters
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Initialization
  tau <- rep(mean(y),H)
  beta_mixing <- matrix(0.1, H - 1, p)  # A little jitter is added
  
  # Output
  beta_mixing_out <- array(0, c(R, H - 1, p))
  tau_out <- matrix(0, R, H)
  logpost <- numeric(R + burn_in)
  
  if (method_init == "cluster") {
    if (verbose) 
      cat("Clustering observation before starting Gibbs sampling...\n")
    G <- clara(X, H)$clustering
  } else {
    G <- G_pois_update(y, X, beta_mixing, tau)$G
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
        linh  <- as.numeric(Xh %*% beta_mixing[h, ])
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
      tau[h] <- rgamma(1, a_tau + sum(yG), b_tau + nG)
    }
    
    # Cluster allocation (Step 1)
    Update <- G_pois_update(y, X, beta_mixing,  tau)
    G      <- c(Update$G)

    # Convergence checks
    loglik <- Update$loglik
    logpriors <- sum(ldmvnorm(beta_mixing, b_mixing, eig_B)) + sum(dgamma(tau, a_tau, b_tau, log = TRUE))
    logpost[r] <- loglik + logpriors
    
    # Output
    if (r > burn_in) {
      beta_mixing_out[r - burn_in, , ] <- beta_mixing
      tau_out[r - burn_in, ] <- tau
    }
    
    if (verbose) {
      if (r%%verbose_step == 0) 
        cat(paste("Sampling iteration: ", r, ", logposterior: ", round(logpost[r], 4), ".\n", 
                  sep = ""))
    }
  }
  list(param = list(beta_mixing = beta_mixing_out, tau = tau_out), logposterior = logpost)
}
