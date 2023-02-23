LSBP_VB_univ <- function(y, X, H, prior, maxiter, tol, method_init, verbose, verbose_step) {
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  log2pi <- log(2 * pi) # Compute the logarithm of pi.

  # Hyperparameters
  b_mixing <- prior$b_mixing
  B_mixing <- prior$B_mixing
  P_mixing <- solve(B_mixing)
  Pb_mixing <- P_mixing %*% b_mixing
  eig_B <- eigen(B_mixing)
  logdetB <- sum(log(eig_B$values))
  mu_mu <- as.numeric(prior$b_kernel)
  tau_mu <- as.numeric(1 / prior$B_kernel)
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau

  # Initialization
  tau_tilde <- numeric(H)
  mu_tilde <- rep(median(y), H)
  tau <- as.numeric(rep(1 / diff(quantile(y, c(0.25, 0.75))), H))
  ltau <- a_tilde <- b_tilde <- numeric(H)
  xi <- matrix(1, n, H - 1)
  residual <- matrix(0, n, H)
  linpred <- matrix(0, n, H - 1)
  rho <- matrix(0, n, H - 1)
  mu_mixing <- matrix(0, H - 1, p)
  Sigma_mixing <- array(0, c(H - 1, p, p))


  if (method_init == "cluster") {
    if (verbose) {
      cat("Clustering observation before starting VB algorithm...\n")
    }
    cluster <- clara(X, H)$clustering
    z <- model.matrix(y ~ as.factor(cluster) - 1)
    z <- z + 1e-08
    z <- z / rowSums(z)
    for (h in 1:(H - 1)) {
      rho[, h] <- z[, h] / rowSums(z[, h:H])
    }
  } else {
    if (verbose) {
      cat("Random initialization of the VB algorithm...\n")
    }
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
      # Logistic regressions
      if (h < H) {
        omega <- tanh(xi[, h] / 2) / (2 * xi[, h])
        omega[is.nan(omega)] <- 1 / 4

        # Parameters for beta_mixing
        Sigma_mixing[h, , ] <- solve(crossprod(X * sqrt(omega)) + P_mixing)
        mu_mixing[h, ] <- as.numeric(Sigma_mixing[h, , ] %*% (crossprod(X, rho[, h] - 1 / 2) +
          Pb_mixing))
        linpred[, h] <- X %*% mu_mixing[h, ]

        # Variational step
        xi[, h] <- sqrt(linpred[, h]^2 + rowSums(X %*% Sigma_mixing[h, , ] * X))
      }

      # Parameters for posterior of mu and tau
      tau_tilde[h] <- tau[h] * sum(z[, h]) + tau_mu
      mu_tilde[h] <- (tau[h] * sum(z[, h] * y) + tau_mu * mu_mu) / tau_tilde[h]

      residual[, h] <- y^2 - 2 * y * mu_tilde[h] + mu_tilde[h]^2 + 1 / tau_tilde[h]
      a_tilde[h] <- a_tau + sum(z[, h]) / 2
      b_tilde[h] <- b_tau + sum(z[, h] * residual[, h]) / 2
      tau[h] <- a_tilde[h] / b_tilde[h]
      ltau[h] <- digamma(a_tilde[h]) - log(b_tilde[h])
    }

    # Mixture probabilities
    Variational <- Variational_step(rho, linpred, residual, tau, ltau)
    rho <- pmin(pmax(Variational$rho, 1e-16), 1 - 1e-16) # Added for numerical stability
    z <- Variational$z

    # Lowerbound computatations
    for (h in 1:H) {
      lower1[h] <- sum(z[, h] * (ltau[h] / 2 - log2pi / 2 - tau[h] / 2 * residual[, h]))
      lower4[h] <- -log2pi / 2 + log(tau_mu) / 2 - tau_mu / 2 * (mu_tilde[h]^2 + 1 / tau_tilde[h] + mu_mu^2 -
        2 * mu_mu * mu_tilde[h]) + a_tau * log(b_tau) - lgamma(a_tau) - b_tau * tau[h] + (a_tau -
        1) * ltau[h]
      lower6[h] <- -log2pi / 2 + log(tau_tilde[h]) / 2 - 0.5 + a_tilde[h] * log(b_tilde[h]) - lgamma(a_tilde[h]) -
        b_tilde[h] * tau[h] + (a_tilde[h] - 1) * ltau[h]

      if (h < H) {
        lower2[h] <- sum(rho[, h] * linpred[, h] - (linpred[, h] + xi[, h]) / 2 + log(plogis(xi[
          ,
          h
        ])))
        lower3[h] <- -p * log2pi / 2 - logdetB / 2 - 0.5 * (sum(diag(P_mixing %*% Sigma_mixing[h, , ])) + t(mu_mixing[h, ] - b_mixing) %*% P_mixing %*% (mu_mixing[h, ] - b_mixing))
        lower5[h] <- -p * log2pi / 2 - ldet(Sigma_mixing[h, , ]) / 2 -
          p / 2
        lower7[h] <- sum(rho[, h] * log(rho[, h]) + (1 - rho[, h]) * log(1 - rho[, h]))
      }
    }
    # Convergence checks
    lowerbound_new <- sum(lower1) + sum(lower2) + sum(lower3) + sum(lower4) - sum(lower5) - sum(lower6) -
      sum(lower7)

    # Break the loop at convergence
    if (lowerbound_new - lowerbound < tol) {
      if (verbose) {
        cat(paste("Convergence reached after", r, "iterations."))
      }
      break
    }

    # Otherwise continue
    lowerbound <- lowerbound_new

    # Display status
    if (verbose) {
      if (r %% verbose_step == 0) {
        cat(paste("Lower-bound: ", round(lowerbound, 4), ", iteration: ", r, "\n", sep = ""))
      }
    }
  }

  if (r == maxiter) {
    warning(paste("Convergence has not been reached after", r, "iterations."))
  }

  # Output
  cluster <- apply(z, 1, which.max)
  list(param = list(
    mu_mixing = mu_mixing, Sigma_mixing = Sigma_mixing, mu_kernel = mu_tilde, Sigma_kernel = 1 / tau_tilde,
    a_tilde = a_tilde, b_tilde = b_tilde
  ), cluster = cluster, z = z, lowerbound = lowerbound)
}

LSBP_VB_multi <- function(y, X1, X2, H, prior, maxiter, tol, method_init, verbose, verbose_step) {
  # Fixed quantities
  n <- length(y)
  p_kernel <- NCOL(X1)
  p_mixing <- NCOL(X2)
  log2pi <- log(2 * pi) # Compute the logarithm of pi.

  # Hyperparameters
  b_mixing <- prior$b_mixing
  b_kernel <- prior$b_kernel
  B_mixing <- prior$B_mixing
  P_mixing <- solve(B_mixing)
  Pb_mixing <- P_mixing %*% b_mixing
  eig_B <- eigen(B_mixing)
  logdetB <- sum(log(eig_B$values))
  B_kernel <- prior$B_kernel
  P_kernel <- solve(B_kernel)
  Pb_kernel <- P_kernel %*% b_kernel
  eig_Bkernel <- eigen(B_kernel)
  logdetB_kernel <- sum(log(eig_Bkernel$values))
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau

  # Initialization
  tau <- as.numeric(rep(1 / diff(quantile(y, c(0.25, 0.75))), H))
  ltau <- a_tilde <- b_tilde <- numeric(H)
  xi <- matrix(1, n, H - 1)
  residual <- matrix(0, n, H)
  linpred_mixing <- matrix(0, n, H - 1)
  linpred_kernel <- matrix(0, n, H)
  rho <- matrix(0, n, H - 1)
  mu_mixing <- matrix(0, H - 1, p_mixing)
  Sigma_mixing <- array(0, c(H - 1, p_mixing, p_mixing))
  mu_kernel <- matrix(0, H, p_kernel)
  Sigma_kernel <- array(0, c(H, p_kernel, p_kernel))

  if (method_init == "cluster") {
    if (verbose) {
      cat("Clustering observation before starting VB algorithm...\n")
    }
    cluster <- clara(X2, H)$clustering
    z <- model.matrix(y ~ as.factor(cluster) - 1)
    z <- z + 1e-08
    z <- z / rowSums(z)
    for (h in 1:(H - 1)) {
      rho[, h] <- z[, h] / rowSums(z[, h:H])
    }
  } else {
    if (verbose) {
      cat("Random initialization of the VB algorithm...\n")
    }

    # Random initialization
    mu_mixing <- matrix(rnorm((H - 1) * p_mixing, -0.1, 0.1), H - 1, p_mixing)
    Sigma_mixing <- array(0, c(H - 1, p_mixing, p_mixing))

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
      # Logistic regressions
      if (h < H) {
        omega <- tanh(xi[, h] / 2) / (2 * xi[, h])
        omega[is.nan(omega)] <- 1 / 4

        # Parameters for beta_mixing
        Sigma_mixing[h, , ] <- solve(crossprod(X2 * sqrt(omega)) + P_mixing)
        mu_mixing[h, ] <- Sigma_mixing[h, , ] %*% (crossprod(X2, rho[, h] - 1 / 2) + Pb_mixing)
        linpred_mixing[, h] <- X2 %*% mu_mixing[h, ]

        # Variational step
        xi[, h] <- sqrt(linpred_mixing[, h]^2 + rowSums(X2 %*% Sigma_mixing[h, , ] * X2))
      }

      # Parameters for posterior of mu and tau

      Sigma_kernel[h, , ] <- solve(tau[h] * crossprod(X1 * sqrt(z[, h])) + P_kernel)
      mu_kernel[h, ] <- Sigma_kernel[h, , ] %*% (crossprod(X1 * z[, h] * tau[h], y) + Pb_kernel)

      # tau_tilde[h] <- tau[h]*sum(z[,h]) + tau_mu mu_tilde[h] <- (tau[h]*sum(z[,h]* y) +
      # tau_mu*mu_mu)/tau_tilde[h]

      linpred_kernel[, h] <- X1 %*% mu_kernel[h, ]
      residual[, h] <- y^2 - 2 * y * linpred_kernel[, h] + linpred_kernel[, h]^2 + rowSums(X1 %*%
        Sigma_kernel[h, , ] * X1) # E( (y - xkernel)^2) = y^2 - 2*y*E(eta) + E(eta^2)
      a_tilde[h] <- a_tau + sum(z[, h]) / 2
      b_tilde[h] <- b_tau + sum(z[, h] * residual[, h]) / 2
      tau[h] <- a_tilde[h] / b_tilde[h]
      ltau[h] <- digamma(a_tilde[h]) - log(b_tilde[h])
    }

    # Mixture probabilities
    Variational <- Variational_step(rho, linpred_mixing, residual, tau, ltau)
    rho <- pmin(pmax(Variational$rho, 1e-16), 1 - 1e-16) # Added for numerical stability
    z <- Variational$z

    # Lowerbound computatations
    for (h in 1:H) {
      lower1[h] <- sum(z[, h] * (ltau[h] / 2 - log2pi / 2 - tau[h] / 2 * residual[, h]))
      lower4[h] <- -p_kernel * log2pi / 2 - logdetB_kernel / 2 - 0.5 * (sum(diag(P_kernel %*% Sigma_kernel[h, , ])) + t(mu_kernel[h, ] - b_kernel) %*% P_kernel %*% (mu_kernel[h, ] - b_kernel)) +
        a_tau * log(b_tau) - lgamma(a_tau) - b_tau * tau[h] + (a_tau - 1) * ltau[h]
      lower6[h] <- -p_kernel * log2pi / 2 - ldet(Sigma_kernel[h, , ]) / 2 -
        p_kernel / 2 + a_tilde[h] * log(b_tilde[h]) - lgamma(a_tilde[h]) - b_tilde[h] * tau[h] +
        (a_tilde[h] - 1) * ltau[h]

      if (h < H) {
        lower2[h] <- sum(rho[, h] * linpred_mixing[, h] - (linpred_mixing[, h] + xi[, h]) / 2 +
          log(plogis(xi[, h])))
        lower3[h] <- -p_mixing * log2pi / 2 - logdetB / 2 - 0.5 * (sum(diag(P_mixing %*% Sigma_mixing[h, , ])) + t(mu_mixing[h, ] - b_mixing) %*% P_mixing %*% (mu_mixing[h, ] - b_mixing))
        lower5[h] <- -p_mixing * log2pi / 2 - ldet(Sigma_mixing[h, , ]) / 2 -
          p_mixing / 2
        lower7[h] <- sum(rho[, h] * log(rho[, h]) + (1 - rho[, h]) * log(1 - rho[, h]))
      }
    }
    # Convergence checks
    lowerbound_new <- sum(lower1) + sum(lower2) + sum(lower3) + sum(lower4) - sum(lower5) - sum(lower6) -
      sum(lower7)

    # Break the loop at convergence
    if (lowerbound_new - lowerbound < tol) {
      if (verbose) {
        cat(paste("Convergence reached after", r, "iterations."))
      }
      break
    }

    # Otherwise continue
    lowerbound <- lowerbound_new

    # Display status
    if (verbose) {
      if (r %% verbose_step == 0) {
        cat(paste("Lower-bound: ", round(lowerbound, 4), ", iteration: ", r, "\n", sep = ""))
      }
    }
  }

  if (r == maxiter) {
    warning(paste("Convergence has not been reached after", r, "iterations."))
  }

  # Output
  cluster <- apply(z, 1, which.max)
  list(param = list(
    mu_mixing = mu_mixing, Sigma_mixing = Sigma_mixing, mu_kernel = mu_kernel, Sigma_kernel = Sigma_kernel,
    a_tilde = a_tilde, b_tilde = b_tilde
  ), cluster = cluster, z = z, lowerbound = lowerbound)
}
