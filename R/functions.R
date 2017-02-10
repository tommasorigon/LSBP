#'@importFrom Rcpp evalCpp sourceCpp
#'@importFrom BayesLogit rpg.devroye
#'@importFrom mvtnorm dmvnorm rmvnorm
#'@importFrom fields rdist Exponential Matern RadialBasis
#'@importFrom splines spline.des
#'@import Matrix
#'@useDynLib DLSBP

sb <- function(nu) {
  nu    <- c(nu,1)
  H     <- length(nu)
  prob  <- nu * c(1,cumprod(1 - nu[-H]))
  return(prob)
}


#' EM algorithm for the DLSBP model
#'
#' The dependent logit stick-breaking process (DLSBP) model estimated trhough the E(C)M algorithm, which provides the posterior mode.
#' 
#' @param y a vector containing the response vector
#' @param X a n x p design matrix containing the covariates
#' @param H an integer indicating the number of mixture components
#' @param maxiter an integer indicating the maximum number of iterations for the EM algorithm
#' @param prior a list containing prior hyperparameters (See details for a detailed description). If not provided, the MLE is computed.
#' @param tol a real number controlling the convergence criterium
#' @param random a logical value indicating whether a random initialization should be used. Default is FALSE.
#' @param verbose Logical: Should the logposterior be displayed while the algorithm is running?
#' 
#' @details  The prior argument is a list which should contains the following elements
#' \itemize{
#' \item \verb{b}. A p dimensional vector containing the prior mean of the Gaussian beta coefficients
#' \item \verb{B}. A p x p matrix representing the prior covariance of the Gaussian beta coefficients.
#' \item \verb{mu_mu}. A real number representing the prior mean for the kernel mean.
#' \item \verb{tau_mu}. A positive number representing the precision for the kernel mean.
#' \item \verb{a_tau}, \verb{b_tau}. The hyperparameters of a Gamma prior distribution for the kernel precision.
#' }
#' 
#' 
#' @return The output is a list containing the following quantities
#' \itemize{
#' \item \verb{beta}. A (H - 1) x p matrix containing the beta coefficients.
#' \item \verb{mu}. A H dimensional vector containing the mu coefficients.
#' \item \verb{tau}. A H dimensional vector containing the tau coefficients.
#' \item \verb{pred}. A n dimensional vector containing the the predicted values.
#' \item \verb{cluster}. A n dimensional vector containing, for each observation, the mixture component having with the highest probability.
#' \item \verb{z}. A n x H matrix containing the probabilities of belonging to each of the mixture components.
#' \item \verb{logposterior}. The logposterior of the DLSBP model at convergence.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2017), Logit stick-breaking priors for Bayesian density regression, ArXiv.
#' @examples 
#' library(DLSBP)
#' data(geyser)
#' x <- geyser$waiting
#' y <- geyser$duration
#' X <- cbind(1,x*I(x<=68) + 68*I(x > 68),((x-68)*I(x>68))) # Predictors
#' p <- NCOL(X)
#' prior   <- list(b = rep(0,p),        
#'                 B = diag(100,p), 
#'                 mu_mu =0,            
#'                 tau_mu = 0.001,     
#'                 a_tau = 2,           
#'                 b_tau = 0.001)
#' fit_EM   <- DLSBP_EM(y=y, X=X, H=3, prior=prior)
#' 
#' @export
#' 

DLSBP_EM <- function(y, X, H = 2, maxiter = 1000, prior = NULL, tol=1e-4, random=FALSE, verbose = TRUE) {
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  X <- Matrix(X)
  verbose_step = ceiling(maxiter / 50)
  
  # If prior hyperparameters are not defined, then use the MLE
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
  B <- Matrix(prior$B); P <- solve(B); Pb <- P %*% b
  tau_mu <- prior$tau_mu
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Initialization
  tau     <- rep(1/diff(quantile(y,c(.25,0.75))), H)
  mu      <- rep(0,H)
  beta    <- Matrix(1e-2, H - 1, p)    # A little jitter is added
  logpost <- -Inf
  
  if(random){
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
        
        # (Maximization) Step 3
        Sigma_beta1 <- crossprod(X*sqrt(omega)) + P #Equivalent to t(X) %*% Omega %*% X + P, but faster
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
    if((logpost_new - logpost)<tol) {
      cat(paste("Convergence reached after",r,"iterations."))
      break
    }
    
    # Otherwise continue!
    logpost <- logpost_new
    
    # Display status
    if (verbose) {
      if(r%%verbose_step==0) cat(paste("log-posterior:",round(logpost,4),", iteration:", r, "\n",sep=""))
    }
  }
  if(r==maxiter) warning(paste("Convergence has not been reached after",r,"iterations."))
  
  # Output
  cluster <- apply(z,1,which.max)
  pred    <- prediction(as.matrix(X), as.matrix(beta), mu, tau)
  list(beta = as.matrix(beta), mu = mu, tau=tau, pred = pred, cluster=cluster,z=z, logposterior=logpost)
}

#' Gibbs sampling for the DLSBP model
#'
#' The dependent logit stick-breaking process (DLSBP) model estimated through the Gibbs sampling.
#' 
#' @param y a vector containing the response vector.
#' @param X a n x p design matrix containing the covariates
#' @param H an integer indicating the number of mixture components
#' @param R an integer indicating the number of MCMC replications, after the burn-in period
#' @param burn.in an integer indicating the number of MCMC discarded as burn in period
#' @param prior a list containing prior hyperparameters (See details for a detailed description). 
#' @param verbose Logical: Should the MCMC iteration be displayed while the algorithm is running?
#' @param store_parameters If FALSE, for faster inference only predictions are stored
#' @param EM_initialization If TRUE, the starting values of the MCMC chain are those obtained by the EM algorithm.
#' 
#' @details  The prior argument is a list which should contains the following elements
#' \itemize{
#' \item \verb{b}. A p dimensional vector containing the prior mean of the Gaussian beta coefficients
#' \item \verb{B}. A p x p matrix representing the prior covariance of the Gaussian beta coefficients.
#' \item \verb{mu_mu}. A real number representing the prior mean for the kernel mean.
#' \item \verb{tau_mu}. A positive number representing the precision for the kernel mean.
#' \item \verb{a_tau}, \verb{b_tau}. The hyperparameters of a Gamma prior distribution for the kernel precision.
#' }
#' 
#' 
#' @return The output is a list containing the following quantities
#' \itemize{
#' \item \verb{beta}. An  R x (H - 1) x p array containing the beta coefficients at each step of the chain
#' \item \verb{mu}. A R x H matrix containing the mu coefficients at each step of the chain.
#' \item \verb{tau}. A R x H matrix containing the tau coefficients at each step of the chain.
#' \item \verb{pred}. A R x n matrix containing a sample from predictive distribution fore each observation at each step of the chain.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2017), Logit stick-breaking priors for Bayesian density regression, ArXiv.
#' @examples 
#' library(DLSBP)
#' data(geyser)
#' x <- geyser$waiting
#' y <- geyser$duration
#' X <- cbind(1,x*I(x<=68) + 68*I(x > 68),((x-68)*I(x>68))) # Predictors
#' p <- NCOL(X)
#' R <- 1000; burn.in <- 1000
#' prior   <- list(b = rep(0,p),        
#'                 B = diag(100,p), 
#'                 mu_mu =0,            
#'                 tau_mu = 0.001,     
#'                 a_tau = 2,           
#'                 b_tau = 0.001)
#' fit_Gibbs   <- DLSBP_Gibbs(y=y, X=X, H=3, prior=prior, R=R,burn.in=burn.in)
#' 
#' @export
#' 
DLSBP_Gibbs <- function(y, X, H = 2, R = 5000, burn.in = 2000, prior = NULL, verbose = TRUE, store_param=TRUE, EM_initialization=FALSE){
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  X <- Matrix(X)
  verbose_step = ceiling(R / 200)
  
  # Hyperparameters
  b <- prior$b
  mu_mu  <- prior$mu_mu
  B <- Matrix(prior$B); P <- solve(B); Pb <- P %*% b
  tau_mu <- prior$tau_mu
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Output 
  if(store_param) {
    beta_out   <- array(0, c(R, H - 1, p)) 
    mu_out     <- matrix(0, R, H)
    tau_out    <- matrix(0, R, H)
  }
  prediction <- matrix(0, R, n)
  
  # Initialization
  tau  <- rep(1/diff(quantile(y,c(.25,0.75))), H)
  mu   <- rep(0, H)
  beta <- Matrix(0, H - 1, p)
  
  if(EM_initialization) {
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
        
        Sigma_beta <- solve(crossprod(Xh*sqrt(omega)) + P) #Faster than solve(t(Xh) %*% Omega %*% Xh + P)
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
      if(store_param){ 
        beta_out[r - burn.in, , ] <- as.matrix(beta)
        mu_out[r - burn.in, ]     <- mu
        tau_out[r - burn.in, ]    <- tau
      }
      # Prediction
      prediction[r - burn.in, ] <- rnorm(n, mu[G_mixt], as.numeric(1/sqrt(tau[G_mixt])))
    }
    
    if (verbose) {
      if(r%%verbose_step==0) cat(paste("Sampling iteration: ", r, "\n",sep=""))
    }
  }
  if(store_param) return(list(pred = prediction, beta=beta_out, mu=mu_out, tau=tau_out))
  list(pred=prediction)
}

#' Variational Bayes algorithm for the DLSBP model
#'
#' The dependent logit stick-breaking process (DLSBP) model estimated through  Variational Bayes (VB).
#' 
#' @param y a vector containing the response vector
#' @param X a n x p design matrix containing the covariates
#' @param H an integer indicating the number of mixture components
#' @param maxiter an integer indicating the maximum number of iterations for the EM algorithm
#' @param prior a list containing prior hyperparameters (See details for a detailed description). If not provided, the MLE is computed.
#' @param tol a real number controlling the convergence criterium
#' @param seed set the seed, for reproducibility since the initialization is random
#' @param verbose Logical: Should the logposterior be displayed while the algorithm is running?
#' 
#' @details  The prior argument is a list which should contains the following elements
#' \itemize{
#' \item \verb{b}. A p dimensional vector containing the prior mean of the Gaussian beta coefficients
#' \item \verb{B}. A p x p matrix representing the prior covariance of the Gaussian beta coefficients.
#' \item \verb{mu_mu}. A real number representing the prior mean for the kernel means.
#' \item \verb{tau_mu}. A positive number representing the precision for the kernel means.
#' \item \verb{a_tau}, \verb{b_tau}. The hyperparameters of a gamma prior distribution for the kernel precisions.
#' }
#' 
#' @return The output is a list containing the following quantities
#' \itemize{
#' \item \verb{mu_beta}. A (H - 1) x p matrix containing the mean of the Gaussian variational distribution of the beta coefficients.
#' \item \verb{Sigma_beta}. A list of p x p matrixes containing the covariance of the Gaussian variational distribution of the beta coefficients.
#' \item \verb{mu_tilde}. A p dimensional vector containing the means of the Gaussian variational distribution for the kernel means.
#' \item \verb{tau_tilde}. A p dimensional vector containing the precisions of the Gaussian variational distribution for the kernel means.
#' \item \verb{a_tilde}, \verb{b_tilde}. The parameters of the gamma variational distribution for the kernel precisions.
#' \item \verb{cluster}. A n dimensional vector containing, for each observation, the mixture component having with the highest probability.
#' \item \verb{z}. A n x H matrix containing the probabilities of belonging to each of the mixture components.
#' \item \verb{lowerbound}. The lowerbound of the DLSBP model at convergence.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2017), Logit stick-breaking priors for Bayesian density regression, ArXiv.
#' @examples 
#' library(DLSBP)
#' data(geyser)
#' x <- geyser$waiting
#' y <- geyser$duration
#' X <- cbind(1,x*I(x<=68) + 68*I(x > 68),((x-68)*I(x>68))) # Predictors
#' p <- NCOL(X)
#' prior   <- list(b = rep(0,p),        
#'                 B = diag(100,p), 
#'                 mu_mu =0,            
#'                 tau_mu = 0.001,     
#'                 a_tau = 2,           
#'                 b_tau = 0.001)
#' fit_VB   <- DLSBP_VB(y=y, X=X, H=3, prior=prior, seed=99)
#' 
#' @export
#' 
DLSBP_VB <- function(y, X, H = 2, maxiter = 1000, prior = NULL, tol=1e-4, seed=NULL, verbose = TRUE) {
  
  if(!is.null(seed)) set.seed(seed)
  
  # Fixed quantities
  n <- length(y)
  p <- NCOL(X)
  X <- Matrix(X)
  verbose_step = ceiling(maxiter / 50)
  log2pi <- log(2*pi) # Compute the logarithm of pi.
  
  # Hyperparameters
  b <- prior$b
  mu_mu  <- prior$mu_mu
  B <- Matrix(prior$B); P <- solve(B); Pb <- P%*%b
  tau_mu <- prior$tau_mu
  a_tau <- prior$a_tau
  b_tau <- prior$b_tau
  
  # Initialization
  tau_tilde <- numeric(H)
  mu_tilde  <- numeric(H)
  tau       <- as.numeric(rep(1/diff(quantile(y,c(.25,0.75))),H))
  ltau      <- a_tilde <- b_tilde <- numeric(H) 
  xi        <- matrix(1,n,H-1)
  
  # Random initialization, with value very close to 0.
  mu_beta    <- matrix(rnorm((H-1)*p,0,0.01), H - 1, p) 
  Sigma_beta <- list(); for(h in 1:(H-1)) Sigma_beta[[h]] <- Matrix(0,p,p)
  
  # Random and uniform starting probabilities
  rho        <- matrix(runif(n*(H-1)),n,H-1)
  z          <- t(apply(rho,1,sb))
  
  # Initialization of different pieces of the lowerbound
  lower1 <- lower4 <- lower6 <- numeric(H)
  lower2 <- lower3 <- lower5 <- lower7 <- numeric(H-1)
  lowerbound <- -Inf
  
  # VB Algorithm
  for (r in 1:maxiter) {
  
    for (h in 1:H) {
      # Logistic regressions
      if (h < H) {
        
        omega <- tanh(xi[,h]/2)/(2 * xi[,h])
        omega[is.nan(omega)] <- 1/4
        
        # Parameters for beta
        Sigma_beta[[h]] <- solve(crossprod(X*sqrt(omega)) + P)
        mu_beta[h, ]    <- as.numeric(Sigma_beta[[h]] %*% (crossprod(X,rho[, h] - 1/2) + Pb))
        
        # Variational step
        xi[,h] <- sqrt(as.numeric((X %*% mu_beta[h, ])^2) + mahalanobis(X,center=0,cov=as.matrix(Sigma_beta[[h]]),inverted=TRUE))
      }
      
      # Parameters for posterior of mu and tau
      tau_tilde[h] <- tau[h]*sum(z[,h]) + tau_mu
      mu_tilde[h]  <- (tau[h]*sum(z[,h]* y) + tau_mu*mu_mu)/tau_tilde[h]
      
      residual  <- y^2 - 2*y*mu_tilde[h] + mu_tilde[h]^2 + 1/tau_tilde[h]
      a_tilde[h]<- a_tau  + sum(z[, h])/2
      b_tilde[h]   <- b_tau + sum(z[, h] * residual)/2
      tau[h]    <- a_tilde[h] / b_tilde[h]
      ltau[h]   <- digamma(a_tilde[h]) - log(b_tilde[h])
    }
    
    # Mixture probabilities
    Variational <- Variational_step(rho, y, as.matrix(X), mu_beta, mu_tilde, tau_tilde, tau, ltau)
    rho <- pmin(pmax(Variational$rho,1e-16),1-1e-16) # Added for numerical stability
    z   <- Variational$z
    
    # Lowerbound computatations
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
    if(lowerbound_new - lowerbound <tol){
      cat(paste("Convergence reached after",r,"iterations."))
      break
    }
    
    # Otherwise continue
    lowerbound <- lowerbound_new
    
    # Display status
    if (verbose) {
      if(r%%verbose_step==0) cat(paste("Lower-bound: ",round(lowerbound,4),", iteration: ", r, "\n",sep=""))
    }
  }
  
  if(r==maxiter) warning(paste("Convergence has not been reached after",r,"iterations."))
  
  # Output
  cluster <- apply(z,1,which.max)
  list(mu_beta = mu_beta, Sigma_beta=as.matrix(Sigma_beta), mu_tilde = mu_tilde, tau_tilde=tau_tilde, a_tilde=a_tilde, b_tilde=b_tilde, cluster=cluster,z=z, lowerbound=lowerbound)
}