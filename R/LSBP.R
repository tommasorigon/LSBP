#'@import stats
#'@importFrom mvtnorm rmvnorm
#'@importFrom Rcpp evalCpp sourceCpp
#'@importFrom BayesLogit rpg.devroye
#'@import Formula
#'@importFrom cluster clara
#'@useDynLib LSBP
#'

sb <- function(nu) {
   nu <- c(nu, 1)
   H <- length(nu)
   prob <- nu * c(1, cumprod(1 - nu[-H]))
   return(prob)
}

ldmvnorm <- function(x, mean, eig) {
   eigsqrt <- sqrt(eig$values)
   Sigma1 <- crossprod(t(eig$vectors)/eigsqrt)
   x <- x - mean
   -sum(log(eigsqrt)) - 0.5 * sum(x %*% Sigma1 * x) - 0.5 * log(2 * pi) * length(x)
}

ldet <- function(x) {
  if(!is.matrix(x)) return(log(x))
  determinant(x,logarithm = TRUE)$modulus
}

#' Control parameters for the ECM algorithm
#'
#' The \code{control_ECM} function can be used for specifying the technical settings (i.e. the maximum number of iterations, the tolerance level, and the initialization method), of the \code{\link[LSBP]{LSBP_ECM}} function.
#'
#' @param maxiter An integer indicating the maximum number of iterations for the \code{\link[LSBP]{LSBP_ECM}} algorithm.
#' @param tol  A real number controlling the convergence of the algorithm. The \code{\link[LSBP]{LSBP_ECM}} algorithm stops when the difference between consecutive values of the log-posterior is smaller than \code{tol}.
#' @param method_init The initialization method. The default \code{method_init='cluster'} partitions the covariates  using the \code{\link[cluster]{clara}} clustering algorithm. Other available options are: \code{method_init='random'} and \code{method_init='deterministic'}.
#' 
#' @return The inputs are converted into a \code{list}. Missing arguments are filled with default values.
#' 
#' @export
#' 

control_ECM <- function(maxiter = 10000, tol = 0.001, method_init = "cluster") {
   if (!(method_init %in% c("cluster", "random", "deterministic"))) 
      stop("Please provide a valid initialization method")
   list(maxiter = maxiter, tol = tol, method_init = method_init)
}

#' Control parameters for the VB algorithm
#'
#' The \code{control_VB} function can be used for specifying the technical settings (i.e. the maximum number of iterations, the tolerance level, and the initialization method), of the \code{\link[LSBP]{LSBP_VB}} main function.
#'
#'
#' @param maxiter An integer indicating the maximum number of iterations for the VB algorithm.
#' @param tol  A real number controlling the convergence of the algorithm. The \code{\link[LSBP]{LSBP_VB}} algorithm stops when the difference between consecutive values of the evidence lower bound (ELBO) is smaller than \code{tol}.
#' @param method_init The initialization method. The default \code{method_init='cluster'} partitions the covariates  using the \code{\link[cluster]{clara}} clustering algorithm. Another available option is \code{method_init='random'}.
#' 
#' @return The inputs are converted into a \code{list}. Missing arguments are filled with default values.
#' 
#' @export
#' 

control_VB <- function(maxiter = 10000, tol = 0.01, method_init = "cluster") {
   if (!(method_init %in% c("cluster", "random"))) 
      stop("Please provide a valid initialization method")
   list(maxiter = maxiter, tol = tol, method_init = method_init)
}

#' Control parameters for the Gibbs sampling algorithm
#'
#' The \code{control_Gibbs} function can be used for specifying the technical settings (i.e. the number of MCMC iterations, the burn-in, and the initialization method), of the \code{\link[LSBP]{LSBP_Gibbs}} function.
#'
#' @param R An integer indicating the number of MCMC iterations to be computed after the \code{burn-in}.
#' @param burn_in  An integer indicating the number of MCMC iterations discarded as burn-in period.
#' @param method_init The initialization method. The default \code{method_init='cluster'} partitions the covariates  using the \code{\link[cluster]{clara}} clustering algorithm.  Other available options are: \code{method_init='random'} and \code{method_init='deterministic'}.
#' 
#' @return The inputs are converted into a \code{list}. Missing arguments are filled with default values.
#' 
#' @export
#' 

control_Gibbs <- function(R = 5000, burn_in = 1000, method_init = "cluster") {
   if (!(method_init %in% c("cluster", "random", "deterministic"))) 
      stop("Please provide a valid initialization method")
   list(R = R, burn_in = burn_in, method_init = method_init)
}

#' Prior specification for the LSBP model
#'
#' This auxiliary function can be used for specifying the prior hyperparameters in the \code{\link[LSBP]{LSBP_Gibbs}}, \code{\link[LSBP]{LSBP_ECM}}, \code{\link[LSBP]{LSBP_VB}} main functions.
#' 
#' @param p_kernel,p_mixing The dimension of the design matrices for the kernel component and the mixing component, respectively.
#' @param b_kernel A \code{p_kernel} dimensional vector representing the prior mean for the  Gaussian kernel coefficients.
#' @param B_kernel A \code{p_kernel x p_kernel} matrix representing the prior covariance of the  Gaussian kernel coefficients.
#' @param b_mixing A \code{p_mixing} dimensional vector containing the prior mean of the Gaussian mixing coefficients
#' @param B_mixing A \code{p_mixing x p_mixing} matrix representing the prior covariance  of the Gaussian mixing coefficients.
#' @param a_tau,b_tau The hyperparameters of a Gamma prior distribution for the kernel precision.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values, although this is NOT recommended in general.
#' 
#' @examples 
#' \dontrun{
#' data(cars)
#' prior  <- prior_LSBP(p_kernel=1, p_mixing=2, a_tau=1.5 ,b_tau=1.5)
#' fit_em <- LSBP_ECM(dist ~ 1 | speed, data=cars, H=4, prior=prior)
#' }
#' @export
#' 

prior_LSBP <- function(p_kernel, p_mixing, b_kernel = rep(0, p_kernel), B_kernel = diag(10^6, p_kernel), 
   b_mixing = rep(0, p_mixing), B_mixing = diag(10^4, p_mixing), a_tau = 0.1, b_tau = 0.1) {
   if (missing(p_kernel) | missing(p_mixing)) 
      stop("Please provide a value for p_kernel and p_mixing")
   prior <- list(b_kernel = b_kernel, B_kernel = B_kernel, b_mixing = b_mixing, B_mixing = B_mixing, 
      a_tau = a_tau, b_tau = b_tau)
   return(prior)
}



#' ECM algorithm for the LSBP model
#'
#' This function is an implementation of the expectation maximization Algorithm 2 in Rigon, T. and Durante, D. (2020). 
#' 
#' @param Formula An object of class \code{\link[Formula]{Formula}}: a symbolic description of the model to be estimated. The details of model specification are given under "Details".
#' @param data A data frame containing the variables described in \code{Formula}. The data frame must be provided.
#' @param H An integer indicating the number of mixture components.
#' @param prior A list of prior hyperparameters as returned by \code{\link[LSBP]{prior_LSBP}}. If missing, default prior values are used, although this is NOT recommended.
#' @param control A list as returned by \code{\link[LSBP]{control_ECM}}.
#' @param verbose A logical value indicating whether additional information should be displayed while the algorithm is running.
#' 
#' @details 
#' The \code{Formula} specification contains the response \code{y}, separated from the covariates with the symbol '\code{~}', and two sets of covariates. The latters are separated by the symbol '\code{|}', indicating the kernel covariates and the mixing covariates, respectively. For example, one could specify \code{y ~ x1 + x2 | x3 + x4}. NOTE: if the second set of covariates is omitted, then it is implicitely assumed that the two sets are the same.
#' 
#' If \code{offsets} or \code{weights} are provided in the \code{Formula} they will be ignored in the current version.
# 
#' A \code{predict} method is available and described at \code{\link[LSBP]{predict.LSBP_ECM}}.
#' 
#' 
#' @return The output is an object of class "\code{LSBP_ECM}" containing the following quantities:
#' \itemize{
#' \item \code{param}. A list containing the maximum a posteriori, for each set of coefficients: \code{beta_mixing}, \code{beta_kernel}, and  \code{tau}.
#' \item \code{cluster}. A \code{n} dimensional vector containing, for each observation, the mixture component having with the highest probability.
#' \item \code{z}. A \code{n x H} matrix containing the probabilities of belonging to each of the mixture components, where \code{n} denotes the number of observations.
#' \item \code{logposterior}. The \code{log-posterior} of the model at convergence. NOTE: the \code{log-posterior} is reported up to an additive constant.
#' \item \code{call}. The input Formula.
#' \item \code{data}. The input data frame.
#' \item \code{control}. The control list provided as input.
#' \item \code{H}. The input number of mixture components.
#' \item \code{prior}. The input prior hyperparameters.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2020), Tractable Bayesian density regression via logit stick-breaking priors. Journal of Statistical Planning and Inference.
#' @examples 
#' \dontrun{
#' data(cars)
#' 
#' # A model with constant kernels
#' fit_em <- LSBP_ECM(dist ~  1 | speed, data=cars, H=4)
#' plot(cars) 
#' lines(cars$speed,predict(fit_em))
#' 
#' # A model with linear kernels
#' fit_em <- LSBP_ECM(dist ~ speed | speed, data=cars, H=2)
#' plot(cars) 
#' lines(cars$speed,predict(fit_em))
#' }
#' @export
#' 

LSBP_ECM <- function(Formula, data, H, prior, control = control_ECM(), verbose = TRUE) {
   
  if(is.null(data)) stop("The data argument can not be NULL and must be specified")
   Formula <- as.Formula(Formula)
   y <- model.frame(Formula, data = data, rhs = 0)[, 1]
   X1 <- model.matrix(Formula, data = data, rhs = 1)
   if (length(Formula)[2] == 1) {
      X2 <- X1
   } else {
      X2 <- model.matrix(Formula, data = data, rhs = 2)
   }
   
   p_kernel <- NCOL(X1)
   p_mixing <- NCOL(X2)
   
   if (missing(prior)) 
      prior <- prior_LSBP(p_kernel, p_mixing)
   
   if(any(c(p_kernel != ncol(prior$B_kernel),p_mixing != ncol(prior$B_mixing)))) stop("The dimension of the prior distribution must coincide with the dimension originated from the Formula")
   
   # Settings
   verbose_step = 100
   
   if (NCOL(X1) > 1) {
      out <- LSBP_ECM_multi(y = y, X1 = X1, X2 = X2, H = H, prior = prior, maxiter = control$maxiter, 
         tol = control$tol, method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   } else {
      out <- LSBP_ECM_univ(y = y, X = X2, H = H, prior = prior, maxiter = control$maxiter, tol = control$tol, 
         method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   }
   attr(out, "class") <- "LSBP_ECM"
   out$call <- Formula
   out$data <- list(y = y, X1 = X1, X2 = X2)
   out$control <- control
   out$H <- H
   out$prior <- prior
   return(out)
}

#' Gibbs sampling algorithm for the LSBP model
#'
#' This function is an implementation of the Gibbs sampling Algorithm 1 in Rigon, T. and Durante, D. (2020). 
#'  
#' @param Formula An object of class \code{\link[Formula]{Formula}}: a symbolic description of the model to be estimated. The details of model specification are given under "Details".
#' @param data A data frame containing the variables described in \code{Formula}. The data frame must be provided.
#' @param H An integer indicating the number of mixture components.
#' @param prior A list of prior hyperparameters as returned by \code{\link[LSBP]{prior_LSBP}}. If missing, default prior values are used, although this is NOT recommended.
#' @param control A list as returned by \code{\link[LSBP]{control_Gibbs}}.
#' @param verbose A logical value indicating whether additional information should be displayed while the algorithm is running.
#' 
#' @details 
#' The \code{Formula} specification contains the response \code{y}, separated from the covariates with the symbol '\code{~}', and two sets of covariates. The latters are separated by the symbol '\code{|}', indicating the kernel covariates and the mixing covariates, respectively. For example, one could specify \code{y ~ x1 + x2 | x3 + x4}. NOTE: if the second set of covariates is omitted, then it is implicitely assumed that the two sets are the same.
#' 
#' If \code{offsets} or \code{weights} are provided in \code{Formula}, they will be IGNORED in the current version.
# 
#' A \code{predict} method is available and described at \code{\link[LSBP]{predict.LSBP_Gibbs}}.
#' 
#' 
#' @return The output is an object of class "\code{LSBP_Gibbs}" containing the following quantities:
#' \itemize{
#' \item \code{param}. A list containing MCMC replications for each set of coefficients: \code{beta_mixing, beta_kernel, tau}.
#' \item \code{logposterior}. The \code{log-posterior} of the model at each MCMC iteration.  NOTE: the \code{log-posterior} is reported up to an additive constant.
#' \item \code{call}. The input Formula.
#' \item \code{data}. The input data frame.
#' \item \code{control}. The control list provided as input.
#' \item \code{H}. The input number of mixture components.
#' \item \code{prior}. The input prior hyperparameters.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2020), Tractable Bayesian density regression via logit stick-breaking priors. Journal of Statistical Planning and Inference.
#' @examples 
#' \dontrun{
#' data(cars)
#' 
#' # A model with constant kernels
#' fit_gibbs <- LSBP_Gibbs(dist ~  1 | speed, data=cars, H=4)
#' plot(cars) 
#' lines(cars$speed,colMeans(predict(fit_gibbs))) # Posterior mean
#' 
#' # A model with linear kernels
#' fit_gibbs <- LSBP_Gibbs(dist ~ speed | speed, data=cars, H=2)
#' plot(cars) 
#' lines(cars$speed,colMeans(predict(fit_gibbs))) # Posterior mean
#' }
#' 
#' @export
#' 
LSBP_Gibbs <- function(Formula, data, H , prior, control = control_Gibbs(), verbose = TRUE) {
   
   if(is.null(data)) stop("The data argument can not be NULL and must be specified")
   Formula <- as.Formula(Formula)
   y <- model.frame(Formula, data = data, rhs = 0)[, 1]
   X1 <- model.matrix(Formula, data = data, rhs = 1)
   if (length(Formula)[2] == 1) {
      X2 <- X1
   } else {
      X2 <- model.matrix(Formula, data = data, rhs = 2)
   }
   
   p_kernel <- NCOL(X1)
   p_mixing <- NCOL(X2)
   
   # Settings
   verbose_step = ceiling(control$R/50)
   
   if (missing(prior)) 
      prior <- prior_LSBP(p_kernel, p_mixing)
   
   if(any(c(p_kernel != ncol(prior$B_kernel),p_mixing != ncol(prior$B_mixing)))) stop("The dimension of the prior distribution must coincide with the dimension originated from the Formula")
   
   
   if (NCOL(X1) > 1) {
      out <- LSBP_Gibbs_multi(y = y, X1 = X1, X2 = X2, H = H, prior = prior, R = control$R, burn_in = control$burn_in, 
         method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   } else {
      out <- LSBP_Gibbs_univ(y = y, X = X2, H = H, prior = prior, R = control$R, burn_in = control$burn_in, 
         method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   }
   attr(out, "class") <- "LSBP_Gibbs"
   out$call <- Formula
   out$data <- list(y = y, X1 = X1, X2 = X2)
   out$control <- control
   out$H <- H
   out$prior <- prior
   return(out)
}

#' Variational Bayes algorithm for the LSBP model
#'
#' This function is an implementation of the variational Bayes Algorithm 3 in Rigon, T. and Durante, D. (2020). 
#' 
#' @param Formula An object of class \code{\link[Formula]{Formula}}: a symbolic description of the model to be estimated. The details of model specification are given under "Details".
#' @param data A data frame containing the variables described in \code{Formula}. The data frame must be provided.
#' @param H An integer indicating the number of mixture components.
#' @param prior A list of prior hyperparameters as returned by \code{\link[LSBP]{prior_LSBP}}. If missing, default prior values are used, although this is NOT recommended.
#' @param control A list as returned by \code{\link[LSBP]{control_VB}}.
#' @param verbose A logical value indicating whether additional information should be displayed while the algorithm is running.
#' 
#' @details 
#' The \code{Formula} specification contains the response \code{y}, separated from the covariates with the symbol '\code{~}', and two sets of covariates. The latters are separated by the symbol '\code{|}', indicating the kernel covariates and the mixing covariates, respectively. For example, one could specify \code{y ~ x1 + x2 | x3 + x4}. NOTE: if the second set of covariates is omitted, then it is implicitely assumed that the two sets are the same.
#' 
#' If \code{offsets} or \code{weights} are provided in \code{Formula}, they will be IGNORED in the current version.
# 
#' A \code{predict} method is available and described at \code{\link[LSBP]{predict.LSBP_VB}}.
#' 
#' 
#' @return The output is an object of class "\code{LSBP_VB}" containing the following quantities:
#' \itemize{
#' \item \code{param}. A list containing the parameters for the variational approximation of each distribution: \code{mu_mixing}, \code{Sigma_mixing}, \code{mu_kernel}, \code{Sigma_kernel}, \code{a_tilde}, \code{b_tilde}.
#' \item \code{cluster}. A \code{n} dimensional vector containing, for each observation, the mixture component having with the highest probability.
#' \item \code{z}. A \code{n x H} matrix containing the probabilities of belonging to each of the mixture components, where \code{n} denotes the number of observations.
#' \item \code{lowerbound}. The \code{lowerbound} is the evidence lower bound (ELBO) of the model at convergence.  NOTE: the \code{lowerbound} is reported up to an additive constant.
#' \item \code{call}. The input Formula.
#' \item \code{data}. The input data frame.
#' \item \code{control}. The control list provided as input.
#' \item \code{H}. The input number of mixture components.
#' \item \code{prior}. The input prior hyperparameters.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2020), Tractable Bayesian density regression via logit stick-breaking priors. Journal of Statistical Planning and Inference.
#' 
#' @examples 
#' data(cars)
#' 
#' # A model with constant kernels
#' fit_vb <- LSBP_VB(dist ~  1 | speed, data=cars, H=4)
#' plot(cars) 
#' lines(cars$speed,colMeans(predict(fit_vb)))
#' 
#' # A model with linear kernels
#' fit_vb <- LSBP_VB(dist ~ speed | speed, data=cars, H=2)
#' plot(cars) 
#' lines(cars$speed,colMeans(predict(fit_vb)))
#' 
#' @export
#' 

LSBP_VB <- function(Formula, data, H , prior, control = control_VB(), verbose = TRUE) {
   
   if(is.null(data)) stop("The data argument can not be NULL and must be specified")
   Formula <- as.Formula(Formula)
   y <- model.frame(Formula, data = data, rhs = 0)[, 1]
   X1 <- model.matrix(Formula, data = data, rhs = 1)
   if (length(Formula)[2] == 1) {
      X2 <- X1
   } else {
      X2 <- model.matrix(Formula, data = data, rhs = 2)
   }
   
   p_kernel <- NCOL(X1)
   p_mixing <- NCOL(X2)
   
   if (missing(prior)) 
      prior <- prior_LSBP(p_kernel, p_mixing)
   
   if(any(c(p_kernel != ncol(prior$B_kernel),p_mixing != ncol(prior$B_mixing)))) stop("The dimension of the prior distribution must coincide with that originated from the Formula")
   
   # Settings
   verbose_step = 100
   
   if (NCOL(X1) > 1) {
      out <- LSBP_VB_multi(y = y, X1 = X1, X2 = X2, H = H, prior = prior, maxiter = control$maxiter, 
         tol = control$tol, method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   } else {
      out <- LSBP_VB_univ(y = y, X = X2, H = H, prior = prior, maxiter = control$maxiter, tol = control$tol, 
         method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   }
   attr(out, "class") <- "LSBP_VB"
   out$call <- Formula
   out$data <- list(y = y, X1 = X1, X2 = X2)
   out$control <- control
   out$H <- H
   out$prior <- prior
   return(out)
}