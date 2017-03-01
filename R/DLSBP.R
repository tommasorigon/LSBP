#'@import stats
#'@importFrom Rcpp evalCpp sourceCpp
#'@importFrom BayesLogit rpg.devroye
#'@import Formula
#'@importFrom cluster clara
#'@useDynLib DLSBP

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
#'
#' @param maxiter An integer indicating the maximum number of iterations for the ECM algorithm.
#' @param tol  A real number controlling the convergence of the algorithm. It is defined as the (absolute) difference of the log-posterior between consecutive iterations that has to be reached.
#' @param method_init The initialization criterium. By default, \verb{method_init='cluster'} preallocates covariates into groups using \code{\link[cluster]{clara}}. Other available possibilities are: \verb{method_init='random'} and \verb{method_init='deterministic'}.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
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
#'
#' @param maxiter An integer indicating the maximum number of iterations for the ECM algorithm.
#' @param tol  A real number controlling the convergence of the algorithm. It is defined as the (absolute) difference of the log-posterior between consecutive iterations that has to be reached.
#' @param method_init The initialization criterium. By default, \verb{method_init='cluster'} preallocates covariates into groups using \code{\link[cluster]{clara}}. Another available possibility is \verb{method_init='random'}.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
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
#'
#' @param R An integer indicating the number of replications to be computed after the burn-in.
#' @param burn_in  An integer indicating the number of replication left out as burn-in period.
#' @param method_init The initialization criterium. By default, \verb{method_init='cluster'} preallocates covariates into groups using \code{\link[cluster]{clara}}. Other available possibilities are: \verb{method_init='random'} and \verb{method_init='deterministic'}.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
#' 
#' @export
#' 
control_Gibbs <- function(R = 5000, burn_in = 1000, method_init = "cluster") {
   if (!(method_init %in% c("cluster", "random", "deterministic"))) 
      stop("Please provide a valid initialization method")
   list(R = R, burn_in = burn_in, method_init = method_init)
}
#' Prior specification for the DLSBP model
#'
#' The prior argument is a list which contains the following elements:
#' 
#' @param p_kernel,p_mixing The dimension of the design matrices for the kernel component and the mixing component, respectively.
#' @param b_kernel A \verb{p_kernel} dimensional vector representing the prior mean for the  Gaussian kernel coefficients.
#' @param B_kernel A \verb{p_kernel x p_kernel} matrix representing the prior covariance of the  Gaussian kernel coefficients.
#' @param b_mixing A \verb{p_mixing} dimensional vector containing the prior mean of the Gaussian mixing coefficients
#' @param B_mixing A \verb{p_mixing x p_mixing} matrix representing the prior covariance  of the Gaussian mixing coefficients.
#' @param a_tau,b_tau The hyperparameters of a Gamma prior distribution for the kernel precision.
#' 
#' @return The function returns a list having the same entries provided as argument. Missing arguments are filled with default values.
#' 
#' @examples 
#' data(cars)
#' prior  <- prior_DLSBP(p_kernel=1, p_mixing=2, a_tau=1.5 ,b_tau=1.5)
#' fit_em <- DLSBP_ECM(dist ~ 1 | speed,data=cars, H=4, prior=prior)
#' 
#' @export

prior_DLSBP <- function(p_kernel, p_mixing, b_kernel = rep(0, p_kernel), B_kernel = diag(10^6, p_kernel), 
   b_mixing = rep(0, p_mixing), B_mixing = diag(10^4, p_mixing), a_tau = 0.1, b_tau = 0.1) {
   if (missing(p_kernel) | missing(p_mixing)) 
      stop("Please provide a value for p_kernel and p_mixing")
   prior <- list(b_kernel = b_kernel, B_kernel = B_kernel, b_mixing = b_mixing, B_mixing = B_mixing, 
      a_tau = a_tau, b_tau = b_tau)
   return(prior)
}



#' ECM algorithm for the DLSBP model
#'
#' The dependent logit stick-breaking process (DLSBP) model is estimated using the E(C)M algorithm, which provides the posterior mode.
#' 
#' @param Formula An object of class \code{\link[Formula]{Formula}}: a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#' @param data A data frame containing the variables of \verb{Formula}.
#' @param H An integer indicating the number of mixture components.
#' @param prior A list of prior hyperparameters as returned by \code{\link[DLSBP]{prior_DLSBP}}. If missing, default prior values are used.
#' @param control A list as returned by \code{\link[DLSBP]{control_ECM}}.
#' @param verbose A logical value indicating whether additional information should be displayed while the algorithm is running.
#' 
#' @details 
#' The \verb{Formula} specification contains the response, separated from the covariates with the symbol '\verb{~}', and two sets of covariates. The latters are separated by the symbol '\verb{|}', indicating the kernel covariates and the mixing covariates, respectively. For example, one could specify \verb{y ~ x1 + x2 | x3 + x4}. NOTE: if the second set of covariates is omitted it is assumed that the two sets are the same.
#' 
#' If \verb{offsets} or \verb{weights} are provided in the \verb{Formula} they will be ignored in the current version.
# 
#' A \verb{predict} method is available and described at \code{\link[DLSBP]{predict.DLSBP_ECM}}.
#' 
#' 
#' @return The output is an object of class '\verb{DLSBP_ECM}' containing the following quantities:
#' \itemize{
#' \item \verb{param}. A list containing the MAP (Maximum A Posteriori), for each set of coefficients: \verb{beta_mixing,beta_kernel,tau}.
#' \item \verb{cluster}. A n dimensional vector containing, for each observation, the mixture component having with the highest probability.
#' \item \verb{z}. A \verb{n x H} matrix containing the probabilities of belonging to each of the mixture components, where \verb{n} denotes the number of observations.
#' \item \verb{logposterior}. The log-posterior of the model at convergence.
#' \item \verb{call}. The input Formula 
#' \item \verb{data}. The input data frame.
#' \item \verb{control}. The control list provided as input.
#' \item \verb{H}. The input number of mixture components.
#' \item \verb{prior}. The input prior hyperparameters.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2017), Logit stick-breaking priors for Bayesian density regression, ArXiv.
#' @examples 
#' data(cars)
#' 
#' # A model with constant kernels
#' fit_em <- DLSBP_ECM(dist ~  1 | speed, data=cars, H=4)
#' plot(cars) 
#' lines(cars$speed,predict(fit_em))
#' 
#' # A model with linear kernels
#' fit_em <- DLSBP_ECM(dist ~ speed | speed, data=cars, H=2)
#' plot(cars) 
#' lines(cars$speed,predict(fit_em))
#' 
#' @export
#' 

DLSBP_ECM <- function(Formula, data = NULL, H, prior, control = control_ECM(), verbose = TRUE) {
   
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
      prior <- prior_DLSBP(p_kernel, p_mixing)
   
   # Settings
   verbose_step = 100
   
   if (NCOL(X1) > 1) {
      out <- DLSBP_ECM_multi(y = y, X1 = X1, X2 = X2, H = H, prior = prior, maxiter = control$maxiter, 
         tol = control$tol, method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   } else {
      out <- DLSBP_ECM_univ(y = y, X = X2, H = H, prior = prior, maxiter = control$maxiter, tol = control$tol, 
         method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   }
   attr(out, "class") <- "DLSBP_ECM"
   out$call <- Formula
   out$data <- list(y = y, X1 = X1, X2 = X2)
   out$control <- control
   out$H <- H
   out$prior <- prior
   return(out)
}

# #' @export
# print.DLSBP_ECM <- function(x) {
#    cat(paste("Convergence reached at logposterior: ", x$logposterior, ".\n", sep = ""))
#    print(lapply(x$param, function(y) round(y, 2)))
# }


#' Gibbs sampling algorithm for the DLSBP model
#'
#' The dependent logit stick-breaking process (DLSBP) model estimated through the Gibbs sampling.
#' 
#' @param Formula An object of class \code{\link[Formula]{Formula}}: a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#' @param data A data frame containing the variables of \verb{Formula}.
#' @param H An integer indicating the number of mixture components.
#' @param prior A list of prior hyperparameters as returned by \code{\link[DLSBP]{prior_DLSBP}}. If missing, default prior values are used.
#' @param control A list as returned by \code{\link[DLSBP]{control_Gibbs}}.
#' @param verbose A logical value indicating whether additional information should be displayed while the algorithm is running.
#' 
#' @details 
#' The \verb{Formula} specification contains the response, separated from the covariates with the symbol '\verb{~}', and two sets of covariates. The latters are separated by the symbol '\verb{|}', indicating the kernel covariates and the mixing covariates, respectively. For example, one could specify \verb{y ~ x1 + x2 | x3 + x4}. NOTE: if the second set of covariates is omitted it is assumed that the two sets are the same.
#' 
#' If \verb{offsets} or \verb{weights} are provided in the \verb{Formula} they will be ignored in the current version.
# 
#' A \verb{predict} method is available and described at \code{\link[DLSBP]{predict.DLSBP_Gibbs}}.
#' 
#' 
#' @return The output is an object of class '\verb{DLSBP_Gibbs}' containing the following quantities:
#' \itemize{
#' \item \verb{param}. A list containing MCMC replications for each set of coefficients: \verb{beta_mixing, beta_kernel, tau}.
#' \item \verb{logposterior}. The log-posterior of the model at each MCMC iteration.
#' \item \verb{call}. The input Formula.
#' \item \verb{data}. The input data frame.
#' \item \verb{control}. The control list provided as input.
#' \item \verb{H}. The input number of mixture components.
#' \item \verb{prior}. The input prior hyperparameters.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2017), Logit stick-breaking priors for Bayesian density regression, ArXiv.
#' @examples 
#' \dontrun{
#' data(cars)
#' 
#' # A model with constant kernels
#' fit_gibbs <- DLSBP_Gibbs(dist ~  1 | speed, data=cars, H=4)
#' plot(cars) 
#' lines(cars$speed,colMeans(predict(fit_gibbs)))
#' 
#' # A model with linear kernels
#' fit_gibbs <- DLSBP_Gibbs(dist ~ speed | speed, data=cars, H=2)
#' plot(cars) 
#' lines(cars$speed,colMeans(predict(fit_gibbs)))
#' }
#' 
#' @export
#' 

DLSBP_Gibbs <- function(Formula, data = NULL, H , prior, control = control_Gibbs(), verbose = TRUE) {
   
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
      prior <- prior_DLSBP(p_kernel, p_mixing)
   
   if (NCOL(X1) > 1) {
      out <- DLSBP_Gibbs_multi(y = y, X1 = X1, X2 = X2, H = H, prior = prior, R = control$R, burn_in = control$burn_in, 
         method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   } else {
      out <- DLSBP_Gibbs_univ(y = y, X = X2, H = H, prior = prior, R = control$R, burn_in = control$burn_in, 
         method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   }
   attr(out, "class") <- "DLSBP_Gibbs"
   out$call <- Formula
   out$data <- list(y = y, X1 = X1, X2 = X2)
   out$control <- control
   out$H <- H
   out$prior <- prior
   return(out)
}

# #' @export
# print.DLSBP_Gibbs <- function(x, plot = FALSE) {
#    param <- x$logposterior
#    if (plot == TRUE) 
#       plot(param)
#    print(summary(param))
# }



#' Variational Bayes algorithm for the DLSBP model
#'
#' The dependent logit stick-breaking process (DLSBP) model is estimated through a Variational Bayes (VB) algorithm.
#' 
#' @param Formula An object of class \code{\link[Formula]{Formula}}: a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#' @param data A data frame containing the variables of \verb{Formula}.
#' @param H An integer indicating the number of mixture components.
#' @param prior A list of prior hyperparameters as returned by \code{\link[DLSBP]{prior_DLSBP}}. If missing, default prior values are used.
#' @param control A list as returned by \code{\link[DLSBP]{control_VB}}.
#' @param verbose A logical value indicating whether additional information should be displayed while the algorithm is running.
#' 
#' @details 
#' The \verb{Formula} specification contains the response, separated from the covariates with the symbol '\verb{~}', and two sets of covariates. The latters are separated by the symbol '\verb{|}', indicating the kernel covariates and the mixing covariates, respectively. For example, one could specify \verb{y ~ x1 + x2 | x3 + x4}. NOTE: if the second set of covariates is omitted it is assumed that the two sets are the same.
#' 
#' If \verb{offsets} or \verb{weights} are provided in the \verb{Formula} they will be ignored in the current version.
# 
#' A \verb{predict} method is available and described at \code{\link[DLSBP]{predict.DLSBP_VB}}.
#' 
#' 
#' @return The output is an object of class '\verb{DLSBP_VB}' containing the following quantities:
#' \itemize{
#' \item \verb{param}. A list containing the parameters for the variational approximation of each distribution: \verb{mu_mixing, Sigma_mixing, mu_kernel, Sigma_kernel, a_tilde, b_tilde}.
#' \item \verb{cluster}. A n dimensional vector containing, for each observation, the mixture component having with the highest probability.
#' \item \verb{z}. A \verb{n x H} matrix containing the probabilities of belonging to each of the mixture components, where \verb{n} denotes the number of observations.
#' \item \verb{logposterior}. The log-posterior of the model at convergence.
#' \item \verb{call}. The input Formula 
#' \item \verb{data}. The input data frame.
#' \item \verb{control}. The control list provided as input.
#' \item \verb{H}. The input number of mixture components.
#' \item \verb{prior}. The input prior hyperparameters.
#' }
#' 
#' @references Rigon, T. and Durante, D., (2017), Logit stick-breaking priors for Bayesian density regression, ArXiv.
#' @examples 
#' data(cars)
#' 
#' # A model with constant kernels
#' fit_vb <- DLSBP_VB(dist ~  1 | speed, data=cars, H=4)
#' plot(cars) 
#' lines(cars$speed,colMeans(predict(fit_vb)))
#' 
#' # A model with linear kernels
#' fit_vb <- DLSBP_VB(dist ~ speed | speed, data=cars, H=2)
#' plot(cars) 
#' lines(cars$speed,colMeans(predict(fit_vb)))
#' 
#' @export
#' 

DLSBP_VB <- function(Formula, data = NULL, H , prior, control = control_VB(), verbose = TRUE) {
   
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
      prior <- prior_DLSBP(p_kernel, p_mixing)
   
   # Settings
   verbose_step = 100
   
   if (NCOL(X1) > 1) {
      out <- DLSBP_VB_multi(y = y, X1 = X1, X2 = X2, H = H, prior = prior, maxiter = control$maxiter, 
         tol = control$tol, method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   } else {
      out <- DLSBP_VB_univ(y = y, X = X2, H = H, prior = prior, maxiter = control$maxiter, tol = control$tol, 
         method_init = control$method_init, verbose = verbose, verbose_step = verbose_step)
   }
   attr(out, "class") <- "DLSBP_VB"
   out$call <- Formula
   out$data <- list(y = y, X1 = X1, X2 = X2)
   out$control <- control
   out$H <- H
   out$prior <- prior
   return(out)
}


# #' @export
# print.DLSBP_VB <- function(x) {
#    cat(paste("Convergence reached at lower-bound: ", x$lowerbound, ".\n", sep = ""))
#    print(lapply(x$param, function(y) round(y, 2)))
# }
