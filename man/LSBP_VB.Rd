% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LSBP.R
\name{LSBP_VB}
\alias{LSBP_VB}
\title{Variational Bayes algorithm for the LSBP model}
\usage{
LSBP_VB(Formula, data, H, prior, control = control_VB(), verbose = TRUE)
}
\arguments{
\item{Formula}{An object of class \code{\link[Formula]{Formula}}: a symbolic description of the model to be estimated. The details of model specification are given under "Details".}

\item{data}{A data frame containing the variables described in \code{Formula}. The data frame must be provided.}

\item{H}{An integer indicating the number of mixture components.}

\item{prior}{A list of prior hyperparameters as returned by \code{\link[LSBP]{prior_LSBP}}. If missing, default prior values are used, although this is NOT recommended.}

\item{control}{A list as returned by \code{\link[LSBP]{control_VB}}.}

\item{verbose}{A logical value indicating whether additional information should be displayed while the algorithm is running.}
}
\value{
The output is an object of class "\code{LSBP_VB}" containing the following quantities:
\itemize{
\item \code{param}. A list containing the parameters for the variational approximation of each distribution: \code{mu_mixing}, \code{Sigma_mixing}, \code{mu_kernel}, \code{Sigma_kernel}, \code{a_tilde}, \code{b_tilde}.
\item \code{cluster}. A \code{n} dimensional vector containing, for each observation, the mixture component having with the highest probability.
\item \code{z}. A \code{n x H} matrix containing the probabilities of belonging to each of the mixture components, where \code{n} denotes the number of observations.
\item \code{lowerbound}. The \code{lowerbound} is the evidence lower bound (ELBO) of the model at convergence.  NOTE: the \code{lowerbound} is reported up to an additive constant.
\item \code{call}. The input Formula.
\item \code{data}. The input data frame.
\item \code{control}. The control list provided as input.
\item \code{H}. The input number of mixture components.
\item \code{prior}. The input prior hyperparameters.
}
}
\description{
This function is an implementation of the variational Bayes Algorithm 3 in Rigon, T. and Durante, D. (2020).
}
\details{
The \code{Formula} specification contains the response \code{y}, separated from the covariates with the symbol '\code{~}', and two sets of covariates. The latters are separated by the symbol '\code{|}', indicating the kernel covariates and the mixing covariates, respectively. For example, one could specify \code{y ~ x1 + x2 | x3 + x4}. NOTE: if the second set of covariates is omitted, then it is implicitely assumed that the two sets are the same.

If \code{offsets} or \code{weights} are provided in \code{Formula}, they will be IGNORED in the current version.
A \code{predict} method is available and described at \code{\link[LSBP]{predict.LSBP_VB}}.
}
\examples{
data(cars)

# A model with constant kernels
fit_vb <- LSBP_VB(dist ~ 1 | speed, data = cars, H = 4)
plot(cars)
lines(cars$speed, colMeans(predict(fit_vb)))

# A model with linear kernels
fit_vb <- LSBP_VB(dist ~ speed | speed, data = cars, H = 2)
plot(cars)
lines(cars$speed, colMeans(predict(fit_vb)))

}
\references{
Rigon, T. and Durante, D., (2020), Tractable Bayesian density regression via logit stick-breaking priors. Journal of Statistical Planning and Inference.
}
