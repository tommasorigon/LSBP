#' Predict method for the LSBP
#'
#'
#' Predict method for a LSBP model, estimated using the \code{\link[LSBP]{LSBP_ECM}} function.
#'
#' @param object An object of class \code{\link[LSBP]{LSBP_ECM}}.
#' @param type String indicating the type of prediction. The available options are: \code{type="mean"}, \code{type="variance"}, or \code{type="cdf"}. See "Details".
#' @param newdata A new data frame containing the same variables declared in \code{Formula}. If missing, the dataset provided for estimation is used.
#' @param threshold Only needed if \code{type="cdf"} is selected. See "Details".
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method \code{predict.LSBP_ECM} produces predicted values, obtained by evaluating the conditional mean (if \code{type="mean"}), the conditional variance (if \code{type="variance"}) or the conditional cumulative distribution function (if \code{type="cdf"}) at a given \code{threshold}, after plugging-in the maximum a posteriori, and using the observations contained in the \code{newdata} data frame.
#'
#' @export
#'

predict.LSBP_ECM <- function(object, type = "mean", newdata = NULL, threshold = NULL, ...) {
  if (!(type %in% c("mean", "variance", "cdf"))) {
    stop("Please provide a valid predictive type")
  }
  if (is.null(newdata)) {
    X1 <- object$data$X1
    X2 <- object$data$X2
  } else {
    f <- object$call
    X1 <- model.matrix(f, data = newdata, rhs = 1)
    if (length(f)[2] == 1) {
      X2 <- X1
    } else {
      X2 <- model.matrix(f, data = newdata, rhs = 2)
    }
  }
  if (type == "mean") {
    if (NCOL(X1) > 1) {
      pred <- pred_mean_multi(X1, X2, object$param$beta_mixing, object$param$beta_kernel)
    } else {
      # Univariate case
      pred <- pred_mean(X2, object$param$beta_mixing, object$param$beta_kernel)
    }
  }

  if (type == "variance") {
    if (NCOL(X1) > 1) {
      pred <- pred_var_multi(X1, X2, object$param$beta_mixing, object$param$beta_kernel, object$param$tau)
    } else {
      # Univariate case
      pred <- pred_var(X2, object$param$beta_mixing, object$param$beta_kernel, object$param$tau)
    }
  }

  if (type == "cdf") {
    if (NCOL(X1) > 1) {
      pred <- pred_cdf_multi(X1, X2, object$param$beta_mixing, object$param$beta_kernel, object$param$tau, threshold)
    } else {
      # Univariate case
      pred <- pred_cdf(X2, object$param$beta_mixing, object$param$beta_kernel, object$param$tau, threshold)
    }
  }
  return(as.numeric(pred))
}

#' Predict method for the LSBP
#'
#'
#' Predict method for a LSBP model estimated using the \code{\link[LSBP]{LSBP_Gibbs}} function.
#'
#' @param object An object of class \code{\link[LSBP]{LSBP_Gibbs}}.
#' @param type String indicating the type of prediction. The available options are \code{type="mean"},\code{type="predictive"}, \code{type="variance"}, or \code{type="cdf"}. See "Details".
#' @param newdata A new data frame containing the same variables declared in \code{Formula}. If missing, the dataset provided for estimation is used.
#' @param threshold Only needed if \code{type="cdf"} is selected. See "Details".
#' @param ... Further arguments passed to or from other methods.
#' @details The method \code{predict.LSBP_Gibbs} produces a sample of predicted values, obtained by evaluating the conditional mean of the LSBP model or the predictive distribution, using the observations contained in the \code{newdata} data frame.
#'
#' If \code{type="mean"}, then a sample from the posterior of the mean of a LSBP model is returned. If \code{type="predictive"} is selected, then a sample from the predictive distribution is returned.  If \code{type="variance"}, then a sample from the posterior distribution of the LSBP variance is returned.  If \code{type="cdf"}, then a sample from the posterior distribution of the LSBP cumulative distribution function is returned, evaluated at \code{threshold}.
#' @export
#'
predict.LSBP_Gibbs <- function(object, type = "mean", newdata = NULL, threshold = NULL, ...) {
  if (!(type %in% c("mean", "predictive", "variance", "cdf"))) {
    stop("Please provide a valid predictive type")
  }

  if (is.null(newdata)) {
    X1 <- object$data$X1
    X2 <- object$data$X2
  } else {
    f <- object$call
    X1 <- model.matrix(f, data = newdata, rhs = 1)
    if (length(f)[2] == 1) {
      X2 <- X1
    } else {
      X2 <- model.matrix(f, data = newdata, rhs = 2)
    }
  }

  n <- NROW(X2)
  p_kernel <- NCOL(X1)
  p_mixing <- NCOL(X2)

  if (type == "mean") {
    pred_mean <- matrix(0, object$control$R, n)
    if (NCOL(X1) > 1) {
      # Multivariate case
      for (r in 1:object$control$R) {
        pred_mean[r, ] <- pred_mean_multi(
          X1, X2, matrix(object$param$beta_mixing[r, , ], object$H - 1, p_mixing),
          object$param$beta_kernel[r, , ]
        )
      }
    } else {
      # Univariate case
      for (r in 1:object$control$R) {
        pred_mean[r, ] <- pred_mean(X2, matrix(object$param$beta_mixing[r, , ], object$H - 1, p_mixing), object$param$beta_kernel[r, ])
      }
    }
    return(pred_mean)
  }

  if (type == "predictive") {
    predictive <- matrix(0, object$control$R, n)
    if (NCOL(X1) > 1) {
      # Multivariate case
      for (r in 1:object$control$R) {
        predictive[r, ] <- predictive_multi(
          X1, X2, matrix(object$param$beta_mixing[r, , ], object$H - 1, p_mixing),
          object$param$beta_kernel[r, , ], object$param$tau[r, ]
        )
      }
    } else {
      # Univariate case
      for (r in 1:object$control$R) {
        predictive[r, ] <- predictive(X2, matrix(object$param$beta_mixing[r, , ], object$H - 1, p_mixing), object$param$beta_kernel[r, ], object$param$tau[r, ])
      }
    }
    return(predictive)
  }

  if (type == "variance") {
    variance <- matrix(0, object$control$R, n)
    if (NCOL(X1) > 1) {
      # Multivariate case
      for (r in 1:object$control$R) {
        variance[r, ] <- pred_var_multi(X1, X2, matrix(object$param$beta_mixing[r, , ], object$H - 1, p_mixing), object$param$beta_kernel[r, , ], object$param$tau[r, ])
      }
    } else {
      # Univariate case
      for (r in 1:object$control$R) {
        variance[r, ] <- pred_var(X2, matrix(object$param$beta_mixing[r, , ], object$H - 1, p_mixing), object$param$beta_kernel[r, ], object$param$tau[r, ])
      }
    }
    return(variance)
  }

  if (type == "cdf") {
    cdf <- matrix(0, object$control$R, n)
    if (NCOL(X1) > 1) {
      # Multivariate case
      for (r in 1:object$control$R) {
        cdf[r, ] <- pred_cdf_multi(X1, X2, matrix(object$param$beta_mixing[r, , ], object$H - 1, p_mixing), object$param$beta_kernel[r, , ], object$param$tau[r, ], threshold)
      }
    } else {
      # Univariate case
      for (r in 1:object$control$R) {
        cdf[r, ] <- pred_cdf(X2, matrix(object$param$beta_mixing[r, , ], object$H - 1, p_mixing), object$param$beta_kernel[r, ], object$param$tau[r, ], threshold)
      }
    }
    return(cdf)
  }
}

#' Predict method for the LSBP
#'
#'
#' Predict method for a LSBP model estimated using the \code{\link[LSBP]{LSBP_VB}} function.
#'
#' @param object An object of class \code{\link[LSBP]{LSBP_VB}}.
#' @param type String indicating the type of prediction: \code{type="mean"},\code{type="predictive"},  \code{type="variance"} or \code{type="cdf"}. See "Details".
#' @param R An integer indicating the number of replications for the returned sample.
#' @param newdata A new data frame containing the same variables declared in \code{Formula}. If missing, the dataset provided for estimation is used.
#' @param threshold Only needed if \code{type="cdf"} is selected. See "Details".
#' @param ... Further arguments passed to or from other methods.
#' @details The method \code{predict.LSBP_VB} produces a sample of predicted values, obtained by evaluating the conditional mean of the LSBP model or the predictive distribution, using the observations contained in the \code{newdata} data frame.
#'
#'  If \code{type="mean"}, then a sample from the (variational) posterior of the mean of a LSBP model is returned. If \code{type="predictive"} is selected, then a sample from the (variational) predictive distribution is returned.  If \code{type="variance"}, then a sample from the (variational) posterior distribution of the LSBP variance is returned.  If \code{type="cdf"}, then a sample from the (variational) posterior distribution of the LSBP cumulative distribution function is returned, evaluated at \code{threshold}.
#' @export
#'
predict.LSBP_VB <- function(object, type = "mean", R = 5000, newdata = NULL, threshold = NULL, ...) {
  if (!(type %in% c("mean", "predictive", "variance", "cdf"))) {
    stop("Please provide a valid predictive type")
  }
  if (is.null(newdata)) {
    X1 <- object$data$X1
    X2 <- object$data$X2
  } else {
    f <- object$call
    X1 <- model.matrix(f, data = newdata, rhs = 1)
    if (length(f)[2] == 1) {
      X2 <- X1
    } else {
      X2 <- model.matrix(f, data = newdata, rhs = 2)
    }
  }

  n <- NROW(X2)
  p_kernel <- NCOL(X1)
  p_mixing <- NCOL(X2)
  beta_mixing <- array(0, c(R, object$H - 1, p_mixing))
  tau <- matrix(0, R, object$H)

  # Generating the parameters

  if (NCOL(X1) > 1) {
    beta_kernel <- array(0, c(R, object$H, p_kernel))
    for (h in 1:object$H) {
      if (h < object$H) {
        eig <- eigen(object$param$Sigma_mixing[h, , ], symmetric = TRUE)
        A1 <- t(eig$vectors) * sqrt(eig$values)
        beta_mixing[, h, ] <- t(object$param$mu_mixing[h, ] + t(matrix(rnorm(R * p_mixing), R, p_mixing) %*%
          A1))
      }

      eig <- eigen(object$param$Sigma_kernel[h, , ], symmetric = TRUE)
      A1 <- t(eig$vectors) * sqrt(eig$values)
      beta_kernel[, h, ] <- t(object$param$mu_kernel[h, ] + t(matrix(rnorm(R * p_kernel), R, p_kernel) %*%
        A1))
      tau[, h] <- rgamma(R, object$param$a_tilde[h], object$param$b_tilde[h])
    }
  } else {
    beta_kernel <- matrix(0, R, object$H)
    for (h in 1:object$H) {
      if (h < object$H) {
        eig <- eigen(object$param$Sigma_mixing[h, , ], symmetric = TRUE)
        A1 <- t(eig$vectors) * sqrt(eig$values)
        beta_mixing[, h, ] <- t(object$param$mu_mixing[h, ] + t(matrix(rnorm(R * p_mixing), R, p_mixing) %*%
          A1))
      }
      beta_kernel[, h] <- rnorm(R, object$param$mu_kernel[h], sqrt(object$param$Sigma_kernel[h]))
      tau[, h] <- rgamma(R, object$param$a_tilde[h], object$param$b_tilde[h])
    }
  }

  if (type == "mean") {
    pred <- matrix(0, R, n)
    if (NCOL(X1) > 1) {
      for (r in 1:R) {
        pred[r, ] <- pred_mean_multi(X1, X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), beta_kernel[r, , ])
      }
    } else {
      for (r in 1:R) {
        pred[r, ] <- pred_mean(X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), beta_kernel[r, ])
      }
    }
  }

  if (type == "predictive") {
    pred <- matrix(0, R, n)
    if (NCOL(X1) > 1) {
      for (r in 1:R) {
        pred[r, ] <- predictive_multi(
          X1, X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing),
          beta_kernel[r, , ], tau[r, ]
        )
      }
    } else {
      for (r in 1:R) {
        pred[r, ] <- predictive(X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), beta_kernel[r, ], tau[r, ])
      }
    }
  }

  if (type == "variance") {
    pred <- matrix(0, R, n)
    if (NCOL(X1) > 1) {
      for (r in 1:R) {
        pred[r, ] <- pred_var_multi(
          X1, X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing),
          beta_kernel[r, , ], tau[r, ]
        )
      }
    } else {
      for (r in 1:R) {
        pred[r, ] <- pred_var(X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), beta_kernel[r, ], tau[r, ])
      }
    }
  }

  if (type == "cdf") {
    pred <- matrix(0, R, n)
    if (NCOL(X1) > 1) {
      for (r in 1:R) {
        pred[r, ] <- pred_cdf_multi(
          X1, X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing),
          beta_kernel[r, , ], tau[r, ], threshold
        )
      }
    } else {
      for (r in 1:R) {
        pred[r, ] <- pred_cdf(X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), beta_kernel[r, ], tau[r, ], threshold)
      }
    }
  }

  return(pred)
}

#' Conditional density of a LSBP model
#'
#'
#' Evaluate the conditional density \eqn{f_x(y)} of a LSBP model, given the parameters and the covariates.
#'
#' @param y The value at which the conditional density must be evaluated.
#' @param X1 A \code{n x p_kernel} design matrix for the kernel.
#' @param X2 A \code{n x p_mixing} design matrix for the stick-breaking weights.
#' @param beta_mixing A \code{H-1 x p_mixing} dimensional matrix of coefficients for the linear predictor of the stick-breaking weights.
#' @param beta_kernel A \code{H x p_kernel} dimensional matrix of coefficients for the linear predictor of the kernel.
#' @param tau A \code{H} dimensional vector of coefficients for the kernel precision.
#' @details The function \code{LSBP_density} evaluates the conditional density \eqn{f_x(y)}. The number of mixture components \code{H} is inferred from the dimensions of \code{beta_mixing} and \code{beta_kernel}.
#' @export
#'
LSBP_density <- function(y, X1, X2, beta_mixing, beta_kernel, tau) {
  LSBP_density_C(y, X1, X2, beta_mixing, beta_kernel, tau)
}
