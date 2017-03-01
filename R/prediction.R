#' Predict method for the DLSBP
#' 
#' 
#' Predict method for the DLSBP estimated using the ECM algorithm.
#' 
#' @param object An object of class \verb{DLSBP_EM}.
#' @param newdata A new data frame containing the same variables declared in \verb{Formula}. If missing, the dataset provided for estimation is used.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return The method \verb{predict.DLSBP_ECM} produces predicted values, obtained by evaluating the conditional mean of the DLSBP model, after plugging-in the MAP, and using the observations contained in the \verb{newdata} data frame.
#' @export
predict.DLSBP_ECM <- function(object, newdata=NULL,  ...) {
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
   
   if (NCOL(X1) > 1) {
      pred <- pred_mean_multi(X1, X2, object$param$beta_mixing, object$param$beta_kernel, object$param$tau)
      
   } else {
      # Univariate case
      pred <- pred_mean(X2, object$param$beta_mixing, object$param$beta_kernel, object$param$tau)
   }
   return(as.numeric(pred))
}

#' Predict method for the DLSBP
#' 
#' 
#' Predict method for the DLSBP estimated using the Gibbs sampling algorithm.
#' 
#' @param object An object of class \verb{DLSBP_Gibbs}.
#' @param type String indicating the type of prediction, either \verb{type="mean"} or \verb{type="predictive"}. See details.
#' @param newdata A new data frame containing the same variables declared in \verb{Formula}. If missing, the dataset provided for estimation is used.
#' @param ... Further arguments passed to or from other methods.
#' @details The method \verb{predict.DLSBP_Gibbs} produces a sample of predicted values, obtained by evaluating the conditional mean of the DLSBP model or the predictive distribution, using the observations contained in the \verb{newdata} data frame. 
#' 
#' If \verb{type="mean"} a sample from the posterior distribution of the DLSBP mean is returned. If \verb{type="predictive"} is selected, then a sample from the predictive distribution is returned.
#' @export
predict.DLSBP_Gibbs <- function(object, type = "mean", newdata = NULL,...) {
  if (!(type %in% c("mean", "predictive"))) 
    stop("Please provide a valid predictive type")
  
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
   
   n <- NROW(X2); p_kernel <- NCOL(X1); p_mixing <- NCOL(X2)
   
   if (type == "mean") {
      pred_mean <- matrix(0, object$control$R, n)
      if (NCOL(X1) > 1) {
         # Multivariate case
         for (r in 1:object$control$R) {
            pred_mean[r, ] <- pred_mean_multi(X1, X2, matrix(object$param$beta_mixing[r, , ],object$H-1, p_mixing), 
              object$param$beta_kernel[r, , ], object$param$tau[r, ])
         }
      } else {
         # Univariate case
         for (r in 1:object$control$R) {
            pred_mean[r, ] <- pred_mean(X2, matrix(object$param$beta_mixing[r, , ], object$H-1, p_mixing), object$param$beta_kernel[r, 
              ], object$param$tau[r, ])
         }
      }
      return(pred_mean)
   }
   
   if (type == "predictive") {
      predictive <- matrix(0, object$control$R, n)
      if (NCOL(X1) > 1) {
         # Multivariate case
         for (r in 1:object$control$R) {
            predictive[r, ] <- predictive_multi(X1, X2, matrix(object$param$beta_mixing[r, , ],object$H-1,  p_mixing), 
              object$param$beta_kernel[r, , ], object$param$tau[r, ])
         }
      } else {
         # Univariate case
         for (r in 1:object$control$R) {
            predictive[r, ] <- predictive(X2, matrix(object$param$beta_mixing[r, , ],object$H-1,  p_mixing), object$param$beta_kernel[r, 
              ], object$param$tau[r, ])
         }
      }
      return(predictive)
   }
}

#' Predict method for the DLSBP
#' 
#' 
#' Predict method for the DLSBP estimated using the Variational Bayes algorithm.
#' 
#' @param object An object of class \verb{DLSBP_VB}.
#' @param type String indicating the type of prediction, either \verb{type="mean"} or \verb{type="predictive"}. See details.
#' @param R An integer indicating the number of replications for the returned sample.
#' @param newdata A new data frame containing the same variables declared in \verb{Formula}. If missing, the dataset provided for estimation is used.
#' @param ... Further arguments passed to or from other methods.
#' @details The method \verb{predict.DLSBP_VB} produces a sample of predicted values, obtained by evaluating the conditional mean of the DLSBP model or the predictive distribution, using the observations contained in the \verb{newdata} data frame. 
#' 
#' If \verb{type="mean"} a sample from the posterior distribution of the DLSBP mean is returned. If \verb{type="predictive"} is selected, then a sample from the predictive distribution is returned.
#' @export
predict.DLSBP_VB <- function(object, type = "mean", R = 5000, newdata = NULL, ...) {
  if (!(type %in% c("mean", "predictive"))) 
    stop("Please provide a valid predictive type")
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
   
   # Generating te parameters
   
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
            pred[r, ] <- pred_mean_multi(X1, X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), beta_kernel[r, 
              , ], tau[r, ])
         }
      } else {
         for (r in 1:R) {
            pred[r, ] <- pred_mean(X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), beta_kernel[r, 
              ], tau[r, ])
         }
      }
   }
   
   if (type == "predictive") {
      pred <- matrix(0, R, n)
      if (NCOL(X1) > 1) {
         for (r in 1:R) {
            pred[r, ] <- predictive_multi(X1, X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), 
              beta_kernel[r, , ], tau[r, ])
         }
      } else {
         for (r in 1:R) {
            pred[r, ] <- predictive(X2, matrix(beta_mixing[r, , ], object$H - 1, p_mixing), beta_kernel[r, 
              ], tau[r, ])
         }
      }
   }
   return(pred)
}
