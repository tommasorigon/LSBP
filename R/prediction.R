#' Hello
#' 
#' 
#' @export
predict.DLSBP_ECM <- function(x, newdata=NULL){
  if(is.null(newdata)) {
    X1 <- x$data$X1; X2 <- x$data$X2
  } else {
    f <- x$call
    X1 <- model.matrix(f, data=newdata, rhs = 1)
    if(length(f)[2]==1) {X2 <- X1} else {X2 <- model.matrix(f, data=newdata, rhs = 2)}
  }
  
  if(NCOL(X1)>1){
    pred    <- pred_mean_multi(X1, X2, x$param$beta_mixing, x$param$beta_kernel, x$param$tau)
    
  } else {
    # Univariate case
    pred    <- pred_mean(X2, x$param$beta_mixing, x$param$beta_kernel, x$param$tau)
  }
  return(as.numeric(pred))
}

#' @export
predict.DLSBP_Gibbs <- function(x, type="mean", newdata=NULL) {
  
  if(is.null(newdata)) {
    X1 <- x$data$X1; X2 <- x$data$X2
  } else {
    f <- x$call
    X1 <- model.matrix(f, data=newdata, rhs = 1)
    if(length(f)[2]==1) {X2 <- X1} else {X2 <- model.matrix(f, data=newdata, rhs = 2)}
  }
  
  n         <- NROW(X2)
  
  if(type=="mean"){
    pred_mean <- matrix(0,x$control$R,n)
    if(NCOL(X1)>1){
      #Multivariate case
      for(r in 1:x$control$R){pred_mean[r,] <- pred_mean_multi(X1,X2,x$param$beta_mixing[r,,],
                                                               x$param$beta_kernel[r,,], x$param$tau[r,])}
    } else {
      # Univariate case
      for(r in 1:x$control$R){
        pred_mean[r,] <- pred_mean(X2,x$param$beta_mixing[r,,], x$param$beta_kernel[r,], x$param$tau[r,])
      }
    }
    return(pred_mean)
    }
  
  if(type=="predictive"){
    predictive <- matrix(0,x$control$R,n)
    if(NCOL(X1)>1){
      #Multivariate case
      for(r in 1:x$control$R){predictive[r,] <- predictive_multi(X1,X2,x$param$beta_mixing[r,,],
                                                               x$param$beta_kernel[r,,], x$param$tau[r,])}
    } else {
      # Univariate case
      for(r in 1:x$control$R){
        predictive[r,] <- predictive(X2,x$param$beta_mixing[r,,], x$param$beta_kernel[r,], x$param$tau[r,])
      }
    }
    return(predictive)
  }
}

#' @export
predict.DLSBP_VB <- function(x, type="mean", R=5000, newdata=NULL){
  
  if(is.null(newdata)) {
    X1 <- x$data$X1; X2 <- x$data$X2
  } else {
    f <- x$call
    X1 <- model.matrix(f, data=newdata, rhs = 1)
    if(length(f)[2]==1) {X2 <- X1} else {X2 <- model.matrix(f, data=newdata, rhs = 2)}
  }
  
  n         <- NROW(X2); p_kernel <- NCOL(X1); p_mixing <- NCOL(X2)
  beta_mixing      <- array(0,c(R,x$H-1,p_mixing))
  tau       <- matrix(0,R, x$H)
  
  # Generating te parameters
  
  if(NCOL(X1)>1){
    beta_kernel      <- array(0,c(R,x$H,p_kernel))
    for(h in 1:x$H){
      if(h < x$H) {
        eig       <- eigen(x$param$Sigma_beta_mixing[h,,],symmetric = TRUE)
        A1        <- t(eig$vectors)*sqrt(eig$values)
        beta_mixing[, h, ] <- t(x$param$mu_beta_mixing[h,] + t(matrix(rnorm(R*p_mixing),R,p_mixing)%*%A1))
      }
      
      eig         <- eigen(x$param$Sigma_beta_kernel[h,,],symmetric = TRUE)
      A1          <- t(eig$vectors)*sqrt(eig$values)
      beta_kernel[, h, ]<- t(x$param$mu_beta_kernel[h,] + t(matrix(rnorm(R*p_kernel),R,p_kernel)%*%A1))
      
      tau[r,]   <- rgamma(R,x$param$a_tilde[h],x$param$b_tilde[h])
    }
    } else {
    beta_kernel       <- matrix(0,R, x$H)
    for(h in 1:x$H){
      if(h < x$H) {
        eig       <- eigen(x$param$Sigma_beta_mixing[h,,],symmetric = TRUE)
        A1        <- t(eig$vectors)*sqrt(eig$values)
        beta_mixing[, h, ] <- t(x$param$mu_beta_mixing[h,] + t(matrix(rnorm(R*p_mixing),R,p_mixing)%*%A1))
      }
      beta_kernel[,h]  <- rnorm(R,x$param$mu_tilde[h],1/sqrt(x$param$tau_tilde[h]))
      tau[,h] <- rgamma(R,x$param$a_tilde[h],x$param$b_tilde[h])
    }
    
  }
  
  if(type=="mean"){
    pred <- matrix(0,R,n)
    if(NCOL(X1)>1){
      for(r in 1:R){pred[r,] <- pred_mean_multi(X1,X2,beta_mixing[r,,], beta_kernel[r,,], tau[r,])}} else {
      for(r in 1:R){pred[r,] <- pred_mean(X2, beta_mixing[r,,],beta_kernel[r,], tau[r,])}
      }
  }
  
  if(type=="predictive"){
    pred <- matrix(0,R,n)
    if(NCOL(X1)>1){
      for(r in 1:R){pred[r,] <- predictive_multi(X1,X2,beta_mixing[r,,], beta_kernel[r,,], tau[r,])}} else {
      for(r in 1:R){pred[r,] <- predictive(X2, beta_mixing[r,,],beta_kernel[r,], tau[r,])}
      }
  }
  return(pred)
}
