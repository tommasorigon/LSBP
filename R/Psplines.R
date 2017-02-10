linvgamma <- function (x, shape, scale) {
  return(- (shape + 1) * log(x) - (scale/x))
}

#' @export
Pspline_fit <- function(x, y, inner_knots = min(length(y)/4,40), degree=3, dif = 1, prior=NULL, maxiter = 50, verbose = FALSE) {
  
  # If prior hyperparameters are not defined, then use the following vague prior
  if (is.null(prior)) {
    prior <- list(
      a_sigma2 = 0.001, 
      b_sigma2 = 0.001,
      a_tau2=0.001,
      b_tau2=0.001)
  }
  
  # Convergence tolerance
  tol=1e-5
  
  # Knots placement
  n     <- length(x)
  xl    <- min(x); xr <- max(x); dx <- (xr - xl) / (inner_knots-1)
  knots <- seq(xl - degree * dx, xr + degree * dx, by = dx)
  
  # Prior settings
  a_sigma2    <- prior$a_sigma2
  b_sigma2    <- prior$b_sigma2
  a_tau2      <- prior$a_tau2
  b_tau2      <- prior$b_tau2
  
  # Fixed quantities
  B     <- spline.des(knots, x, degree + 1, 0 * x, outer.ok=TRUE, sparse=TRUE)$design
  BtB   <- crossprod(B)
  rankB <- NCOL(B) - dif
  
  if(dif==0) {D <- Diagonal(NCOL(B))}
  else{D   <- Matrix(diff(diag(NCOL(B)),dif=dif),sparse=TRUE)}
  
  DtD <- crossprod(D)
  Bty <- crossprod(B,y)
  
  # Initialization
  sigma2 <- var(y)/4
  tau2   <- var(y)/4
  logpost <- -Inf
  
  # Start CM algorithm
  for(r in 1:maxiter){
    
    L          <- Cholesky(BtB/sigma2 + DtD/tau2)
    beta       <- solve(L,Bty/sigma2)

    # Beta quantities
    mahalanob     <- as.numeric(crossprod(D%*%beta)) 
    pred          <- as.numeric(B%*%beta) 
    
    # Smoothing component
    a_tau2_tilde  <- a_tau2 + rankB/2
    b_tau2_tilde  <- b_tau2 + mahalanob/2
    tau2          <- b_tau2_tilde/ (a_tau2_tilde+1)

    # Variance
    a_sigma2_tilde <- a_sigma2 + n/2
    b_sigma2_tilde <- b_sigma2 + sum((y-pred)^2)/2
    sigma2 <-  b_sigma2_tilde / (a_sigma2_tilde+1)
    
    # Convergence checks
    loglik         <-  sum(dnorm(y,pred,sqrt(sigma2),log=TRUE))
    logpriors      <-  -0.5*mahalanob/tau2 - 0.5*rankB*log(tau2) + 
      linvgamma(sigma2,a_sigma2, b_sigma2) + 
      linvgamma(tau2,  a_tau2,   b_tau2)
    
    logpost_new    <- loglik + logpriors
    
    # Break the loop at convergence
    if((logpost_new - logpost)<tol) {
      if(verbose) cat(paste("Convergence reached after",r,"iterations."))
      break
    }
    logpost <- logpost_new
    
    # Display status
    if (verbose) {
      cat(paste("log-posterior:",round(logpost,4),", iteration:", r, "\n",sep=""))
    }
  }
  out <- list(param=list(beta=beta, sigma2=sigma2, tau2=tau2), data=list(x=x,y=y),fit=list(knots=knots,degree=degree,dif=dif,iter=r,logpost=logpost))
  attr(out,"class") <- "PSpline"
  return(out)
}

#' @export
print.PSpline <- function(x) cat(paste("Smoothing parameter: tau2=",round(x$param$tau2,6),".\nVariance parameter: sigma2=",round(x$param$sigma2,6),".\nConvergence reached after ",x$fit$iter," iterations.", sep=""))


#' @export
predict.PSpline <- function(x,newx=NULL) {
  if(!is.null(newx)) {x_grid <- newx} 
  else{x_grid <- x$data$x}
  
  B     <- spline.des(x$fit$knots, x_grid, x$fit$degree + 1, 0 * x_grid, outer.ok=TRUE, sparse=TRUE)$design
  pred  <- as.numeric(B %*% x$param$beta)
  return(pred)
}

#' @export
plot.PSpline <- function(x,ngrid=500) {
  x_grid <- seq(from=min(x$data$x),to=max(x$data$x),length=ngrid)
  fit <- predict(x,newx=x_grid)
  plot(x$data$x,x$data$y,xlab="x",ylab="y",cex=.6);
  lines(x_grid,fit,col="red",lwd=2)
}







