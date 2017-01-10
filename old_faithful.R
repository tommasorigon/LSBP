library(DLSBP)
rm(list=ls())
data(geyser)

x <- geyser$waiting
y <- geyser$duration
X <- cbind(1,x*I(x<=68) + 68*I(x > 68),((x-68)*I(x>68))) 
p <- NCOL(X)
R <- 30000; burn.in <- 5000

prior   <- list(b = rep(0,p),        
                B = diag(100,p),     
                mu_mu =0,            
                tau_mu = 0.001,      
                a_tau = 2,           
                b_tau = 0.001)

fit_EM      <- DLSBP_EM(y, X, H=3, prior=prior,maxiter=10)
fit_VB      <- DLSBP_VB(y=y, X=X, H=3, prior=prior, seed=99)

set.seed(99)
fit_Gibbs   <- DLSBP_Gibbs(y=y, X=X, H=3, prior=prior, R=R,burn.in = burn.in,EM_initialization=TRUE)

fit_EM$beta
