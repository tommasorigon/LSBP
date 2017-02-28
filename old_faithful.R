rm(list=ls())
library(DLSBP)
data(geyser,package = "MASS")

x <- geyser$waiting
y <- geyser$duration
geyser <- data.frame(y=y, x=x, x1=x*I(x<=68) + 68*I(x > 68),x2=((x-68)*I(x>68))) 

R <- 30000; burn_in <- 5000

fit_em    <- DLSBP_EM(y ~ 1 | x1 + x2, data=geyser,  H=3)
set.seed(10)
fit_vb    <- DLSBP_VB(y ~ 1 | x1 + x2, data=geyser,  H=3)
fit_gibbs <- DLSBP_Gibbs(y ~ 1 | x1 + x2, data=geyser, H=3)

plot(x,y,col=fit_em$cluster,cex=.8)

fit_VB      <- DLSBP_VB(y ~ 1 | x1 + x2, data=geyser,  H=3)
fit_EM      <- DLSBP_EM(f= y ~ x | x1 + x2, data=geyser,  H=3, control=control_EM(tol=1e-5))
fit_VB      <- DLSBP_VB(f= y ~ x | x1 + x2, data=geyser,  H=3, control=control_EM(tol=1e-5))

set.seed(99)
fit_Gibbs   <- DLSBP_Gibbs(y ~ 1 | x1 + x2, data=geyser, H=3,control = control_Gibbs(R=R,burn_in=burn_in))

print(fit_Gibbs,plot=TRUE)

plot(x,y,cex=0.8)
points(x,predict(fit_EM),col=fit_EM$cluster,pch=16)
points(x,colMeans(fit_Gibbs$pred),col="red",pch=16)
points(x,apply(fit_Gibbs$pred,2, function(x) quantile(x,0.99)),col="blue",pch=16)
points(x,apply(fit_Gibbs$pred,2, function(x) quantile(x,0.1)),col="blue",pch=16)
