# DLSBP: Dependent Logit Stick-Breaking Process

R functions and scripts to implement the Dependent Logit Stick-Breaking process model of

* Rigon, T. and Durante, D., (2017), [Logit stick-breaking priors for Bayesian density regression](https://arxiv.org/abs/1701.02969), ArXiv.

The `DLSBP` package can be installed as follow

```R
# If the devtools R package is not already installed
# install.packages("devtools")

devtools::install_github("tommasorigon/DLSBP")
```

Alternatively, install the R package from the source file `DLSBP_VERSION.tar.gz` with `R CMD INSTALL DLSBP_VERSION.tar.gz`. The analysis of the Old Faithful Geyser can be reproduced running the R script `old_faithful.R`. For instance, the Expectation Maximization (EM) algorithm is run through

```R
library(DLSBP)
data(geyser)

x <- geyser$waiting
y <- geyser$duration
X <- cbind(1,x*I(x<=68) + 68*I(x > 68),((x-68)*I(x>68))) # Predictors
p <- NCOL(X)

prior   <- list(b = rep(0,p),        
                B = diag(100,p), 
                mu_mu =0,            
                tau_mu = 0.001,      
                a_tau = 2,           
                b_tau = 0.001)

fit_EM   <- DLSBP_EM(y = y, X = X, H = 3, prior = prior)
```
