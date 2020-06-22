# LSBP: Bayesian density regression via logit stick-breaking

The repository includes the `R` functions to implement the **logit stick-breaking process model** of the paper:

* Rigon, T. and Durante, D. (2020). [Tractable Bayesian Density Regression via Logit Stick-Breaking Priors](https://doi.org/10.1016/j.jspi.2020.05.009), Journal of Statistical Planning and Inference.

The `LSBP` package can be installed by running the following commands:

```R
# If the devtools R package is not already installed
# install.packages("devtools")

devtools::install_github("tommasorigon/LSBP")
```

Alternatively, install the `R` package from the source file `LSBP_VERSION.tar.gz` with `R CMD INSTALL LSBP_VERSION.tar.gz`. 

The documentation of the `LSBP` package is available in [this repository](https://github.com/tommasorigon/LSBP/blob/master/LSBP_0.0.4.pdf). The file [`tutorial.md`](https://github.com/tommasorigon/LSBP/blob/master/Tutorial/tutorial.md) contains detailed instruction for reproducing the epidemiological application of **Section 4** in the paper.
