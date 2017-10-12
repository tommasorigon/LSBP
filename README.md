# LSBP: Bayesian density regression via **L**ogit **S**tick-**B**reaking **P**rocess

The repository includes the R functions to implement the Logit Stick-Breaking process model of the paper

* Rigon, T. and Durante, D., (2017), [Tractable Bayesian Density Regression via Logit Stick-Breaking Priors](https://arxiv.org/abs/1701.02969), ArXiv.

The `LSBP` package can be installed by running the following commands

```R
# If the devtools R package is not already installed
# install.packages("devtools")

devtools::install_github("tommasorigon/LSBP")
```

Alternatively, install the R package from the source file `LSBP_VERSION.tar.gz` with `R CMD INSTALL LSBP_VERSION.tar.gz`. 

The file [`tutorial.md`](https://github.com/tommasorigon/LSBP/blob/master/Tutorial/tutorial.md) contains detailed instruction for reproducing the epidemiological application of Section 4 in the paper.