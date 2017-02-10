// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
List G_update(arma::vec y, arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau){

  // Initialization
  int H = mu.n_elem;
  int n = X.n_rows;

  arma::vec G(n);                      // Output: cluster index
  arma::vec G_mixt(n);                 // Output: mixture index
  IntegerVector clusters = Range(1,H); // Range of possibile clusters
  arma::vec lprob(H);                  // Conditional log-probabilities
  arma::vec prob(H);                   // Conditional probabilities
  arma::vec pi(H);                     // Mixture weights
  arma::vec nu(H);                     // Stick-breaking weights
  arma::vec cum_nu(H);                 // Cumulate product of stick-breaking weights
  
  // Other initializations
  double mean;
  double sd;
  
  // Cycle the observations
  for(int i=0; i < n; i++) { 
    
    // First mixture component
    nu[0]     = 1/(1+ exp(- (dot(X.row(i),beta.row(0)))));
    pi[0]     = nu[0];
    cum_nu[0] = 1-nu[0];
    
    mean = mu[0];
    sd   = sqrt(1/tau[0]);
    lprob[0]    = log(pi[0]) + R::dnorm(y(i),mean,sd,TRUE);

    // Start from h = 1 (second mixture component), compute the probabilities until H - 2
    for(int h = 1; h < H - 1; h++) {

      nu[h]    = 1/(1 + exp(- (dot(X.row(i),beta.row(h)))));
      pi[h]    = nu[h] * cum_nu[h-1];
      cum_nu[h] = (1-nu[h])*cum_nu[h-1];
      
      mean = mu[h];
      sd   = sqrt(1/tau[h]);
      lprob[h] = log(pi[h]) + R::dnorm(y(i),mean,sd,TRUE);
    }
    
    // Last mixture component
    nu[H-1] = 1;                      // Setted equal to 1.
    pi[H-1] = nu[H-1] * cum_nu[H-2];  
    
    mean = mu[H-1];
    sd   = sqrt(1/tau[H-1]);
    lprob[H-1]    = log(pi[H-1]) + R::dnorm(y(i),mean,sd,TRUE);
    
    // Steps for numerical stabilization - normalization
    lprob = lprob - max(lprob);
    prob  = exp(lprob);
    prob  = prob/sum(prob);
    
    // Sampling the cluster and mixture index
    G[i]      = RcppArmadillo::sample(clusters, 1, TRUE, prob)[0];
    G_mixt[i] = RcppArmadillo::sample(clusters, 1, TRUE, pi)[0];
  }
  List out = List::create(
    Named("G") = G,
    Named("G_mixt") = G_mixt
  );
  return out;
}

// [[Rcpp::export]]
List Expectation_step(arma::vec y, arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau){
  
  // Initialization
  int H = mu.n_elem;
  int n = X.n_rows;
  
  arma::mat z(n,H);                    // Output probabilities
  IntegerVector clusters = Range(1,H); // Range of possibile clusters
  arma::vec lprob(H);                  // Conditional log-probabilities
  arma::vec prob(H);                   // Conditional probabilities
  arma::vec pi(H);                     // Mixture weights
  arma::vec nu(H);                     // Stick-breaking weights
  arma::vec cum_nu(H);                 // Cumulate product of stick-breaking weights
  
  double mean;
  double sd;
  double loglik = 0;
  
  // Cycle the observations
  for(int i=0; i < n; i++) { 
    
    // First mixture component
    nu[0]     = 1/(1+ exp(- (dot(X.row(i),beta.row(0)))));
    pi[0]     = nu[0];
    cum_nu[0] = 1-nu[0];
    
    mean = mu[0];
    sd   = sqrt(1/tau[0]);
    lprob[0]    = log(pi[0]) + R::dnorm(y(i),mean,sd,TRUE);
    
    // Start from h = 1 and arrives to H - 2, that is, exclude the first and the last component
    for(int h = 1; h < H - 1; h++) {
      
      nu[h]    = 1/(1 + exp(- (dot(X.row(i),beta.row(h)))));
      pi[h]    = nu[h] * cum_nu[h-1];
      cum_nu[h] = (1-nu[h])*cum_nu[h-1];
      
      mean = mu[h];
      sd   = sqrt(1/tau[h]);
      lprob[h] = log(pi[h]) + R::dnorm(y(i),mean,sd,TRUE);
    }
    
    // Last component
    nu[H-1] = 1; // Set it equal to 1.
    pi[H-1] = nu[H-1] * cum_nu[H-2];
    mean = mu[H-1];
    sd   = sqrt(1/tau[H-1]);
    lprob[H-1]    = log(pi[H-1]) + R::dnorm(y(i),mean,sd,TRUE);
    
    loglik = loglik + log(sum(exp(lprob))); //Here it is not necessarily numerically stable.
    
    // Steps for numerical stabilization - normalization
    lprob = lprob - max(lprob);
    prob  = exp(lprob);
    prob  = prob/sum(prob);
    // Sampling
    z.row(i) = prob.t();
  }
  List out = List::create(
    Named("z") = z,
    Named("loglik") = loglik
  );
  return out;
}

// [[Rcpp::export]]
List Variational_step(arma::mat rho, arma::vec y, arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau_tilde, arma::vec tau, arma::vec ltau){
  
  // Initialization
  int H = mu.n_elem;
  int n = X.n_rows;
  
  arma::mat zl(n,H); // Weights probabilities
  arma::mat z(n,H);
  arma::vec cum_nu(H);
  
  rho.insert_cols(H-1,arma::ones<arma::mat>(n,1)); // Add a fake column having all ones
  
  // Auxiliary quantities
  double eta = 0;
  arma::uvec indexes;

 for(int i=0; i < n; i++) { 
    
    // First loop, with h=0,1,2,...,H-2
    for(int h = 0; h < H - 1; h++) {
    // Second loop, i=0,1,2,..,n-1
    cum_nu(0) = 1;
  
     arma::uvec row;
     row = i;
     // Third loop, l = 0,1,...,H-1
     for(int l=0; l < H; l++) {
       
       // Positive contributions
       if(l==h){
         if(h == 0) {
           zl(i,l) = 1;
         }
         else{
           indexes    =  arma::regspace<arma::uvec>(0,  l - 1);  // 0,  1, ...,  l-1
           zl(i,l)    =  arma::prod(1 - arma::vectorise(rho(row,indexes)));
         }
       }
       // Negative contributions
       else if(l > h) {
         indexes = arma::regspace<arma::uvec>(0,  l - 1);  // 0,  1, ...,  l-1
         indexes = arma::find(indexes != h);
         if(h == 0 && l == 1){
           zl(i,l) = - rho(i,l);
         }
         else{
           zl(i,l) = - rho(i,l)*arma::prod(1 - arma::vectorise(rho(row,indexes)));
         }
       }
       else {
         zl(i,l) = 0;
       }
     }
     
     // Linear predictor
     eta      = dot(X.row(i),beta.row(h)) + sum(zl.row(i).t() % (ltau/2  - tau/2 % (pow(y(i),2) - 2*y(i)*mu + pow(mu,2) + 1/tau_tilde)));
     rho(i,h) = 1/(1+ exp(-eta));
     
     // Stick-breaking
     z(i,h)  = rho(i,h)*cum_nu(h);
     cum_nu(h+1) = cum_nu(h)*(1 - rho(i,h));
     }
    z(i,H-1) = cum_nu(H-1);
    }
 rho.shed_col(H-1);

 List out = List::create(
   Named("rho") = rho,
   Named("z")   = z
   );
 return out;
}

// [[Rcpp::export]]
arma::vec prediction(arma::mat X, arma::mat beta, arma::vec mu, arma::vec tau){

  // Initialization
  int H = mu.n_elem;
  int n = X.n_rows;
  
  arma::vec pred = arma::zeros<arma::vec>(n);   // Vector with the predictions
  arma::vec pi(H);                              // Mixture weights
  arma::vec nu(H);                              // Stick-breaking weights
  arma::vec cum_nu(H);                          // Cumulate product of stick-breaking weights
  
  // Cycle the observations
  for(int i=0; i < n; i++) { 
    // First mixture component
    nu[0]     = 1/(1+ exp(- (dot(X.row(i),beta.row(0)))));
    pi[0]     = nu[0];
    cum_nu[0] = 1-nu[0];
    pred[i]   = pi[0]*mu[0];
    
    // Start from h = 1 and arrives to H - 2, that is, exclude the first and the last component
    for(int h = 1; h < H - 1; h++) {
      nu[h]      = 1/(1 + exp(- (dot(X.row(i),beta.row(h)))));
      pi[h]     = nu[h] * cum_nu[h-1];
      cum_nu[h]  = (1-nu[h])*cum_nu[h-1];
      pred[i]   = pred[i] + pi[h]*mu[h];
    }
    
    // Last component
    nu[H-1]   = 1; // Set it equal to 1.
    pi[H-1] = nu[H-1] * cum_nu[H-2];
    pred[i]   = pred[i] + pi[H-1]*mu[H-1];
  }

return pred;
}
