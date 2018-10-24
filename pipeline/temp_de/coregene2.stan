
data {
  int<lower = 1> E; // number of experiments
  int<lower = 1> G; // number of genes
  
  real beta_obs[G,E];
  real<lower = 0> beta_se[G,E];
  // int<lower = 0, upper = 1> is_missing[G,E];
}

parameters {
  real phi;
  real<lower = 0> lambda;
  
  real beta[G];
  real<lower = 0> sigma[G];
  

}

model {
  
  phi ~ normal(0, 1);
  
  lambda ~ gamma(2,1);
  sigma ~ gamma(2, 1);
  
  for(g in 1:G) {
    beta ~ normal(phi, lambda);
    for(e in 1:E) {
      //if(is_missing[g,e] == 0) {
        beta_obs[g,e] ~ normal(beta[g], sigma[g] + beta_se[g,e]);
      //}
    }
  }
  
}

