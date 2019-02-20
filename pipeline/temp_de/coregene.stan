
data {
  int<lower = 1> E; // number of experiments
  
  real beta_obs[E];
  real<lower = 0> beta_se[E];
}

parameters {
  real beta;
  real<lower = 0> sigma;
}

model {
  
  beta ~ normal(0, 1);
  sigma ~ gamma(2, 1);
  
  for(e in 1:E) {
    beta_obs[e] ~ normal(beta, sigma + beta_se[e]);
  }
  
}

