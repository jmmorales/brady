
functions { 
  
  // copied from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_zero_inflated_binomial.stan
  /* zero-inflated binomial log-PDF of a single response 
  * Args: 
  *   y: the response value 
  *   trials: number of trials of the binomial part
  *   p: probability parameter of the binomial part
  *   theta: zero-inflation probability
  * Returns:  
  *   a scalar to be added to the log posterior 
  */ 
  real zero_inflated_binomial_lpmf(int y, int trials, 
  real p, real theta) {
    if (y == 0) { 
      return log_sum_exp(bernoulli_lpmf(1 | theta), 
      bernoulli_lpmf(0 | theta) + 
      binomial_lpmf(0 | trials, p)); 
    } else { 
      return bernoulli_lpmf(0 | theta) +  
      binomial_lpmf(y | trials, p); 
    } 
  }
  
  // copied from R. McElreath (2018). 
  // Statistical rethinking: A Bayesian course with examples in R and Stan. 
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }    
} 

data { 
  int<lower=1> N;                 // total number of observations 
  int Y[N];                       // response variable 
  int total[N];                   // binomial trials
  int<lower=1> K;                 // number of covariates
  matrix[N, K] X;                 // obs-level design matrix 
  int<lower=1> M_2;               // one plant r.e.
  int<lower=1> J_2[N];            // ids for plants
  int<lower=1> N_2;               // num plant spp
  real log_seed_mass[N_2];
  matrix[N_2, N_2] DistP;         // cophenetic distance among plants
}

parameters {
  vector[N_2] r_1_2;          // plqnt-level effects
  //vector[N_2] r_1_t;          // plqnt-level effects
  vector[N_2] r_s;
  vector[N_2] r_st;
  real<lower=0> etasqp;
  real<lower=0> rhosqp;
  real<lower=0> deltap;
  //real alpha;
  vector[K] beta;
  vector[K] betaz;
  //real<lower=0, upper=1> theta;
  real<lower=0> sigma;
  real a;
  real b;
  real az;
  real bz;
  real<lower=0> sigmaz;
}

model {
  matrix[N_2,N_2] SIGMAP;
  vector[N_2] mu;
  vector[N_2] muz;
  
  rhosqp ~ exponential( 0.5 );
  etasqp ~ exponential( 2 );
  deltap ~ normal(0, 2.5);
  
  //theta ~ beta(2,2);
  a ~ normal(0,2);
  b ~ normal(0,1);
  az ~ normal(0,2);
  bz ~ normal(0,1);
  //alpha ~ normal(0,2);
  beta ~ normal(0,1);
  betaz ~ normal(0,1);
  
  
  for(i in 1:N_2){
    mu[i] = a + b * log_seed_mass[i];
    muz[i] = az + bz *  log_seed_mass[i];
  } 
  
  SIGMAP = cov_GPL2(DistP, etasqp, rhosqp, deltap);
  r_1_2 ~ multi_normal( rep_vector(0,N_2) , SIGMAP );
  r_s ~ normal(mu, sigma);
  r_st ~ normal(muz, sigmaz);
  {
    vector[N] p;
    vector[N] theta;
    
    for (n in 1:N){
      p[n] = inv_logit(X[n,] * beta + r_1_2[J_2[n]] + r_s[J_2[n]]);
      theta[n] = inv_logit(X[n,] * betaz + r_st[J_2[n]]);
      target += zero_inflated_binomial_lpmf(Y[n]| total[n], p[n], theta[n]);
    }
    
  }
}

generated quantities{
  int y_pred[N];
  vector[N] mus;
  vector[N] thetas;
  vector[N] res;
  vector[N] resp;
  real sr;
  real sr_rep;
  vector[N] log_lik;
  
  for (n in 1:N){
    mus[n] = inv_logit(X[n,] * beta + r_1_2[J_2[n]] + r_s[J_2[n]]);
    thetas[n] = inv_logit(X[n,] * betaz + r_st[J_2[n]]);
    //if(Y[n]==0)
    y_pred[n] = bernoulli_rng(1-thetas[n]) * binomial_rng(total[n], mus[n]);
    //else
    //y_pred[n] = binomial_rng(total[n], mus[n]);
    log_lik[n] = zero_inflated_binomial_lpmf(Y[n] | total[n], mus[n], thetas[n]);
    res[n] = Y[n] - (mus[n] * total[n]);
    resp[n] = y_pred[n] - (mus[n] * total[n]);
  }
  sr = sum(res);
  sr_rep = sum(resp);
}
