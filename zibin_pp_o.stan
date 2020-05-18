
functions { 
  
    /* compute the kronecker product
  * Copied from brms: Paul-Christian Bürkner (2018). 
  * Advanced Bayesian Multilevel Modeling with the R Package brms. 
  * The R Journal, 10(1), 395-411. <doi:10.32614/RJ-2018-017>
  * Args: 
    *   A,B: matrices 
  * Returns: 
    *   kronecker product of A and B
  */ 
    matrix kronecker(matrix A, matrix B) { 
      matrix[rows(A)*rows(B), cols(A)*cols(B)] kron; 
      for (i in 1:cols(A)) { 
        for (j in 1:rows(A)) { 
          kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
        } 
      } 
      return kron; 
    } 
  
  
    // copied from https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_zero_inflated_binomial.stan
    /* zero-inflated binomial log-PDF of a single response 
   * Args: 
   *   y: the response value 
   *   trials: number of trials of the binomial part
   *   theta: probability parameter of the binomial part
   *   zi: zero-inflation probability
   * Returns:  
   *   a scalar to be added to the log posterior 
   */ 
  real zero_inflated_binomial_lpmf(int y, int trials, 
                                   real theta, real zi) {
    if (y == 0) { 
      return log_sum_exp(bernoulli_lpmf(1 | zi), 
                         bernoulli_lpmf(0 | zi) + 
                         binomial_lpmf(0 | trials, theta)); 
    } else { 
      return bernoulli_lpmf(0 | zi) +  
             binomial_lpmf(y | trials, theta); 
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
  //int<lower=1> M_2;               // one plant r.e.
  //int<lower=1,upper=N_1> J_1[N];  // ids for sites
  int<lower=1> J;                 // num of groups (species)
  int<lower=1,upper=J> jj[N];     // group id 
  int<lower=1> L;                 // num group level predictors
  matrix[J,L] TT;                 // group-level traits
  matrix[J,J] C;                  // phylogenetic correlation matrix
  vector[J] ones;                 // vector on 1s
  //matrix[N_2, N_2] DistP;         // cophenetic distance among plants
  //matrix[N_1,N_1] Dmat;           // sites distance matrix
}

parameters {
  //corr_matrix[K] Omega;
  cholesky_factor_corr[K] Lcorr; 
  vector<lower=0>[K] tau;
  vector[J * K] beta;
  real<lower=0,upper=1> rho;
  vector[L * K] z;
  real<lower=0, upper=1> zi;
}

transformed parameters {
  matrix[K, K] Omega = multiply_lower_tri_self_transpose(Lcorr);
  matrix[K, K] Sigma = quad_form_diag(Omega, tau);
  matrix[J*K, J*K] S = kronecker(Sigma, rho * C + (1-rho) * diag_matrix(ones));
  matrix[L, K] Z = to_matrix(z, L, K);
  vector[J * K] m = to_vector(TT * Z);
  matrix[J, K] b_m = to_matrix(beta, J, K);
} 

model {
  //matrix[N_1,N_1] SIGMA;
  //matrix[N_2,N_2] SIGMAP;
  
  //rhosq ~ exponential( 0.5 );
  //etasq ~ exponential( 2 );
  //delta ~ normal(0, 2.5);
  //rhosqp ~ exponential( 0.5 );
  //etasqp ~ exponential( 2 );
  //deltap ~ normal(0, 2.5);
    
  //Omega ~ lkj_corr(2);
  Lcorr ~ lkj_corr_cholesky(1);
  tau ~ student_t(3,0,10); // cauchy(0, 2.5);
  beta ~ multi_normal(m, S);
  zi ~ beta(2,2);

 
 {
  vector[N] theta;
  
  target += log_sum_exp(log(0.5) +  beta_lpdf(rho|1, 100), log(0.5) +  beta_lpdf(rho|2, 2));

  for (n in 1:N){
    theta[n] = inv_logit(b_m[jj[n],1]+X[n,2]*b_m[jj[n],2]*X[n,3]*b_m[jj[n],3]+X[n,4]*b_m[jj[n],4]);
    //theta[n] = inv_logit(sum(X[n,] * b_m[jj[n],]));
    target += zero_inflated_binomial_lpmf(Y[n]| total[n], theta[n], zi);
    
  }

  }
}

// generated quantities{
//   int y_pred[N];
//   vector[N] mus;
//   vector[N] res;
//   vector[N] resp;
//   real sr;
//   real sr_rep;
//   vector[N] log_lik;
// 
//   for (n in 1:N){
//      mus[n] = inv_logit(b_m[jj[n],1]+X[n,2]*b_m[jj[n],2]*X[n,3]*b_m[jj[n],3]+X[n,4]*b_m[jj[n],4]);
//      if(Y[n]==0)
//        y_pred[n] = bernoulli_rng(1-zi) * binomial_rng(total[n], mus[n]);
//      else
//        y_pred[n] = binomial_rng(total[n], mus[n]);
//      log_lik[n] = zero_inflated_binomial_lpmf(Y[n] | total[n], mus[n], zi);
//      res[n] = Y[n] - mus[n];
//      resp[n] = y_pred[n] - mus[n];
//   }
//   sr = sum(res);
//   sr_rep = sum(resp);
// }
