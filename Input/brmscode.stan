data { 
  int<lower=1> N;  // total number of observations 
  int Y[N];  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // population-level design matrix 
  // data for group-level effects of ID 1 
  int<lower=1> J_1[N]; 
  int<lower=1> N_1; 
  int<lower=1> M_1; 
  // cholesky factor of known covariance matrix 
  matrix[N_1, N_1] Lcov_1; 
  vector[N] Z_1_1; 
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc; 
  matrix[N, K - 1] Xc;  // centered version of X 
  vector[K - 1] means_X;  // column means of X before centering 
  Kc = K - 1;  // the intercept is removed from the design matrix 
  for (i in 2:K) { 
    means_X[i - 1] = mean(X[, i]); 
    Xc[, i - 1] = X[, i] - means_X[i - 1]; 
  } 
} 
parameters { 
  vector[Kc] b;  // population-level effects 
  real temp_Intercept;  // temporary intercept 
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations 
  vector[N_1] z_1[M_1];  // unscaled group-level effects 
} 
transformed parameters { 
  // group-level effects 
  vector[N_1] r_1_1; 
  r_1_1 = sd_1[1] * (Lcov_1 * z_1[1]); 
} 
model { 
  vector[N] mu; 
  mu = Xc * b + temp_Intercept; 
  for (n in 1:N) { 
    mu[n] = mu[n] + (r_1_1[J_1[n]]) * Z_1_1[n]; 
  } 
  // prior specifications 
  b ~ normal(0, 5); 
  temp_Intercept ~ normal(0, 5); 
  sd_1 ~ student_t(3, 0, 5); 
  z_1[1] ~ normal(0, 1); 
  // likelihood contribution 
  if (!prior_only) { 
    Y ~ bernoulli_logit(mu); 
  } 
} 
generated quantities { 
  real b_Intercept;  // population-level intercept 
  b_Intercept = temp_Intercept - dot_product(means_X, b); 
} 
