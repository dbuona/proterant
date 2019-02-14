data { 
  int<lower=0> N; 
  int<lower=0,upper=1> y[N];
  int<lower=0,upper=1> pol[N];
} 
parameters {
  real<lower=0,upper=1> alpha;
  real b_pol;
  
} 
model {
  alpha~ beta(1,1);
  b_pol ~ beta(1,1);
  for (n in 1:N) 
    y[n] ~ bernoulli(alpha+b_pol*pol[N]);
}
