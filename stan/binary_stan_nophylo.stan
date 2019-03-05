data { 
  int<lower=0> N; 
  int<lower=0,upper=1> y[N];
  vector[N] pol;
  vector[N] flo;
  vector[N] minp;
} 
parameters {
  real alpha;
  real  b_pol;
  real  b_flo;
  real  b_drought;
} 


    

model {
  alpha~ uniform(0,1);
  b_pol ~ normal(0,1);
  b_flo ~ normal(0,1);
  b_drought ~normal(0,1);
  
  y~bernoulli_logit(alpha+b_pol*pol+b_flo*flo+b_drought*minp);
    
    }
