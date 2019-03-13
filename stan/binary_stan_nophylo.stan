data { 
  int<lower=0> N; 
  int<lower=0,upper=1> y[N];
  vector[N] pol;
  vector[N] flotime;
  vector[N] minP;

} 
parameters {
  real alpha;
  real  b_pol;
  real b_flotime;
  real b_minP;
}

transformed parameters {
real p[N];
for (i in 1:N)
p[i]=1/(1+exp(-(alpha+b_pol*pol[i]+b_flotime*flotime[i]+b_minP*minP[i])));
}
    
model {
  alpha~ normal(0,10);
  b_pol ~ normal(0,10);
  b_flotime~normal(0,10);
  b_minP~normal(0,10);

  
  y~bernoulli(p);
    
    }
