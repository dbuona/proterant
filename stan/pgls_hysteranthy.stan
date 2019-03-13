//make a pgls bernoulli model 3-13-2019

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
  matrix Vsigma;
}

//transformed parameters {
//real p[N];

//for (i in 1:N)
//p[i]=1/(1+exp(-(alpha+b_pol*pol[i]+b_flotime*flotime[i]+b_minP*minP[i])));

//}
    
model {
  real yhat[N];
  alpha~ normal(0,1);
  b_pol ~ normal(0,5);
  b_flotime~normal(0,1);
  b_minP~normal(0,1);


yhat= 1/(1+exp(-(alpha+b_pol*pol+b_flotime*flotime+b_minP*minP)));

y~multi_normal(yhat, Vsigma);
    }

