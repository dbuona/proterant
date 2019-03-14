//make a pgls bernoulli model 3-13-2019

data { 
  int<lower=0> N; 
  int<lower=0,upper=1> y[N];
  vector[N] pol;
  vector[N] flotime;
  vector[N] minP;

	matrix[N, N] V; // will be VCV?
	matrix[N, N] Lmat; // some sort of pre-designed matrix
} 

transformed data{
	real Ndiv; // a real number (holder)
	matrix[N, N] Ident; // a N by N matrix (holder)

	Ndiv = N; // fill out Ndiv
	Ident = diag_matrix(rep_vector(1, N)); // matrix of 0s with 1s on diagonal
}

parameters {
  real alpha;
  real  b_pol;
  real b_flotime;
  real b_minP;
  
  real<lower=0> sigma;
	real<lower=0, upper=1> lambda;
}

transformed parameters {

matrix[N,N] Vlambda; // VCV x lambda??
matrix[N,N] Vsigma; // (VCVxlambda) x sigma??
vector[N] yhat; // holder for predicted values
real detV; // determinant of the VCV
real p[N];
real theta[N];


Vlambda = (lambda*Lmat + Ident) .* V;
Vsigma = sigma^2*Vlambda;

for (i in 1:N)
p[i]=1/(1+exp(-(alpha+b_pol*pol[i]+b_flotime*flotime[i]+b_minP*minP[i])));
for (i in 1:N)
theta[i]=p[i]+Vsigma //with the sigma of the phylogeny
}
    
model {
  real yhat[N];
  alpha~ normal(0,1);
  b_pol ~ normal(0,5);
  b_flotime~normal(0,1);
  b_minP~normal(0,1);


yhat= 1/(1+exp(-(alpha+b_pol*pol+b_flotime*flotime+b_minP*minP)));

y~multi_normal(theta);
    }

