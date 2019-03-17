//make a pgls bernoulli model 3-13-2019
//https://gist.github.com/mvuorre/839c382e1add301d75cd

data { 
  int<lower=0> N; 
  int<lower=0,upper=1> y;
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
vector[N ]yhat;
vector[N] theta;
matrix[N,N] Vlambda; // VCV x lambda??
matrix[N,N] Vsigma; // (VCVxlambda) x sigma??
real detV; // determinant of the VCV

Vlambda = (lambda*Lmat + Ident) .* V;
Vsigma = sigma^2*Vlambda;


yhat= 1/(1+exp(-(alpha+b_pol*pol+b_flotime*flotime+b_minP*minP)));


 //with the sigma of the phylogeny

}
    
model {
  
  alpha~ normal(0,5);
  b_pol ~ normal(0,5);
  b_flotime~normal(0,5);
  b_minP~normal(0,5);



theta~multi_normal(yhat,Vsigma);
y~bernoulli_logit(theta);
    }

