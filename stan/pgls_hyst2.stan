data { 
  int<lower=0> N; // rows of data
  int<lower=0,upper=1> y[N]; //0,1 respose variable
  real pol[N]; // polination syndrome 0,1
  real flotime[N]; // flowering time
  real minP[N]; // minimum repcitp
  
  matrix[N, N] V; //  VCV
	matrix[N, N] Lmat; // some sort of pre-designed matrix

} 

transformed data{
	real Ndiv; // a real number (holder)
	matrix[N, N] Ident; // a N by N matrix (holder)

	Ndiv = N; // fill out Ndiv
	Ident = diag_matrix(rep_vector(1, N)); // matrix of 0s with 1s on diagonal
}

parameters {
 real alpha; // probability if all treatments are 0
  real b_pol; // slope on above with pol syndrom
  real b_flotime; //slope on above
  real b_minP; //slope on above
  real<lower=0> sigma; //sigma
	real<lower=0, upper=1> lambda; //phylogenetic signal. This probably needs to be phyo D.
}

transformed parameters {
vector[N] yhat; // the theota of bernoulli (prob y=1) with out accounting for phylogeny


for (i in 1:N)
yhat[i]= 1/(1+exp(-(alpha+b_pol*pol[i]+b_flotime*flotime[i]+b_minP*minP[i]))); //(logistic transformation of linear model)
}

model {
  vector[N] theta;
  matrix[N,N] Vlambda; // VCV x lambda??
	matrix[N,N] Vsigma; // (VCVxlambda) x sigma??
	real detV; // determinant of the VCV
	
		
  alpha~ normal(0,5);
  b_pol ~ normal(0,5);
  b_flotime~normal(0,5);
  b_minP~normal(0,5);

  Vlambda = (lambda*Lmat + Ident) .* V;
	Vsigma = sigma^2*Vlambda;

theta~multi_normal(yhat,Vsigma);//I think is makes a new probability that y=1 account for the phylogeny
target += multi_normal_lpdf(theta |yhat, Vsigma);
//Chain 1: Exception: multi_normal_lpdf: Random variable[1] is nan, but must not be nan!  (in 'model16b11c3826a7_pgls_hyst2' at line 53)

y~bernoulli_logit(theta); //final model
}
