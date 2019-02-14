data{
	int N; // number of spp 
	int K; // 1 (WTF?)
	int<lower=0,upper=1> y[N]; // hyseranthy
  
	matrix[N, N] V; // will be VCV?
	matrix[N, N] Lmat; // some sort of pre-designed matrix
	
	//predictors
	vector[N] flo_time; // try with just flowering time
}
transformed data{
	real Ndiv; // a real number (holder)
	matrix[N, N] Ident; // a N by N matrix (holder)

	Ndiv = N; // fill out Ndiv
	Ident = diag_matrix(rep_vector(1, N)); // matrix of 0s with 1s on diagonal
}
parameters{
	real<lower=0> sigma;
	real<lower=0, upper=1> lambda;
	vector[2] B;
}
transformed parameters{
}
model{
	matrix[N,N] Vlambda; // VCV x lambda??
	matrix[N,N] Vsigma; // (VCVxlambda) x sigma??
	vector[N] yhat; // holder for predicted values
	real detV; // determinant of the VCV
//	real cdf;  // cumulative distribution f(x)??
//	real logLike_PGLS; // a log likelihood!

	Vlambda = (lambda*Lmat + Ident) .* V;
	Vsigma = sigma^2*Vlambda;

  for(i in 1:N)
	yhat[i] = B[1] + B[2]*flo_time[i];


	//detV <- log_determinant(Vsigma);
	//cdf <- ((y-yhat)'*inverse(Vsigma)*(y-yhat));
	//logLike_PGLS <-  -0.5*(detV + cdf);
	//increment_log_prob(logLike_PGLS);
	
	y ~ bernoulli(yhat,Vsigma); // Cant parse becauseof vsig

	B ~ normal(0, 1);
	sigma ~ cauchy(0, 2.5);
	lambda ~ beta(1, 1);
}
