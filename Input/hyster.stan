data{
	int N; // number of spp 
	int K; // 5 (WTF?) number of predictors?
	row_vector[N] y; // response pro 
	matrix[N, N] V; // VCV
	matrix[N, N] Lmat; // some sort of pre-designed matrix with ones and zero-diagonal
	row_vector[N] pol; // predictor1 
	row_vector[N] class2; // predictor2
	row_vector[N] shade_bin; // predictor3
	row_vector[N] fruit_bin; // predictor4 
	row_vector[N] flo_type; // predictor5
}

transformed data{
	real Ndiv; // a real number (holder)
	matrix[N, N] Ident; // a N by N matrix (holder)
	matrix[N,N] X;
  
	Ndiv = N; // fill out Ndiv
	Ident = diag_matrix(rep_vector(1, N)); // matrix of 0s with 1s on diagonal
	X[1] = pol;
	X[2] = class2;
	X[3] = shade_bin;
	X[4] = fruit_bin;
	X[5] = flo_type;
}

parameters{
	real<lower=0> sigma;
	real<lower=0, upper=1> lambda;
	vector[N] B;
	real alpha;
}


model{
	matrix [N,N] L_Vsigma; // (VCVxlambda) x sigma?? I think this cant be a mtrix
	vector [N] yhat; // holder for predicted values
//	real detV; // determinant of the VCV
//	real cdf;  // cumulative distribution f(x)??
//	real logLike_PGLS; // a log likelihood!

  {
    matrix[N, N] Vsigma;
    for (i in 1:(N-1)) {
      for (j in (i + 1):N) 
        Vsigma[i,j] = V[i,j] * lambda; //changed this from little v since it say there is no V
      Vsigma[i,i] = sigma^2;
    }
    L_Vsigma = (Vsigma);//cholesky_decompose
  }
	yhat = alpha + X*B;


	//detV <- log_determinant(Vsigma);
	//cdf <- ((y-yhat)'*inverse(Vsigma)*(y-yhat));
	//logLike_PGLS <-  -0.5*(detV + cdf);
	//increment_log_prob(logLike_PGLS);
	
	y ~ multi_normal_cholesky(yhat,L_Vsigma); // is this right? 
	
	B ~ normal(0, 1); //deal with these once it runs
	sigma ~ cauchy(0, 2.5);
	lambda ~ beta(1, 1);
	alpha ~ normal(0, 1);
	
}
