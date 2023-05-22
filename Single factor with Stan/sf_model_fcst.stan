functions {
	real integral_func(real x,   
				real xc,    
                    array[] real theta,
			data array[] real x_r,
			data array[] int x_i){

		real pd = theta[1];
		real corr = theta[2];
		int df = x_i[1];
		int nc = x_i[2];
		real pi;
		
		pi = Phi((inv_Phi(pd)-sqrt(corr)*x)/sqrt(1-corr));
		return exp(binomial_lpmf(df|nc, pi)*normal_lpdf(x|0,1));
	}

	real binom_vasicek_lpmf(int df, int creds, real pd, real corr, data array[] real x_r, data array[] int x_i){

		real aux1;

		array[2] int x_i_int;

		x_i_int[1] = df;
		x_i_int[2] = creds;



		aux1= integrate_1d(integral_func,
					-4.0,
					4.0,
					{pd,corr},
					x_r,
					x_i_int);

		return log(aux1);
	}

	real vasicek_lpdf(real x, real pd, real corr){
		return 1.0/2 * (log(1-corr) - log(corr)) + 0.5*inv_Phi(x)^2 - 1.0/(2*corr)*(inv_Phi(x)*sqrt(1-corr)-inv_Phi(pd))^2;
	}
	
	real vasicek_rng(real pd, real corr){
		real z = normal_rng(0,1);
		return Phi((inv_Phi(pd)-sqrt(corr)*z)/sqrt(1-corr));
	}


	real binom_vasicek2_lpmf(int df, int creds, real pd, real corr){
		real pi;

		
		return binomial_lpmf(df|creds,pi);
		
	}

	
	

}

data {
  int<lower=0> N_observations;
  array[N_observations] int defaults; 
  array[N_observations] int credits; 
  int<lower=0> credit_fcst;
  real a;
  real mu_r;
  real mu_p;
  real phi;
}

transformed data {
  array[0] real x_r;
  array[2] int x_i;
}

parameters {
  real<lower = 0, upper = 1> pd;
  real<lower = 0.01, upper = 0.99> corr;
  array[N_observations] real<lower=0,upper=1> pi;
}

model {

	corr ~ beta_proportion(mu_r,phi);
	pd ~ beta_proportion(mu_p, corr*a);
	

	for(n in 1:N_observations){
		
		target += vasicek_lpdf(pi[n]|pd,corr);
	}

	defaults ~ binomial(credits,pi);
}

generated quantities {
	real pi_aux;
	array[N_observations] int<lower=0,upper=credits> Dreg;
	real pi_fcst;
	int<lower=0,upper=credit_fcst> D_fcst;

	for(n in 1:N_observations){
		
		Dreg[n] = binomial_rng(credits[n],pi[n]);
	}

	pi_fcst = vasicek_rng(pd,corr);
	D_fcst = binomial_rng(credit_fcst, pi_fcst);
}

