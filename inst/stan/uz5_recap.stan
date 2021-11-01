//
// This Stan program defines the Underhill-Zucchini Type 5 repeated measures model
//  covariates need to be on the individual level not the observation model for logical consistncy

data {
  //responses
  int<lower=0> N_ind;//number of individuals
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  int<lower=0> Nobs_replicated;//number of observations from individuals with repeat measures
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0>individual[N_moult+N_old]; //individual identifier
  int<lower=0>individual_first_index[N_ind];//row first occurrence of each individual in the model frame
  int<lower=0>replicated[Nobs_replicated];//indices of obs from individuals with repeat measures
  int<lower=0>not_replicated[(N_moult + N_old) - Nobs_replicated];//indices of obs from individuals without repeat measures
  int<lower=0>is_replicated[N_ind];
  //int<lower=0> N_new;//K
  //vector[N_new] new_dates;//v_k
  //predictors
  int N_pred_mu;//number of predictors for start date
  matrix[N_old+N_moult,N_pred_mu] X_mu;//design matrix for start date NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_tau;//number of predictors for duration
  matrix[N_old+N_moult,N_pred_tau] X_tau;//design matrix for duration NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_sigma;//number of predictors for start date sigma
  matrix[N_old+N_moult,N_pred_sigma] X_sigma;//design matrix for sigma start NB: when forming design matrix must paste together responses in blocks old, moult, new
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {

  //real<lower=0> tau;//moult duration
  //real<lower=0> sigma;//population varianmoult start date
  vector[N_pred_mu] beta_mu;//regression coefficients for start date
  vector[N_pred_tau] beta_tau;//regression coefficients for duration
  vector[N_pred_sigma] beta_sigma;//regression coefficients for sigma start date
  vector[N_ind] mu_ind;//individual effect on sigma start date
  real<lower=0> sigma_mu_ind;//?residual variance in regression of score on date within individuals
}

transformed parameters{
  real sigma_intercept = exp(beta_sigma[1]);
  //post-sweep random effects
  real beta_star = beta_mu[1] + mean(mu_ind);
  vector[N_ind] mu_ind_star = mu_ind - mean(mu_ind);
  real finite_sd = sd(mu_ind_star);
}

// The model to be estimated.
model {
  vector[N_old] Rt;
  vector[N_old] P;
  vector[N_moult] Ru;
  vector[N_moult] q;
  vector[N_old+N_moult] mu;//start date lin pred
  vector[N_old+N_moult] tau;//duration lin pred
  vector[N_old+N_moult] sigma;//duration lin pred

  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);//use log link for variance lin pred

for (i in 1:N_old) {
	if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	  P[i] = 1 - Phi((old_dates[i] - mu[i] + mu_ind[individual[i]])/sigma_mu_ind);
	} else {//standard likelihood for Type 5 model
    P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
    Rt[i] = Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
	}
}
for (i in 1:N_moult){
  if (is_replicated[individual[i]] == 1) {
   q[i] = normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i]) | mu[i] + mu_ind[individual[i]], sigma_mu_ind);//replicated individuals
  } else {
   Ru[i] = Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old] + + mu_ind[individual[i + N_old]])/sigma[i + N_old]);
   q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old] + mu_ind[individual[i + N_old]], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
  }
}
mu_ind ~ normal(0, sigma[individual_first_index]);//

//print(sum(q));
//print(sum(log(P)));
//print(sum(log(1-Rt)));
//print(sum(log(1-Ru)));
target += sum(log(P) - log1m(Rt)) + sum(q - log1m(Ru)) ;
//priors
beta_mu[1] ~ uniform(0,366);
beta_tau[1] ~ uniform(0,366);
beta_sigma[1] ~ normal(0,5);
sigma_mu_ind ~ normal(0,1);
}

generated quantities{

//NB: code duplication for the likelihood calculation is less than ideal - refactor to a use stan function?
//real end_date;
// end_date = mu + tau;
//
//vector[N_old + N_moult] log_lik;
//vector[N_old] Rt;
//vector[N_old] P;
//vector[N_moult] Ru;
//vector[N_moult] q;
//vector[N_old+N_moult] mu;//start date lin pred
//vector[N_old+N_moult] tau;//duration lin pred
//vector[N_old+N_moult] sigma;//duration lin pred
//
//  mu = X_mu * beta_mu;
//  print(mu);
//  tau = X_tau * beta_tau;
//  print(tau);
//  sigma = exp(X_sigma * beta_sigma);
//
//for (i in 1:N_old) {
  //P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
  //Rt[i] = Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
//}
//for (i in 1:N_moult){
  // Ru[i] = Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old])/sigma[i + N_old]);
  // q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
//}
//log_lik = append_row((log(P) - log1m(Rt)), (q - log1m(Ru)));
}
