//
// This Stan program defines the Underhill-Zucchini Type 5 repeated measures model
//  covariates need to be on the individual level not the observation model for logical consistncy

data {
  //responses
  int<lower=0,upper=1> flat_prior;//translated logical use flat prior on start and duration?
  real beta_sd;//sd for non-flat prior on regression coeficients
  int<lower=0> N_ind;//number of individuals
  int<lower=0> N_ind_rep;//number of recaptured individuals
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  int<lower=0> N_new;//K
  vector[N_new] new_dates;//v_k
  int<lower=0> Nobs_replicated;//number of observations from individuals with repeat measures
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0>individual[N_moult+N_old+N_new]; //individual identifier
  int<lower=0>individual_first_index[N_ind];//row first occurrence of each individual in the model frame
  //int<lower=0>replicated[Nobs_replicated];//indices of obs from individuals with repeat measures//not used
  //int<lower=0>not_replicated[(N_moult + N_old + N_new) - Nobs_replicated];//indices of obs from individuals without repeat measures//not use
  int<lower=0>is_replicated[N_ind];
  int<lower=0>replicated_individuals[N_ind_rep];//individual id's that are replicated - i.e. indices of the random_effect intercept
  //int<lower=0>replicated_id[N_ind];//new running id from 1:N_ind_rep to index the individual start date estimates. Unreplicated individuals receive a zero value
  int<lower=0>replicated_id_obs[N_moult+N_old+N_new];//new running id from 1:N_ind_rep to index the individual start date estimates. Unreplicated individuals receive a zero value
  int<lower=0,upper=1> lumped; //indicator variable for t2l likelihood
  int<lower=0,upper=1> llik; //indicator variable for calculating and returning pointwise log-likelihood
  int<lower=0,upper=1> use_phi_approx; //indicator variable for likelihood
  int<lower=0,upper=1> same_sigma; //indicator variable
  //predictors
  int N_pred_mu;//number of predictors for start date
  matrix[N_old+N_moult+N_new,N_pred_mu] X_mu;//design matrix for start date NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_tau;//number of predictors for duration
  matrix[N_old+N_moult+N_new,N_pred_tau] X_tau;//design matrix for duration NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_sigma;//number of predictors for start date sigma
  matrix[N_old+N_moult+N_new,N_pred_sigma] X_sigma;//design matrix for sigma start NB: when forming design matrix must paste together responses in blocks old, moult, new
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {

  //real<lower=0> tau;//moult duration
  //real<lower=0> sigma;//population varianmoult start date
  vector[N_pred_mu] beta_mu;//regression coefficients for start date
  vector[N_pred_tau] beta_tau;//regression coefficients for duration
  vector[N_pred_sigma] beta_sigma;//regression coefficients for sigma start date
  vector[N_ind_rep] mu_ind;//individual effect on sigma start date
  real<lower=0> sigma_mu_ind;//?residual variance in regression of score on date within individuals
}

transformed parameters{
  //all moved to generated quantities
}

// The model to be estimated.
model {
  vector[N_old] P;
  vector[N_moult] q;
  vector[N_new] R;
  vector[N_old+N_moult+N_new] mu;//start date lin pred
  vector[N_old+N_moult+N_new] tau;//duration lin pred
  vector[N_old+N_moult+N_new] sigma;//duration lin pred

  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);//use log link for variance lin pred

if (lumped == 0){
  for (i in 1:N_old) {
	  if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	  if(use_phi_approx == 0){
	    P[i] = 1 - Phi((old_dates[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind);
	  } else {
	    P[i] = 1 - Phi_approx((old_dates[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind);
	  }

	  } else {//standard likelihood for Type 2 model
      P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
	  }
  }
} else {//lumped likelihood
  for (i in 1:N_old) {
	  if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	    P[i] = (1 - Phi((old_dates[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind)) + Phi((old_dates[i] - tau[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind);
	  } else {//standard likelihood for Type 2 model
      P[i] = (1 - Phi((old_dates[i] - mu[i])/sigma[i])) + Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
	  }
  }
}

for (i in 1:N_moult){
  if (is_replicated[individual[i + N_old]] == 1) {
   q[i] = normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old] + mu_ind[replicated_id_obs[i + N_old]], sigma_mu_ind);//replicated individuals. NB - indexing looks messy because i runs from 1:N_moult, but the function uses both vectors of the total dataset (1:(N_old+N_moult) and the moult dataset (1:N_moult))
  } else {
      q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
  }
}

if (lumped == 0){
  for (i in 1:N_new) {
	  if (is_replicated[individual[i + N_old + N_moult]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	  //print(i);
	  //print(is_replicated[individual[i]]);
	  //print(replicated_id_obs[i + N_old + N_moult]);
	    R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - (mu[i + N_old + N_moult] + mu_ind[replicated_id_obs[i + N_old + N_moult]]))/sigma_mu_ind);
	  } else {//standard likelihood for Type 2 model
      R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]);
	  }
  }
} else {//lumped likelihood
  for (i in 1:N_new) {
	  if (is_replicated[individual[i + N_old + N_moult]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	    R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - (mu[i + N_old + N_moult] + mu_ind[replicated_id_obs[i + N_old + N_moult]]))/sigma_mu_ind) + (1 - Phi((new_dates[i] - (mu[i + N_moult + N_old] + mu_ind[replicated_id_obs[i + N_moult + N_old]]))/sigma_mu_ind));
	  } else {//standard likelihood for Type 2 model
      R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]) + (1 - Phi((new_dates[i] - (mu[i + N_moult + N_old]))/sigma[i + N_moult + N_old]));
	  }
  }
}

mu_ind ~ normal(0, sigma[individual_first_index][replicated_individuals]);//

target += sum(log(P))+sum(q)+sum(log(R));

//priors
if (flat_prior == 1) {
 beta_mu[1] ~ uniform(-366,366);
 beta_tau[1] ~ uniform(0,366);
} else {
 beta_mu[1] ~ normal(150,50)T[-366,366];
 beta_tau[1] ~ normal(100,30)T[0,366];
}
 if (beta_sd > 0){//messy implementation, better to do flat priors by default, non-flat priors by explicit values only?
  if (N_pred_mu > 1){
   for (i in 2:N_pred_mu){
     beta_mu[i] ~ normal(0,beta_sd);
   }
  }
  if (N_pred_tau > 1){
   for (i in 2:N_pred_tau){
     beta_tau[i] ~ normal(0,beta_sd);
   }
  }
  if (N_pred_sigma > 1){
   for (i in 2:N_pred_sigma){
     beta_sigma[i] ~ normal(0,beta_sd);
   }
  }
}
beta_sigma[1] ~ normal(0,2);// on log link scale!
sigma_mu_ind ~ normal(0,10);
}

generated quantities{
  vector[N_pred_mu] beta_mu_out;//regression coefficients for start date beta_mu_out
  vector[N_ind_rep] mu_ind_out;//individual intercepts for output
  real sigma_intercept = exp(beta_sigma[1]);//transformed sigma intercept for output

  //post-sweep random effects
  real beta_star = beta_mu[1] + mean(mu_ind);
  vector[N_ind_rep] mu_ind_star = mu_ind - mean(mu_ind);
  real finite_sd = sd(mu_ind_star);
  vector[(N_old+N_moult+N_new)*llik] log_lik;//log_lik - make zero length when not requested



  if (N_pred_mu > 1){
    beta_mu_out = append_row(beta_star,beta_mu[2:N_pred_mu]);// collate post-swept intercept with remaining
  } else {
    beta_mu_out[1] = beta_star;// intercept-only model
  }

  mu_ind_out = mu_ind_star + beta_star;



 if (llik == 1){
  vector[N_old] P;
  vector[N_moult] q;
  vector[N_new] R;
  vector[N_old+N_moult+N_new] mu;//start date lin pred
  vector[N_old+N_moult+N_new] tau;//duration lin pred
  vector[N_old+N_moult+N_new] sigma;//duration lin pred


  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);//use log link for variance lin pred

if (lumped == 0){
  for (i in 1:N_old) {
	  if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	  if(use_phi_approx == 0){
	    P[i] = 1 - Phi((old_dates[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind);
	  } else {
	    P[i] = 1 - Phi_approx((old_dates[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind);
	  }

	  } else {//standard likelihood for Type 2 model
      P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
	  }
  }
} else {//lumped likelihood
  for (i in 1:N_old) {
	  if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	    P[i] = (1 - Phi((old_dates[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind)) + Phi((old_dates[i] - tau[i] - (mu[i] + mu_ind[replicated_id_obs[i]]))/sigma_mu_ind);
	  } else {//standard likelihood for Type 2 model
      P[i] = (1 - Phi((old_dates[i] - mu[i])/sigma[i])) + Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
	  }
  }
}

for (i in 1:N_moult){
  if (is_replicated[individual[i + N_old]] == 1) {
   q[i] = normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old] + mu_ind[replicated_id_obs[i + N_old]], sigma_mu_ind);//replicated individuals. NB - indexing looks messy because i runs from 1:N_moult, but the function uses both vectors of the total dataset (1:(N_old+N_moult) and the moult dataset (1:N_moult))
  } else {
      q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
  }
}

if (lumped == 0){
  for (i in 1:N_new) {
	  if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	    R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - (mu[i + N_old + N_moult] + mu_ind[replicated_id_obs[i + N_old + N_moult]]))/sigma_mu_ind);
	  } else {//standard likelihood for Type 2 model
      R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]);
	  }
  }
} else {//lumped likelihood
  for (i in 1:N_new) {
	  if (is_replicated[individual[i]] == 1) {//longitudinal tobit-like likelihood (this only makes sense if within year recaptures contain at least one active moult score?!)
	    R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - (mu[i + N_old + N_moult] + mu_ind[replicated_id_obs[i + N_old + N_moult]]))/sigma_mu_ind) + (1 - Phi((new_dates[i] - (mu[i + N_moult + N_old] + mu_ind[replicated_id_obs[i + N_moult + N_old]]))/sigma_mu_ind));
	  } else {//standard likelihood for Type 2 model
      R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]) + (1 - Phi((new_dates[i] - (mu[i + N_moult + N_old]))/sigma[i + N_moult + N_old]));
	  }
  }
}

//mu_ind[replicated_individuals] ~ normal(0, sigma[individual_first_index][replicated_individuals]);//

log_lik = append_row(append_row(log(P),q),log(R));
    //log_lik = append_row((log(P) - log1m(Rt)), (q - log1m(Ru)));
  }
}
