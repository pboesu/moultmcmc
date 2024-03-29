//
// This Stan program defines the Underhill-Zucchini Type 3 repeated measures model
// covariates need to be on the individual level, not the observation model
//
data {
  int<lower=0,upper=1> flat_prior;//translated logical use flat prior on start and duration?
  int<lower=0> N_ind;//number of recaptured individuals
  int<lower=0> N_ind_rep;//number of recaptured individuals
  int<lower=0> N_moult;//J
  int<lower=0> Nobs_replicated;//number of observations from individuals with repeat measures
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0>individual[N_moult];//individual identifier
  int<lower=0>individual_first_index[N_ind];//row first occurrence of each individual in the model frame
  int<lower=0>replicated[Nobs_replicated];//indices of obs from individuals with repeat measures
  int<lower=0>not_replicated[N_moult - Nobs_replicated];//indices of obs from individuals without repeat measures
  int<lower=0>is_replicated[N_ind];
  int<lower=0>replicated_individuals[N_ind_rep];//individual id's that are replicated
  int<lower=0>replicated_id_obs[N_moult];//new running id from 1:N_ind_rep to index the individual start date estimates. Unreplicated individuals receive a zero value
  //int<lower=0> N_new;//K
  //vector[N_new] new_dates;//v_k
  //predictors
  int N_pred_mu;//number of predictors for start date
  matrix[N_moult,N_pred_mu] X_mu;//design matrix for start date NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_tau;//number of predictors for duration
  matrix[N_moult,N_pred_tau] X_tau;//design matrix for duration NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_sigma;//number of predictors for start date sigma
  matrix[N_moult,N_pred_sigma] X_sigma;//design matrix for sigma start NB: when forming design matrix must paste together responses in blocks old, moult, new
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {

  //real<lower=0> tau;//moult duration
  //real<lower=0> sigma;//population varianmoult start date
  vector[N_pred_mu] beta_mu;//regression coefficients for start date
  vector[N_pred_tau] beta_tau;//regression coefficients for duration
  vector[N_pred_sigma] beta_sigma;//regression coefficients for sigma start date
  vector[N_ind_rep] mu_ind;//individual effect on start date
  vector[N_ind] tau_ind;//individual effect on duration
  real<lower=0> sigma_mu_ind;//?residual variance in regression of score on date within individuals
  //real<lower=0> sigma_tau;//pop std dev of duration
}

transformed parameters{
  real sigma_intercept = exp(beta_sigma[1]);
  //post-sweep random effects
  real beta_star = beta_mu[1] + mean(mu_ind);
  vector[N_ind_rep] mu_ind_star = mu_ind- mean(mu_ind);
  real finite_sd = sd(mu_ind_star);
}

// The model to be estimated.
model {
  //vector[N_old] P;
  vector[N_moult] q;
  vector[N_moult] Q;
  vector[N_moult] mu;//start date lin pred
  vector[N_moult] tau;//duration lin pred
  vector[N_moult] sigma;//duration lin pred

  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);//use log link for variance lin pred

//for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
//print(P);
for (i in 1:N_moult) {
  if (is_replicated[individual[i]] == 1) {
    q[i] = normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i] - (mu[i] + mu_ind[replicated_id_obs[i]])) | 0 , sigma_mu_ind);//replicated individuals
  } else {
    q[i] = log(tau[i]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i]) | mu[i], sigma[i]);//unreplicated - don't estimate an individual intercept??
  }
}
//for (i in 1:N_moult) q[i] = normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i]) | mu[i], sigma_mu_ind);//N.B. unlike P and R this returns a log density
for (i in 1:N_moult) {
  //if (is_replicated[individual[i]] == 1) { //this appears to break initialisation
    Q[i] = Phi((moult_dates[i] - (mu[i]))/sigma[i]) - Phi((moult_dates[i] - tau[i] - (mu[i]))/sigma[i]);
  //}
  }
//individual start dates are drawn from the population distribution of start dates - TODO: this is slow, so don't calculate it for replicated obs
mu_ind ~ normal(0, sigma[individual_first_index][replicated_individuals]);//only estimate this for replicated individuals? shouldn't this be mu[i] for replicated individuals??
//tau_ind[replicated_individuals] ~ normal(0, sigma_tau);//only estimate this for replicated individuals? shouldn't this be mu[i] for replicated individuals??
//print(R);
//print(sum(log(P)));
//print(sum(log(Q)));
//print(sum(log(R)));
target += sum(q) - sum(log(Q[not_replicated]));
//priors
if (flat_prior == 1) {
 beta_mu[1] ~ uniform(0,366);
 beta_tau[1] ~ uniform(0,366);
  } else {
 beta_mu[1] ~ normal(150,50)T[0,366];
 beta_tau[1] ~ normal(100,30)T[0,366];
}
beta_sigma[1] ~ normal(0,5);
sigma_mu_ind ~ normal(0,1);
//sigma_tau ~ normal(0,1);
}

generated quantities{
//TODO: calculate individual intercepts for all observed individuals here?
// //NB: code duplication for the likelihood calculation is less than ideal - refactor to a use stan function?
// //real end_date;
// // end_date = mu + tau;
//
// vector[N_moult] log_lik;
// vector[N_moult] q;
// vector[N_moult] Q;
vector[N_pred_mu] beta_mu_out;//regression coefficients for start date beta_mu_out
vector[N_moult] mu;//start date lin pred
vector[N_moult] tau;//duration lin pred
vector[N_moult] sigma;//duration lin pred
vector[N_ind_rep] mu_ind_out;//individual intercepts for output
//vector[N_ind] tau_ind_out;//individual intercepts for output

  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);

mu_ind_out = mu_ind_star + beta_star;

//The following code calculates individual intercepts for all individuals, this should be left to a predict method, and not be done during estimation
//   for (i in 1:N_moult) {
//   if (is_replicated[individual[i]] == 1) {
//      mu_ind_out[individual[i]] = mu_ind[individual[i]] + mu[i];//replicated individuals
//      //tau_ind_out[individual[i]] = tau_ind[individual[i]] + tau[i];//replicated individuals
//   } else {
//     mu_ind_out[individual[i]] = moult_dates[i] - moult_indices[i]*tau[i];//unreplicated
//     //tau_ind_out[individual[i]] = tau[i];
//   }
// }
  if (N_pred_mu > 1){
    beta_mu_out = append_row(beta_star,beta_mu[2:N_pred_mu]);// collate post-swept intercept with remaining
  } else {
    beta_mu_out[1] = beta_star;// intercept-only model
  }


// for (i in 1:N_moult) q[i] = log(tau[i]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i]) | mu[i], sigma[i]);//N.B. unlike Q this returns a log density
// for (i in 1:N_moult) Q[i] = Phi((moult_dates[i] - mu[i])/sigma[i]) - Phi((moult_dates[i] - tau[i] - mu[i])/sigma[i]);
//
// log_lik = (q - log(Q));
}
