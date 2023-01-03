//
// This Stan program defines the Underhill-Zucchini Type 1 model
//

// The input data is a vector 'y' of length 'N'.
data {
  //responses
  int<lower=0,upper=1> flat_prior;//translated logical use flat prior on start and duration?
  real beta_sd;//sd for non-flat prior on regression coeficients
  int<lower=0> N_ind;//number of individuals
  int<lower=0> N_ind_rep;//number of recaptured individuals
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  vector[N_moult] moult_dates;//u_j
  int<lower=0> N_new;//K
  vector[N_new] new_dates;//v_k
  int<lower=0> Nobs_replicated;//number of observations from individuals with repeat measures
  int<lower=0>individual[N_moult+N_old+N_new]; //individual identifier
  int<lower=0>individual_first_index[N_ind];//row first occurrence of each individual in the model frame
  int<lower=0>is_replicated[N_ind];
  int<lower=0>replicated_individuals[N_ind_rep];//individual id's that are replicated - i.e. indices of the random_effect intercept
  //predictors
  int N_pred_mu;//number of predictors for start date
  matrix[N_old+N_moult+N_new,N_pred_mu] X_mu;//design matrix for start date NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_tau;//number of predictors for duration
  matrix[N_old+N_moult+N_new,N_pred_tau] X_tau;//design matrix for duration NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_sigma;//number of predictors for start date sigma
  matrix[N_old+N_moult+N_new,N_pred_sigma] X_sigma;//design matrix for sigma start NB: when forming design matrix must paste together responses in blocks old, moult, new
  int<lower=0,upper=1> lumped; //indicator variable for t2l likelihood
  int<lower=0,upper=1> llik; //indicator variable for calculating and returning pointwise log-likelihood
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
  real sigma_intercept = exp(beta_sigma[1]);
}

// The model to be estimated.
model {
  vector[N_old] P;
  vector[N_moult] Q;
  vector[N_new] R;
  vector[N_old+N_moult+N_new] mu;//start date lin pred
  vector[N_old+N_moult+N_new] tau;//duration lin pred
  vector[N_old+N_moult+N_new] sigma;//duration lin pred

  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);

 if (lumped == 0){
    for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
 } else {
   for (i in 1:N_old) P[i] = (1 - Phi((old_dates[i] - mu[i])/sigma[i])) + Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
 }
//print(P);
for (i in 1:N_moult) Q[i] = Phi((moult_dates[i] - mu[i + N_old])/sigma[i + N_old]) - Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old])/sigma[i + N_old]);
//print(Q);
if (lumped == 0){
      for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]);
    } else {
      for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]) + (1 - Phi((new_dates[i] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]));
    }
//print(R);
//print(sum(log(P)));
//print(sum(log(Q)));
//print(sum(log(R)));
target += sum(log(P))+sum(log(Q))+sum(log(R));
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
sigma_mu_ind ~ normal(0,10);}

generated quantities{
//real end_date;
// end_date = mu + tau;
if (llik == 1){
  vector[N_old+N_moult+N_new] log_lik;
  vector[N_old] P;
  vector[N_moult] Q;
  vector[N_new] R;
  vector[N_old+N_moult+N_new] mu;//start date lin pred
  vector[N_old+N_moult+N_new] tau;//duration lin pred
  vector[N_old+N_moult+N_new] sigma;//duration lin pred

    mu = X_mu * beta_mu;
  //  print(mu);
    tau = X_tau * beta_tau;
  //  print(tau);
    sigma = exp(X_sigma * beta_sigma);

  for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
  //print(P);
  for (i in 1:N_moult) Q[i] = Phi((moult_dates[i] - mu[i + N_old])/sigma[i + N_old]) - Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old])/sigma[i + N_old]);
  //print(Q);
  for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]);

  log_lik = append_row(log(P), append_row(log(Q), log(R)));
  }
}
