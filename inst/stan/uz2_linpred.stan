//
// This Stan program defines the Underhill-Zucchini Type 2 model
//

// The input data is a vector 'y' of length 'N'.
data {
  //responses
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0> N_new;//K
  vector[N_new] new_dates;//v_k
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
}

transformed parameters{
  real sigma_intercept = exp(beta_sigma[1]);
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
    for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
 } else {
   for (i in 1:N_old) P[i] = (1 - Phi((old_dates[i] - mu[i])/sigma[i])) + Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
 }
//print(P);
  for (i in 1:N_moult) q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
  if (lumped == 0){
      for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]);
    } else {
      for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]) + (1 - Phi((new_dates[i] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]));
    }
//print(R);
//print(sum(log(P)));
//print(sum(log(Q)));
//print(sum(log(R)));
target += sum(log(P))+sum(q)+sum(log(R));
//priors
beta_mu[1] ~ uniform(0,366);
beta_tau[1] ~ uniform(0,366);
beta_sigma[1] ~ normal(0,2);
}

generated quantities{

//NB: code duplication for the likelihood calculation is less than ideal - refactor to a use stan function?
//real end_date;
// end_date = mu + tau;
if (llik == 1){
  vector[N_old+N_moult+N_new] log_lik;
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
  sigma = exp(X_sigma * beta_sigma);

  for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
    //print(P);
    //for (i in 1:N_moult) Q[i] = Phi((moult_dates[i] - mu[i + N_old])/sigma[i + N_old]) - Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old])/sigma[i + N_old]);
  for (i in 1:N_moult) q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
    //print(Q);
  for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]);

  log_lik = append_row(log(P), append_row(q, log(R)));
}
}
