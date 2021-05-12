//
// This Stan program defines the combined Underhill-Zucchini Type 1 and 2 model
//

// The input data is a vector 'y' of length 'N'.
data {
  //responses
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0> N_moult_cat;//number of moult observations where no score is available
  vector[N_moult_cat] moult_cat_dates;//dates of moult observations without score
  int<lower=0> N_new;//K
  vector[N_new] new_dates;//v_k
  //predictors
  int N_pred_mu;//number of predictors for start date
  matrix[N_old+N_moult+N_moult_cat+N_new,N_pred_mu] X_mu;//design matrix for start date NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_tau;//number of predictors for duration
  matrix[N_old+N_moult+N_moult_cat+N_new,N_pred_tau] X_tau;//design matrix for duration NB: when forming design matrix must paste together responses in blocks old, moult, new
  int N_pred_sigma;//number of predictors for start date sigma
  matrix[N_old+N_moult+N_moult_cat+N_new,N_pred_sigma] X_sigma;//design matrix for sigma start NB: when forming design matrix must paste together responses in blocks old, moult, new
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
  vector[N_moult_cat] Q;
  vector[N_new] R;
  vector[N_old+N_moult+N_moult_cat+N_new] mu;//start date lin pred
  vector[N_old+N_moult+N_moult_cat+N_new] tau;//duration lin pred
  vector[N_old+N_moult+N_moult_cat+N_new] sigma;//duration lin pred

  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);//use log link for variance lin pred

for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
//print(P);
for (i in 1:N_moult) q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
for (i in 1:N_moult_cat) Q[i] = Phi((moult_cat_dates[i] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]) - Phi((moult_cat_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]);
for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult + N_moult_cat] - mu[i + N_old + N_moult + N_moult_cat])/sigma[i + N_old + N_moult + N_moult_cat]);
//print(R);
//print(sum(log(P)));
//print(sum(log(Q)));
//print(sum(log(R)));
target += sum(log(P))+sum(q)+sum(log(Q))+sum(log(R));
//priors
//mu ~ normal(0,10);
//tau ~ normal(0,10);
//sigma ~ normal(0,3);
}

generated quantities{

//NB: code duplication for the likelihood calculation is less than ideal - refactor to a use stan function?
//real end_date;
// end_date = mu + tau;

  vector[N_old+N_moult+N_moult_cat+N_new] log_lik;
  vector[N_old] P;
  vector[N_moult] q;
  vector[N_moult_cat] Q;
  vector[N_new] R;
  vector[N_old+N_moult+N_moult_cat+N_new] mu;//start date lin pred
  vector[N_old+N_moult+N_moult_cat+N_new] tau;//duration lin pred
  vector[N_old+N_moult+N_moult_cat+N_new] sigma;//duration lin pred

  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);

for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
//print(P);
for (i in 1:N_moult) q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
for (i in 1:N_moult_cat) Q[i] = Phi((moult_cat_dates[i] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]) - Phi((moult_cat_dates[i] - tau[i + N_old + N_moult] - mu[i + N_old + N_moult])/sigma[i + N_old + N_moult]);
for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau[i + N_old + N_moult + N_moult_cat] - mu[i + N_old + N_moult + N_moult_cat])/sigma[i + N_old + N_moult + N_moult_cat]);

log_lik = append_row(log(P), append_row(q, append_row(log(Q), log(R))));
}
