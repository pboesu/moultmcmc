//
// This Stan program defines the Underhill-Zucchini Type 3 repeated measures model
// covariates need to be on the individual level, not the observation model
//
data {
  int<lower=0> N_ind;//number of recaptured individuals
  int<lower=0> N_moult;//J
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0>individual[N_moult];//individual identifier
  int<lower=0>individual_first_index[N_ind];//row first occurrence of each individual in the model frame
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
  vector[N_ind] mu_ind;//individual effect on sigma start date
}

transformed parameters{
  real sigma_intercept = exp(beta_sigma[1]);
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
for (i in 1:N_moult) q[i] = log(tau[i]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i]) | mu[i] + mu_ind[individual[i]], sigma[i]);//N.B. unlike P and R this returns a log density
for (i in 1:N_moult) Q[i] = Phi((moult_dates[i] - (mu[i]+ mu_ind[individual[i]]))/sigma[i]) - Phi((moult_dates[i] - tau[i] - (mu[i]+ mu_ind[individual[i]]))/sigma[i]);
//individual effects are drawn from the population distribution of start dates
mu_ind ~ normal(0, sigma[individual_first_index]);//short cut for now - sampling the first sigma for each individual - input function needs to check these are representative for each individual!
//print(R);
//print(sum(log(P)));
//print(sum(log(Q)));
//print(sum(log(R)));
target += sum(q);// - log(Q));//TODO: Why does this differ from the likelihood for the non-recaptures model?!
//priors
beta_mu[1] ~ uniform(0,366);
beta_tau[1] ~ uniform(0,366);
beta_sigma[1] ~ normal(0,5);
}

generated quantities{

//NB: code duplication for the likelihood calculation is less than ideal - refactor to a use stan function?
//real end_date;
// end_date = mu + tau;

vector[N_moult] log_lik;
vector[N_moult] q;
vector[N_moult] Q;
vector[N_moult] mu;//start date lin pred
vector[N_moult] tau;//duration lin pred
vector[N_moult] sigma;//duration lin pred

  mu = X_mu * beta_mu;
//  print(mu);
  tau = X_tau * beta_tau;
//  print(tau);
  sigma = exp(X_sigma * beta_sigma);

for (i in 1:N_moult) q[i] = log(tau[i]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i]) | mu[i], sigma[i]);//N.B. unlike Q this returns a log density
for (i in 1:N_moult) Q[i] = Phi((moult_dates[i] - mu[i])/sigma[i]) - Phi((moult_dates[i] - tau[i] - mu[i])/sigma[i]);

log_lik = (q - log(Q));
}
