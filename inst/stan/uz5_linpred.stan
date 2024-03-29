//
// This Stan program defines the Underhill-Zucchini Type 5 model
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0,upper=1> flat_prior;//translated logical use flat prior on start and duration?
  //responses
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
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
}

transformed parameters{
  real sigma_intercept = exp(beta_sigma[1]);
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
  P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
  Rt[i] = Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
}
for (i in 1:N_moult){
   Ru[i] = Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old])/sigma[i + N_old]);
   q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
}

//print(sum(q));
//print(sum(log(P)));
//print(sum(log(1-Rt)));
//print(sum(log(1-Ru)));
target += sum(log(P) - log1m(Rt)) + sum(q - log1m(Ru)) ;
//priors
if (flat_prior == 1) {
 beta_mu[1] ~ uniform(0,366);
 beta_tau[1] ~ uniform(0,366);
} else {
 beta_mu[1] ~ normal(150,50)T[0,366];
 beta_tau[1] ~ normal(100,30)T[0,366];
}
beta_sigma[1] ~ normal(0,5);
}

generated quantities{

//NB: code duplication for the likelihood calculation is less than ideal - refactor to a use stan function?
//real end_date;
// end_date = mu + tau;

vector[N_old + N_moult] log_lik;
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
  sigma = exp(X_sigma * beta_sigma);

for (i in 1:N_old) {
  P[i] = 1 - Phi((old_dates[i] - mu[i])/sigma[i]);
  Rt[i] = Phi((old_dates[i] - tau[i] - mu[i])/sigma[i]);
}
for (i in 1:N_moult){
   Ru[i] = Phi((moult_dates[i] - tau[i + N_old] - mu[i + N_old])/sigma[i + N_old]);
   q[i] = log(tau[i + N_old]) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau[i + N_old]) | mu[i + N_old], sigma[i + N_old]);//N.B. unlike P and R this returns a log density
}
log_lik = append_row((log(P) - log1m(Rt)), (q - log1m(Ru)));
}
