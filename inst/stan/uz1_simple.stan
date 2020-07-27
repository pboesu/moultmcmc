//
// This Stan program defines the Underhill-Zucchini Type 1 model
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  vector[N_moult] moult_dates;//u_j
  int<lower=0> N_new;//K
  vector[N_new] new_dates;//v_k
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;//mean start date
  real<lower=0> tau;//moult duration
  real<lower=0> sigma;//population varianmoult start date
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N_old] P;
  vector[N_moult] Q;
  vector[N_new] R;

for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu)/sigma);
for (i in 1:N_moult) Q[i] = Phi((moult_dates[i] - mu)/sigma) - Phi((moult_dates[i] - tau - mu)/sigma);
for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau - mu)/sigma);

target += sum(log(P))+sum(log(Q))+sum(log(R));
//priors
//mu ~ normal(0,10);
//tau ~ normal(0,10);
//sigma ~ normal(0,3);
}

generated quantities{
real end_date;
 end_date = mu + tau;
}
