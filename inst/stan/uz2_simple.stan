//
// This Stan program defines the Underhill-Zucchini Type 1 model
//

data {
  int<lower=0> N_old;//I in original derivation
  vector[N_old] old_dates;//t_i
  int<lower=0> N_moult;//J
  vector[N_moult] moult_dates;//u_j
  vector<lower=0,upper=1>[N_moult] moult_indices;//index of moult
  int<lower=0> N_new;//K
  vector[N_new] new_dates;//v_k
}

parameters {
  real mu;//mean start date
  real<lower=0> tau;//moult duration
  real<lower=0> sigma;//population varianmoult start date
}

model {
  vector[N_old] P;
  vector[N_moult] q;
  vector[N_new] R;

for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu)/sigma);
for (i in 1:N_moult) q[i] = log(tau) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau) | mu, sigma);//N.B. unlike P and R this returns a log density
for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau - mu)/sigma);

target += sum(log(P))+sum(q)+sum(log(R));
//priors
//mu ~ normal(0,10);
//tau ~ normal(0,10);
//sigma ~ normal(0,3);
}

generated quantities{
real end_date;
vector[N_old+N_moult+N_new] log_lik;
vector[N_old] P;
vector[N_moult] q;
vector[N_new] R;



for (i in 1:N_old) P[i] = 1 - Phi((old_dates[i] - mu)/sigma);
for (i in 1:N_moult) q[i] = log(tau) + normal_lpdf((moult_dates[i] - moult_indices[i]*tau) | mu, sigma);//N.B. unlike P and R this returns a log density
for (i in 1:N_new) R[i] = Phi((new_dates[i] - tau - mu)/sigma);

log_lik = append_row(log(P), append_row(q, log(R)));

end_date = mu + tau;
}
