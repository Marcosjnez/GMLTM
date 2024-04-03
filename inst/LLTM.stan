
data {
  int N_subj; // Número de sujetos
  int N_item; // Número de ítems
  array[N_subj*N_item] int ID;
  array[N_subj*N_item] int item;
  int K; // Número de reglas
  matrix[N_item, K] Q;
  array[N_subj*N_item] int y;
}
transformed data {
  int N = N_subj*N_item;
}
parameters {
  vector[N_subj] theta;
  vector[K] eta;
}
transformed parameters {
  vector[N_item] beta;
  beta = Q * eta;
}
model {
  vector[N] p;
  p = inv_logit(theta[ID] - beta[item]);
  target += normal_lpdf(theta | 0, 1);
  target += normal_lpdf(eta | 0, 1);
  target += bernoulli_lpmf(y | p);
}
generated quantities {
  // vector[N] log_lik;
  // for(i in 1:N) {
  //     log_lik[i] = bernoulli_lpmf(y[i] | p[i]);
  // }
}

