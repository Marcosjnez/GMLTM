data {
  int N_subj; // Número de sujetos
  int N_item; // Número de ítems
  array[N_subj * N_item] int ID;
  array[N_subj * N_item] int item;
  array[N_subj * N_item] int ones;
  int M; // Número de componentes
  int K; // Número de reglas
  int n_eta;
  array[n_eta, 2] int indexes_eta;
  matrix[N_item, K] Q;
  matrix[N_item, M] C;
  array[N_subj * N_item] int y;
}
transformed data {
  int N = N_subj * N_item;
  matrix[N, M] Cx = C[item, ];
}
parameters {
  matrix<lower=0>[1, M] alpha;
  vector[n_eta] eta;
  matrix[N_subj, M] theta;
}
transformed parameters {
}
model {
  vector[N] p;
  matrix[N, M] mu;
  matrix[N_item, M] beta;
  matrix[K, M] eta_matrix;

  eta_matrix = rep_matrix(0, K, M);
  for(i in 1:n_eta) {
    eta_matrix[indexes_eta[i, 1], indexes_eta[i, 2]] = eta[i];
  }
  beta = Q * eta_matrix;

  mu = inv_logit(alpha[ones, ] .* (theta[ID, ] - beta[item, ])) .^ Cx;
  for(i in 1:N) {
    p[i] = prod(mu[i, ]);
  }
  target += normal_lpdf(to_vector(theta) | 0, 1);
  target += normal_lpdf(to_vector(alpha) | 0, 1) -
            normal_lccdf(0 | 0, 1);;
  target += normal_lpdf(eta | 0, 1);
  target += bernoulli_lpmf(y | p);
}
generated quantities {
  // vector[N] loglik;
  // for(i in 1:N) {
  //   loglik[i] = bernoulli_lpmf(y[i] | p[i]);
  // }
}

