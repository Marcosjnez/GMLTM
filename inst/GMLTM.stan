
data {
  int N_subj;
  int N_item;
  array[N_subj * N_item] int ID;
  array[N_subj * N_item] int item;
  int M;
  int K;
  int n_eta;
  int n_alpha;
  array[n_eta, 2] int indexes_eta;
  array[n_alpha, 3] int indexes_alpha;
  matrix[N_item, K] Q;
  matrix[N_item, M] C;
  array[N_subj * N_item] int y;
  int cf;
}
transformed data {
  int max_n_alpha = max(indexes_alpha[, 3]);
  int N = N_subj * N_item;
  matrix[N, M] Cx = C[item, ];
}
parameters {
  matrix[N_subj, M] theta;
  vector[n_eta] eta;
  vector<lower=0>[max_n_alpha] alpha;
  vector<lower=0, upper=1>[N_item] c;
}
transformed parameters {
}
model {
  matrix[K, M] eta_matrix;
  matrix[N_item, M] alpha_matrix;
  matrix[N_item, M] beta_matrix;
  matrix[N, M] mu;
  vector[N] p;

  eta_matrix = rep_matrix(0, K, M);
  alpha_matrix = rep_matrix(0, N_item, M);
  for(i in 1:n_eta) {
    eta_matrix[indexes_eta[i, 1], indexes_eta[i, 2]] = eta[i];
  }
  for(i in 1:n_alpha) {
    alpha_matrix[indexes_alpha[i, 1], indexes_alpha[i, 2]] = alpha[indexes_alpha[i, 3]];
  }
  beta_matrix = Q * eta_matrix;
  mu = inv_logit(alpha_matrix[item, ] .* (theta[ID, ] - beta_matrix[item, ])) .^ Cx;
  for(i in 1:N) {
    p[i] = c[item[i]] + (1-c[item[i]]) .* prod(mu[i, ]);
  }
  target += normal_lpdf(to_vector(theta) | 0, 1);
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(alpha | 0, 1) -
            normal_lccdf(0 | 0, 1);
  target += beta_lpdf(c | 3, 20);
  target += bernoulli_lpmf(y | p);
}
generated quantities {
  // vector[N] log_lik;
  // for(i in 1:N) {
  //   log_lik[i] = bernoulli_lpmf(y[i] | p[i]);
  // }
}

