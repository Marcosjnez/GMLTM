#' @title
#' The Multidimensional Latent Trait Model
#' @description
#'
#' Estimate the parameters of the MLTM.
#'
#' @usage
#'
#' GMLTM(data, Q, components, iters = 2000, iter_warmup = 1000, chains = 2, quantiles = c(0.025, 0.50,  0.975))
#'
#' @param data $$n \times p$$ data.frame or data matrix with the individuals in rows and items in columns.
#' @param Q $$p \times q$$ Q matrix.
#' @param components List of $$d$$ elements relating each components to a vector of rules.
#' @param iters Number of samples from the posterior distribution.
#' @param iter_warmup Number of samples to discard then initializing a Markov chain.
#' @param chains Number of Markov chains.
#' @param quantiles Vector of probabilities associated with the posterior quantiles.
#' @param parallel_chains Number of chains to run in parallel cores.
#' @param threads_per_chain Number of cores to run within a chain.
#' @param ... Additional arguments to pass to the \code{sample} function from \code{cmdstanr}.
#'
#' @details \code{MLTM} estimates a ...
#'
#' @return \code{cmdstan} object:
#'
#' @references
#'
#' Ramirez E., Jiménez M., Franco V., Alvarado J. (2024). Delving into... preprint.
#'
#' @export
MLTM <- function(data, Q, components, iters = 2000, chains = 2,
                  iter_warmup = 1000, quantiles = c(0.025, 0.975),
                  parallel_chains = 1, threads_per_chain = 1,
                  ...) {

  eta <- get_eta(components) # 1 indicates the parameter is estimated and 0 otherwise
  C <- get_C(Q, components) # Matrix of components
  binary_eta <- get_eta(components) # Binary matrix of rules into components
  indexes_eta <- locate_eta(binary_eta) # Indexes indicating the position of the parameters in the eta matrix
  M <- length(components) # Number of components
  K <- ncol(Q) # Number of rules + interactions
  y <- unname(unlist(dataset)) # Vector of data
  N_subj <- nrow(dataset) # Number of subjects
  N_item <- ncol(dataset) # Number of items
  ID <- rep(1:N_subj, times = N_item) # Subject identification for each data point
  item <- rep(1:N_item, each = N_subj) # Item identification for each data point
  n_eta <- nrow(indexes_eta)
  ones <- rep(1, N_subj*N_item)

  data_list <- list(N_subj = N_subj, N_item = N_item,
                    ID = ID, item = item, K = K, M = M,
                    indexes_eta = indexes_eta, n_eta = n_eta,
                    Q = Q, C = C, y = y, ones = ones)

  stan_file <- system.file("MLTM.stan", package = "GMLTM")
  model <- cmdstanr::cmdstan_model(stan_file, compile = TRUE,
                                   cpp_options = list(stan_threads = TRUE))

  fit <- model$sample(data = data_list, chains = chains,
                      iter_sampling = iters, iter_warmup = iter_warmup,
                      parallel_chains = parallel_chains,
                      threads_per_chain = threads_per_chain,
                      validate_csv = FALSE, ...)

  theta <- fit$draws("theta", format = "draws_matrix")
  alpha <- fit$draws("alpha", format = "draws_matrix")
  eta <- fit$draws("eta", format = "draws_matrix")

  # Get EAP estimates:
  summary_theta <- colMeans(theta)
  summary_alpha <- colMeans(alpha)
  summary_eta <- colMeans(eta)

  if(is.null(colnames(data))) {
    item_names <- paste("item", 1:ncol(data), sep = "")
  } else {
    item_names <- colnames(data)
  }
  if(is.null(rownames(data))) {
    subject_names <- 1:nrow(data)
  } else {
    subject_names <- rownames(data)
  }
  if(is.null(colnames(Q))) {
    rule_names <- paste("rule", 1:K, sep = "")
  } else {
    rule_names <- colnames(Q)
  }
  rownames(Q) <- item_names
  comp_names <- paste("Component", 1:M, sep = "")
  J <- length(quantiles)

  # thetas:
  theta_EAP <- matrix(summary_theta, nrow = N_subj, ncol = M)
  rownames(theta_EAP) <- subject_names
  colnames(theta_EAP) <- comp_names
  quant_theta <- t(apply(theta, MARGIN = 2, quantile, probs = quantiles))
  quantiles_theta <- list()
  for(j in 1:J) {
    temp <- matrix(quant_theta[, j], nrow = N_subj, ncol = M)
    rownames(temp) <- subject_names
    colnames(temp) <- comp_names
    quantiles_theta[[j]] <- temp
  }
  names(quantiles_theta) <- as.character(quantiles)

  # alphas:
  colnames(alpha) <- comp_names
  alpha_EAP <- c(summary_alpha)
  names(alpha_EAP) <- comp_names
  quantiles_alpha <- t(apply(alpha, MARGIN = 2, quantile, probs = quantiles))
  colnames(quantiles_alpha) <- as.character(quantiles)

  # etas and betas:
  eta_EAP <- matrix(0, K, M)
  for(i in 1:n_eta) {
    eta_EAP[indexes_eta[i, 1], indexes_eta[i, 2]] <- summary_eta[i];
  }
  rownames(eta_EAP) <- rule_names
  colnames(eta_EAP) <- comp_names
  beta_EAP <- Q %*% eta_EAP
  quant_eta <- t(apply(eta, MARGIN = 2, quantile, probs = quantiles))
  quantiles_eta <- quantiles_beta <- list()
  for(j in 1:J) {
    temp <- matrix(0, nrow = K, ncol = M)
    for(i in 1:n_eta) {
      temp[indexes_eta[i, 1], indexes_eta[i, 2]] <- quant_eta[i, j];
    }
    rownames(temp) <- rule_names
    colnames(temp) <- comp_names
    quantiles_eta[[j]] <- temp
    quantiles_beta[[j]] <- Q %*% temp
  }
  names(quantiles_eta) <- as.character(quantiles)
  names(quantiles_beta) <- as.character(quantiles)

  # probabilities and loglik:
  y <- unname(unlist(data)) # Vector of data
  p <- loglik <- array(NA, dim = c(iters, N_subj, N_item))
  theta_matrix <- array(theta, dim = c(iters, N_subj, M))
  eta_matrix <- array(0, dim = c(iters, K, M))
  beta_matrix <- array(0, dim = c(iters, N_item, M))
  comp <- rep(1, N_subj*N_item)
  for(i in 1:iters) {
    for(j in 1:n_eta) {
      eta_matrix[i, indexes_eta[j, 1], indexes_eta[j, 2]] <- eta[i, j];
    }
    beta_matrix[i, , ] <- Q %*% eta_matrix[i, , ]
    mu <- plogis(alpha[i, ][comp, ] * (theta_matrix[i, ID, ] - beta_matrix[i, item, ])) ^ C[item, ];
    p[i, , ] <- apply(mu, MARGIN = 1, prod)
    loglik[i, , ] <- dbinom(y, size = 1, prob = c(p[i, , ]), log = TRUE)
  }

  EAP <- list(theta = theta_EAP, alpha = alpha_EAP, eta = eta_EAP, beta = beta_EAP)
  quantiles <- list(theta = quantiles_theta, alpha = quantiles_alpha, eta = quantiles_eta,
                    beta = quantiles_beta)
  posterior <- list(theta = theta_matrix, alpha = alpha, eta = eta_matrix, beta = beta_matrix,
                    probabilities = p, loglik = loglik)

  result <- list(EAP = EAP, quantiles = quantiles, posterior = posterior,
                 model = model, fit = fit, data = data)
  class(result) <- "MLTM"
  return(result)

}












