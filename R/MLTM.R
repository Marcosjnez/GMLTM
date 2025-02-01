#' @title
#' The Multidimensional Latent Trait Model
#' @description
#'
#' Estimate the parameters of the MLTM-D
#'
#' @usage
#'
#' GMLTM(data, Q, components, iters = 2000, iter_warmup = 1000, chains = 2, quantiles = c(0.025, 0.50,  0.975))
#'
#' @param data \eqn{n \times p} data.frame or data matrix with the individuals in rows and items in columns.
#' @param Q \eqn{p \times q} Q matrix.
#' @param components List of \eqn{d} elements relating each components to a vector of rules.
#' @param iters Number of samples from the posterior distribution.
#' @param iter_warmup Number of samples to discard then initializing a Markov chain.
#' @param chains Number of Markov chains.
#' @param quantiles Vector of probabilities associated with the posterior quantiles.
#' @param parallel_chains Number of chains to run in parallel cores.
#' @param threads_per_chain Number of cores to run within a chain.
#' @param ... Additional arguments to pass to the \code{sample} function from \code{cmdstanr}.
#'
#' @details \code{MLTM} Estimates the Bayesian version of the Multicomponent Latent Trait Model for Diagnosis (MLTM-D) by Embretson & Yang (2013), a noncompensatory latent trait model. MLTM-D specifies a hierarchical relationship between components and attributes, allowing for diagnosis at two levels. This model is applicable to broad trait measures, such as achievement tests, where the component structure varies between items and the items are composed of different cognitive operations.
#'
#' @return \code{cmdstan} object:
#'
#' @references
#'
#' Embretson, S. E., & Yang, X. (2013). A multicomponent latent trait model for diagnosis. \emph{Psychometrika}, 78, 14-36.\cr
#'
#' Ramírez, E.S.; Jiménez, M.; Franco, V.R.; Alvarado, J.M. Delving into the Complexity of Analogical Reasoning: A Detailed Exploration with the Generalized Multicomponent Latent Trait Model for Diagnosis. \emph{J. Intell.} 2024, 12, 67. https://doi.org/10.3390/jintelligence12070067
#'
#'@examples
#'\dontrun{
#' # Load the data
#'data(analogy)
#' # These data correspond to the study:
#'# Blum, D., Holling, H., Galibert, M. S., & Forthmann, B. (2016). Task difficulty prediction of figural analogies. Intelligence, 56, 72-81.
#'# Define the dataset
#'
#'dataset <- analogy
#'
#'# Define the Q matrix with item rules
#'Q <- structure(c(0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1,
#'                 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
#'                 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0,
#'                 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1,
#'                 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
#'                 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1,
#'                 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1), dim = c(27L, 5L), dimnames = list(
#'                   c("item 1", "item 2", "item 3", "item 4", "item 5", "item 6",
#'                     "item 7", "item 8", "item 9", "item 10", "item 11", "item 12",
#'                     "item 13", "item 14", "item 15", "item 16", "item 17", "item 18",
#'                     "item 19", "item 20", "item 21", "item 22", "item 23", "item 24",
#'                     "item 25", "item 26", "item 27"), c("rot_fig", "rot_trap",
#'                                                         "reflection", "subt_seg", "mov_point"))) #The labels correspond to the names of the rules.
#'
#'
#'# Define the components list
#'components <- list(global = c(1, 2, 3), local = c(4, 5))
#'
#'# Define the fit2 object
#'fit2 <- MLTM(dataset, Q, components, iters = 2000, iter_warmup = 500,
#'              quantiles = c(0.05, 0.5, 0.95), chains = 2, parallel_chains = 2)
#'
#'# Print the difficulty parameters for each rule in each component
#'fit2$EAP$eta
#'
#'# Print the difficulty parameter for each item. For items with rules from both components, print both difficulties.
#'fit2$EAP$beta
#'
#'# Print the discriminations of each item. When an item is composed of rules from both components, provide a discrimination for each component.
#'fit2$EAP$alpha
#'
#'# Provides a summary of parameter estimates and diagnoses convergence in a Bayesian model for the eta parameter.
#'fit2$fit$print("eta")
#'
#'# Provides a summary of parameter estimates and diagnoses convergence in a Bayesian model for the alpha parameter.
#'fit2$fit$print("alpha")
#'
#'# Marginal reliability for each component
#'reliability(fit2)
#'
#'# Print all ppchecks of the model
#'check2 <- ppchecks(fit2)
#'}
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
  iters <- iters*chains
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












