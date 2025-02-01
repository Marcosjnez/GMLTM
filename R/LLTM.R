#' @title
#' The Logistic Latent Trait Model
#' @description
#'
#' Estimate the parameters of the LLTM
#'
#' @usage
#'
#' LLTM(data, Q, components, iters = 2000, iter_warmup = 1000, chains = 2, quantiles = c(0.025, 0.50,  0.975))
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
#' @details \code{LLTM} Estimates the Bayesian version of the linear logistic model LLTM (Fisher, 1973), a model that extends the unidimensional Rasch model and is suitable for the analysis of items comprising tasks solved by the examinee through the combination of a certain number of cognitive operations or rules.
#'
#' @return \code{cmdstan} object:
#'
#' @references
#'
#' Fischer, G. H. (1973). The linear logistic test model as an instrument in educational research. \emph{Acta psychologica}, 37(6), 359-374.\cr
#'
#' Lozano, J. H., & Revuelta, J. (2021). Bayesian estimation and testing of a linear logistic test model for learning during the test. Applied Measurement in Education, 34(3), 223-235.\cr
#'
#' Ramírez, E.S.; Jiménez, M.; Franco, V.R.; Alvarado, J.M. Delving into the Complexity of Analogical Reasoning: A Detailed Exploration with the Generalized Multicomponent Latent Trait Model for Diagnosis. \emph{J. Intell.} 2024, 12, 67. https://doi.org/10.3390/jintelligence12070067
#'
#'@examples
#' \dontrun{
#' # Load the data
#' data(analogy)
#' # These data correspond to the study:
#' # Blum, D., Holling, H., Galibert, M. S., & Forthmann, B. (2016). Task difficulty prediction of figural analogies. Intelligence, 56, 72-81.
#' # Define the dataset
#'
#' dataset <- analogy
#'
#' # Define the Q matrix with item rules
#' Q <- structure(c(0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1,
#'                  1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
#'                  0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0,
#'                  0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1,
#'                  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
#'                  0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1,
#'                  1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1), dim = c(27L, 5L), dimnames = list(
#'                   c("item 1", "item 2", "item 3", "item 4", "item 5", "item 6",
#'                     "item 7", "item 8", "item 9", "item 10", "item 11", "item 12",
#'                     "item 13", "item 14", "item 15", "item 16", "item 17", "item 18",
#'                     "item 19", "item 20", "item 21", "item 22", "item 23", "item 24",
#'                     "item 25", "item 26", "item 27"), c("rot_fig", "rot_trap",
#'                                                         "reflection", "subt_seg", "mov_point"))) #The labels correspond to the names of the rules.
#'
#'
#' # Define the fit3 object
#' fit3 <- LLTM(dataset, Q, iters = 2000, iter_warmup = 500,
#'              quantiles = c(0.05, 0.5, 0.95), chains = 2, parallel_chains = 2)
#'
#' # Print the difficulty parameter for each rule.
#' fit2$EAP$eta
#'
#' # Print the difficulty parameter for each item. For items with rules from both components, print both difficulties.
#' fit3$EAP$beta
#'
#' # Provides a summary of parameter estimates and diagnoses convergence in a Bayesian model for the eta parameter.
#' fit3$fit$print("eta")
#'
#' # Provides a summary of parameter estimates and diagnoses convergence in a Bayesian model for the beta parameter.
#' fit3$fit$print("beta")
#'
#' # Marginal reliability for each component
#' reliability(fit3)
#'
#' # Print all ppchecks of the model
#' check3 <- ppchecks(fit3)
#' }
#'
#' @export
LLTM <- function(data, Q, components, iters = 2000, chains = 2,
                 iter_warmup = 1000, quantiles = c(0.025, 0.975),
                 parallel_chains = 1, threads_per_chain = 1,
                 ...) {

  K <- ncol(Q) # Number of rules + interactions
  y <- unname(unlist(dataset)) # Vector of data
  N_subj <- nrow(dataset) # Number of subjects
  N_item <- ncol(dataset) # Number of items
  ID <- rep(1:N_subj, times = N_item) # Subject identification for each data point
  item <- rep(1:N_item, each = N_subj) # Item identification for each data point

  data_list <- list(N_subj = N_subj, N_item = N_item,
                    ID = ID, item = item, K = K,
                    Q = Q, y = y)

  stan_file <- system.file("LLTM.stan", package = "GMLTM")
  model <- cmdstanr::cmdstan_model(stan_file, compile = TRUE,
                                   cpp_options = list(stan_threads = TRUE))

  fit <- model$sample(data = data_list, chains = chains,
                      iter_sampling = iters, iter_warmup = iter_warmup,
                      parallel_chains = parallel_chains,
                      threads_per_chain = threads_per_chain,
                      validate_csv = FALSE, ...)

  theta <- fit$draws("theta", format = "draws_matrix")
  eta <- fit$draws("eta", format = "draws_matrix")
  beta <- fit$draws("beta", format = "draws_matrix")

  # Get EAP estimates:
  summary_theta <- colMeans(theta)
  summary_eta <- colMeans(eta)
  summary_beta <- colMeans(beta)

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
  J <- length(quantiles)

  # thetas:
  colnames(theta) <- subject_names
  theta_EAP <- summary_theta
  names(theta_EAP) <- subject_names
  quantiles_theta <- t(apply(theta, MARGIN = 2, quantile, probs = quantiles))
  rownames(quantiles_theta) <- subject_names
  colnames(quantiles_theta) <- as.character(quantiles)

  # etas and betas:
  colnames(eta) <- rule_names
  colnames(beta) <- item_names
  eta_EAP <- summary_eta
  names(eta_EAP) <- rule_names
  beta_EAP <- Q %*% eta_EAP
  quantiles_eta <- t(apply(eta, MARGIN = 2, quantile, probs = quantiles))
  quantiles_beta <- t(apply(beta, MARGIN = 2, quantile, probs = quantiles))
  rownames(quantiles_eta) <- rule_names
  colnames(quantiles_eta) <- as.character(quantiles)
  rownames(quantiles_beta) <- item_names
  colnames(quantiles_beta) <- as.character(quantiles)

  # probabilities and loglik:
  iters <- iters*chains
  y <- unname(unlist(data)) # Vector of data
  p <- loglik <- array(NA, dim = c(iters, N_subj, N_item))
  for(i in 1:iters) {
    p[i, , ] <- plogis(theta[i, ID] - beta[i, item])
    loglik[i, , ] <- dbinom(y, size = 1, prob = c(p[i, , ]), log = TRUE)
  }

  EAP <- list(theta = theta_EAP, eta = eta_EAP, beta = beta_EAP)
  quantiles <- list(theta = quantiles_theta, eta = quantiles_eta, beta = quantiles_beta)
  posterior <- list(theta = theta, eta = eta, beta = beta,
                    probabilities = p, loglik = loglik)

  result <- list(EAP = EAP, quantiles = quantiles, posterior = posterior,
                 model = model, fit = fit, data = data)
  class(result) <- "LLTM"
  return(result)

}












