get_eta <- function(components) {

  M <- length(components)
  K <- length(unique(unlist(components)))
  eta <- matrix(0, K, M)
  for(m in 1:M) {
    eta[components[[m]], m] <- 1
  }

  return(eta)

}
get_C <- function(Q, components) { # Which items belong to each component

  M <- length(components)
  C <- matrix(0, nrow(Q), M)

  for(i in 1:M) {
    sums1 <- rowSums(Q[, components[[i]]])
    indexes <- which(sums1 > 0)
    C[indexes, i] <- 1
  }

  return(C)

}
locate_eta <- function(eta) { # Which rules belong to each component

  indexes <- which(ifelse(eta == 1, TRUE, FALSE), arr.ind = TRUE)

  return(indexes)

}
locate_alpha <- function(Q, components) {

  C <- get_C(Q, components)
  indexes_2 <- locate_eta(C)
  M <- length(components)
  p <- nrow(C)
  i <- 1
  x <- c()

  for(m in 1:M) {

    X <- Q[, components[[m]]]
    D <- cbind(C[, m], 0)
    Z <- X %*% t(X) - rowSums(X)
    Z2 <- (t(Z) + Z) / 2
    nr <- 1:nrow(X)

    repeat{
      j <- min(nr)
      indexes <- which(Z2[j, ] == 0)
      nr <- nr[-match(indexes, nr)]
      D[indexes, 2] <- i
      if(length(nr) == 0) break
      i <- i+1
    }

    x <- c(x, D[D[, 1] != 0, 2])

  }

  x <- as.numeric(as.factor(x))
  indexes_2 <- cbind(indexes_2, x)

  return(indexes_2)

}

#' @title
#' The Generalized Multidimensional Latent Trait Model
#' @description
#'
#' Estimate the parameters of the GMLTM.
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
#' @details \code{GMLTM} estimates a ...
#'
#' @return \code{cmdstan} object:
#'
#' @references
#'
#' Ramirez E., Jiménez M., Franco V., Alvarado J. (2024). Delving into... preprint.
#'
#' @export
GMLTM <- function(data, Q, components, iters = 2000, chains = 2,
                  iter_warmup = 1000, quantiles = c(0.025, 0.975),
                  parallel_chains = 1, threads_per_chain = 1,
                  ...) {

  eta <- get_eta(components) # 1 indicates the parameter is estimated and 0 otherwise
  C <- get_C(Q, components) # Matrix of components
  binary_eta <- get_eta(components) # Binary matrix of rules into components
  indexes_eta <- locate_eta(binary_eta) # Indexes indicating the position of the parameters in the eta matrix
  indexes_alpha <- locate_alpha(Q, components) # Indexes indicating the position of the alpha parameters in the alpha matrix
  M <- length(components) # Number of components
  K <- ncol(Q) # Number of rules + interactions
  y <- unname(unlist(dataset)) # Vector of data
  N_subj <- nrow(dataset) # Number of subjects
  N_item <- ncol(dataset) # Number of items
  ID <- rep(1:N_subj, times = N_item) # Subject identification for each data point
  item <- rep(1:N_item, each = N_subj) # Item identification for each data point
  n_alpha <- nrow(indexes_alpha)
  n_eta <- nrow(indexes_eta)

  data_list <- list(N_subj = N_subj, N_item = N_item,
                    ID = ID, item = item, K = K, M = M,
                    indexes_eta = indexes_eta, n_eta = n_eta,
                    indexes_alpha = indexes_alpha, n_alpha = n_alpha,
                    Q = Q, C = C, y = y, cf = 1)

  stan_file <- system.file("GMLTM.stan", package = "GMLTM")
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
  guessing <- fit$draws("c", format = "draws_matrix")

  # Get EAP estimates:
  summary_theta <- colMeans(theta)
  summary_alpha <- colMeans(alpha)
  summary_eta <- colMeans(eta)
  summary_guessing <- colMeans(guessing)

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
  alpha_EAP <- matrix(0, N_item, M)
  for(i in 1:n_alpha) {
    alpha_EAP[indexes_alpha[i, 1], indexes_alpha[i, 2]] <- summary_alpha[indexes_alpha[i, 3]];
  }
  rownames(alpha_EAP) <- item_names
  colnames(alpha_EAP) <- comp_names
  quant_alpha <- t(apply(alpha, MARGIN = 2, quantile, probs = quantiles))
  quantiles_alpha <- list()
  for(j in 1:J) {
    temp <- matrix(0, nrow = N_item, ncol = M)
    for(i in 1:n_alpha) {
      temp[indexes_alpha[i, 1], indexes_alpha[i, 2]] <- quant_alpha[indexes_alpha[i, 3], j];
    }
    rownames(temp) <- item_names
    colnames(temp) <- comp_names
    quantiles_alpha[[j]] <- temp
  }
  names(quantiles_alpha) <- as.character(quantiles)

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

  # guessings:
  guessing_EAP <- summary_guessing
  names(guessing_EAP) <- item_names
  quant_guessing <- t(apply(guessing, MARGIN = 2, quantile, probs = quantiles))
  quantiles_guessing <- list()
  for(j in 1:J) {
    temp <- quant_guessing[, j]
    names(temp) <- item_names
    quantiles_guessing[[j]] <- temp
  }
  names(quantiles_guessing) <- as.character(quantiles)

  # probabilities and loglik:
  y <- unname(unlist(data)) # Vector of data
  p <- loglik <- array(NA, dim = c(iters, N_subj, N_item))
  theta_matrix <- array(theta, dim = c(iters, N_subj, M))
  alpha_matrix <- array(0, dim = c(iters, N_item, M))
  eta_matrix <- array(0, dim = c(iters, K, M))
  beta_matrix <- array(0, dim = c(iters, N_item, M))
  for(i in 1:iters) {
    for(j in 1:n_alpha) {
      alpha_matrix[i, indexes_alpha[j, 1], indexes_alpha[j, 2]] <- alpha[i, indexes_alpha[j, 3]];
    }
    for(j in 1:n_eta) {
      eta_matrix[i, indexes_eta[j, 1], indexes_eta[j, 2]] <- eta[i, j];
    }
    beta_matrix[i, , ] <- Q %*% eta_matrix[i, , ]
    mu <- plogis(alpha_matrix[i, item, ] * (theta_matrix[i, ID, ] - beta_matrix[i, item, ])) ^ C[item, ];
    p[i, , ] <- guessing[i, item] + (1 - guessing[i, item]) * apply(mu, MARGIN = 1, prod)
    loglik[i, , ] <- dbinom(y, size = 1, prob = c(p[i, , ]), log = TRUE)
  }

  EAP <- list(theta = theta_EAP, alpha = alpha_EAP, eta = eta_EAP,
              beta = beta_EAP, guessing = guessing_EAP)
  quantiles <- list(theta = quantiles_theta, alpha = quantiles_alpha, eta = quantiles_eta,
                    beta = quantiles_beta, guessing = quantiles_guessing)
  posterior <- list(theta = theta_matrix, alpha = alpha_matrix, eta = eta_matrix, beta = beta_matrix, guessing = guessing,
                    probabilities = p, loglik = loglik)

  result <- list(EAP = EAP, quantiles = quantiles, posterior = posterior,
                 model = model, fit = fit, data = data)
  class(result) <- "GMLTM"
  return(result)

}

#' @title
#' Marginal reliability
#' @description
#'
#' Estimate the the marginal reliability of the GMLTM.
#'
#' @usage
#'
#' reliability(fit)
#'
#' @param fit MLTM object.
#'
#' @details \code{reliability} estimates a ...
#'
#' @return A number denoting the reliability estimate.
#'
#' @references
#'
#' Ramirez E., Jiménez M., Franco V., Alvarado J. (2024). Delving into... preprint.
#'
#' @export
reliability <- function(fit) {

  EAP_theta <- fit$EAP$theta
  VAR_theta <- apply(fit$posterior$theta, MARGIN = c(2, 3), FUN = var)
  n <- nrow(EAP_theta)
  q <- ncol(EAP_theta)
  COV <- cov(EAP_theta)
  e <- colMeans(VAR_theta)
  s <- diag(COV)
  rxx <- s/(s+e)

  return(rxx)

}

#' @title
#' Posterior predictive checks, fitted values, and prediction intervals.
#' @description
#'
#' Plot the frequency distribution of the empirical and simulated total scores and obtain the fitted values and prediction intervals.
#'
#' @usage
#'
#' ppchecks(fit, nsim = 100, interval = 0.95)
#'
#' @param fit MLTM object.
#' @param nsim Number of simulated samples.
#' @param interval Probability associated with the credibility intervals.
#' @param ... Additional arguments to pass to the \code{plot} function.
#'
#' @details \code{ppchecks} estimates a ...
#'
#' @return An array of simulated scores.
#'
#' @references
#'
#' Ramirez E., Jiménez M., Franco V., Alvarado J. (2024). Delving into... preprint.
#'
#' @export
ppchecks <- function(fit, nsim = 100, interval = 0.95, ...) {

  n <- nrow(fit$data)
  p <- ncol(fit$data)
  scores <- rowSums(fit$data)
  ysim <- array(NA, dim = c(nsim, n, p))
  probs <- fit$posterior$probabilities
  for(i in 1:nsim) {
    ysim[i, , ] <- rbinom(n*p, size = 1, prob = c(probs[i, , ]))
  }
  scores_ysim <- apply(ysim, MARGIN = 2, rowSums)
  # max(table(scores)) > max(table(scores_ysim))

  p1 <- hist(scores, breaks = 0:p, plot = FALSE)
  plot(p1, col = rgb(0, 0, 1, 1/4), freq = FALSE,
       xlab = "Total scores", main = "", xlim = c(0, p), ...)
  p2 <- hist(c(scores_ysim), breaks = 0:p, plot = FALSE)
  plot(p2, col = rgb(1, 0, 0, 1/4), freq = FALSE, add = TRUE)

  p <- fit$posterior$probabilities
  EAP_p <- apply(p, MARGIN = c(2, 3), mean)
  p_obs <- colMeans(fit$data)
  p_pred <- colMeans(EAP_p)
  marginals <- apply(p, MARGIN = c(1, 3), mean)
  CI_p <- t(apply(marginals, MARGIN = 2, quantile,
                  probs = c((1-interval)/2, 1-(1-interval)/2)))
  fitted_items <- cbind(CI_p[, 1], p_pred, CI_p[, 2], p_obs)
  colnames(fitted_items) <- c("lower", "predicted", "upper", "observed")

  fitted_subjects <- cbind(rowMeans(EAP_p), rowMeans(fit$data))
  colnames(fitted_subjects) <- c("predicted", "observed")

  result <- list(items = fitted_items, subjects = fitted_subjects,
                 fitted = EAP_p, ysim = ysim)

  # cols2 <- c(rgb(1, 0, 0, 1/2), rgb(0, 0, 1, 1/2), rgb(0.36, 0, 0.64, 0.5))
  # legend(17, 0.10, legend = c("Predicted", "Observed", "Overlap"),
  #        col = c(cols2[1], cols2[2], cols2[3]), lwd = 2, cex = 1,
  #        title = "\n Posterior \n Predictive Check", text.width = 3.75) # title = "Group" box.lty = 0

  return(invisible(result))

}












