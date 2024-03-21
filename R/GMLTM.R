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
#' GMLTM(data, Q, components, iters = 2000, iter_warmup = 1000, chains = 2, seed = 1606)
#'
#' @param data $$n \times p$$ data.frame or data matrix with the individuals in rows and items in columns.
#' @param Q $$p \times q$$ Q matrix.
#' @param components List of $$d$$ elements relating each components to a vector of rules.
#' @param iters Number of samples from the posterior distribution.
#' @param iter_warmup Number of samples to discard then initializing a Markov chain.
#' @param chains Number of Markov chains.
#' @param seed Random seed.
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
                  iter_warmup = 1000, seed = 1606) {

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

  data_list <- list(N_subj = N_subj, N_item = N_item,
                    ID = ID, item = item, K = K, M = M,
                    indexes_eta = indexes_eta, n_eta = nrow(indexes_eta),
                    indexes_alpha = indexes_alpha, n_alpha = nrow(indexes_alpha),
                    Q = Q, C = C, y = y, cf = 1)

  GMLTM_model <- cmdstanr::cmdstan_model("R/models/GMLTM.stan", compile = TRUE,
                                         cpp_options = list(stan_threads = TRUE))

  fit <- GMLTM_model$sample(data = data_list, chains = chains,
                              iter_sampling = iters, iter_warmup = iter_warmup,
                              parallel_chains = chains, validate_csv = FALSE,
                              threads_per_chain = 1, seed = seed,
                              output_dir = "posteriors/", output_basename = "GMLTM")

  theta <- fit$draws("theta", format = "draws_matrix")
  alpha <- fit$draws("alpha", format = "draws_matrix")
  eta <- fit$draws("eta", format = "draws_matrix")
  beta <- fit$draws("beta", format = "draws_matrix")
  guessing <- fit$draws("c", format = "draws_matrix")

  # vector_names <- as.character(fread("posteriors/GMLTM-1.csv", header = FALSE,
  #                                    skip = 45, fill = TRUE, nrows = 1))
  # x1 <- fread("posteriors/GMLTM-1.csv", header = FALSE, skip = 50, fill = TRUE,
  #             col.names = vector_names, nrows = iters)
  # x2 <- fread("posteriors/GMLTM-2.csv", header = FALSE, skip = 50, fill = TRUE,
  #             col.names = vector_names, nrows = iters)
  # x <- do.call(rbind, args = list(x1, x2))

  result <- list(theta = theta, alpha = alpha, eta = eta, beta = beta, guessing = guessing)
  return(result)

}

















