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
