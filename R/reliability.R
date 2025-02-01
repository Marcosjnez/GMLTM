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
#' Ramírez, E.S.; Jiménez, M.; Franco, V.R.; Alvarado, J.M. Delving into the Complexity of Analogical Reasoning: A Detailed Exploration with the Generalized Multicomponent Latent Trait Model for Diagnosis. \emph{J. Intell.} 2024, 12, 67. https://doi.org/10.3390/jintelligence12070067
#'
#' @export
reliability <- function(fit) {

  if(class(fit) == "LLTM") {
    EAP_theta <- matrix(fit$EAP$theta)
    VAR_theta <- matrix(apply(fit$posterior$theta, MARGIN = 2, FUN = var))
  } else {
    EAP_theta <- fit$EAP$theta
    VAR_theta <- apply(fit$posterior$theta, MARGIN = c(2, 3), FUN = var)
  }
  n <- nrow(EAP_theta)
  q <- ncol(EAP_theta)
  COV <- cov(EAP_theta)
  e <- colMeans(VAR_theta)
  s <- diag(COV)
  rxx <- s/(s+e)

  return(rxx)

}
