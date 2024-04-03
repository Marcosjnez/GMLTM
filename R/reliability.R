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
