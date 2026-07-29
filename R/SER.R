#' Single Effect Regression (SER)
#'
#' @param X Observed data - matrix
#' @param y Observed response - vector
#' @param sigma_2 Residual variance - scalar
#' @param priors List of prior parameters
#'
#' @return list(mu_1, sigma_1_2, BF), where:
#' mu_1 is the posterior effect mean vector
#' sigma_1_2 is the posterior effect variance vector
#' BF is the Bayesian factor vector
#' @importFrom stats dnorm
#' @export
#'
#' @examples
SER <- function(X, y, sigma_2, priors){
  n  = nrow(X)

  #Unpack priors
  list2env(lapply(priors, `[[`, 1), envir = environment())

  d   <- attr(X, "d")
  xty <- compute_Xty(X, y)

  # Stable reparameterization (multiply num/den by sigma_0_2):
  den <- sigma_2 + sigma_0_2 * d

  # Conditional posterior mean and variance
  mu_1      <- as.vector( (sigma_2 * mu_0 + sigma_0_2 * xty) / den )
  sigma_1_2 <- (sigma_2 * sigma_0_2) / den


  # Original parameterization------------------
  # #Calculate conditional posterior means
  # mu_1 = as.vector(
  #   (sigma_2/sigma_0_2)*mu_0 + compute_Xty(X,y)
  # )/(
  #   (sigma_2/sigma_0_2) + attr(X, "d"))
  #
  # #Calculate conditional posterior variances
  # sigma_1_2 = sigma_2 /(
  #   (sigma_2/sigma_0_2) + attr(X, "d"))

  #Calculate posterior inclusion probabilities
  log_BF = BF(X, y, sigma_2, sigma_0_2, mu_0, xty = xty, d = d) #log of Bayesian Factor

  lpo = log_BF + log(prior_inclusion_weights + sqrt(.Machine$double.eps)) #log of posterior odds
  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)

  alpha = w_weighted / weighted_sum_w

  #Calculate unconditional posterior effect moments
  b_hat = as.matrix(alpha* mu_1)
  b_2_hat = as.matrix(alpha*(sigma_1_2 + mu_1^2))
  xb_hat = compute_Xb(X, b_hat)

  #Calculate the effect's contribution to KL objective (third term)
  lbf_model = maxlpo + log(weighted_sum_w)
  log_lik = lbf_model + sum(dnorm(y,0,sqrt(sigma_2),log = TRUE))
  KL = -log_lik -
    ((n/2)*log(2*pi*sigma_2) +
    (1/(2*sigma_2))*SER_ERSS(X, y, b_hat, b_2_hat, xb_hat = xb_hat))


  output_list = list(
      alpha = I(list(alpha)),
      b_hat = I(list(drop(b_hat))),
      b_2_hat = I(list(drop(b_2_hat))),
      KL = KL
    )

  return(output_list)
}

#' Bayesian Factor - in logspace
#'
#' @param X Design matrix (n x p) with standardization attributes.
#' @param y Response vector (length n).
#' @param sigma_2 Residual variance (scalar).
#' @param sigma_0_2 Prior effect variance (length p vector for the effect).
#' @param mu_0 Prior effect mean (length p vector for the effect).
#' @param xty Optional cached `X^T y` vector on standardized scale.
#' @param d Optional cached `diag(X^T X)` vector on standardized scale.
#'
#' @return Numeric p-vector of log Bayes factors comparing inclusion vs null.
#' @export
#'
#' @examples
BF <- function(X, y, sigma_2, sigma_0_2, mu_0, xty = NULL, d = NULL){

  if (is.null(xty)) xty <- compute_Xty(X, y)
  if (is.null(d)) d <- attr(X, "d")

  beta_hat = xty / d
  shat2 = sigma_2 / d

  log_BF = as.vector(dnorm(beta_hat,mu_0,sqrt(sigma_0_2 + shat2),log = TRUE) -
                       dnorm(beta_hat,0,sqrt(shat2),log = TRUE))

  return(log_BF)

}

#' Log-Sum-Exp Utility
#'
#' @param x Numeric vector.
#'
#' @return Scalar log(sum(exp(x))) computed stably.
#' @export
#'
#' @examples
log_sum_exp <- function(x) {
  max_x <- max(x)
  log_sum_exp_result <- max_x + log(sum(exp(x - max_x)))
  return(log_sum_exp_result)
}
