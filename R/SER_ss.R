#' Single Effect Regression (SER) with Summary Statistics
#'
#' @param Xty p-vector of marginal associations X^T y (standardized scale).
#' @param dXtX p-vector = diag(X^T X) (standardized scale).
#' @param yty Scalar y^T y.
#' @param sigma_2 Residual variance (scalar).
#' @param priors One-row data.frame/list of prior parameters for the effect:
#'   `mu_0` (length p), `sigma_0_2` (length p), `prior_inclusion_weights`
#'   (length p), and `mu_0_scale_factor` (scalar).
#' @param n Sample size.
#'
#' @return List with:
#' - alpha: length-p vector of posterior inclusion probabilities.
#' - b_hat: length-p vector of posterior mean effect.
#' - b_2_hat: length-p vector of posterior second moment \eqn{E[b^2]}.
#' - KL: KL contribution for this effect.
#' @export
#'
#' @examples
SER_ss <- function(Xty, dXtX, yty, sigma_2, priors, n){
  #Unpack priors
  list2env(lapply(priors, `[[`, 1), envir = environment())
  # mu_0 and sigma_0_2 are assumed pre-scaled at initialization

  #Calculate conditional posterior means
  mu_1 = as.vector(
    (sigma_2/sigma_0_2)*mu_0 + Xty
  )/(
    (sigma_2/sigma_0_2) + dXtX)

  #Calculate conditional posterior variances
  sigma_1_2 = sigma_2 /(
    (sigma_2/sigma_0_2) + dXtX)

  #Calculate posterior inclusion probabilities
  log_BF = BF_ss(Xty, dXtX, sigma_2, sigma_0_2, mu_0) #log of Bayesian Factor

  lpo = log_BF + log(prior_inclusion_weights + sqrt(.Machine$double.eps)) #log of posterior odds
  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)

  alpha = w_weighted / weighted_sum_w

  #Calculate unconditional posterior effect moments
  b_hat = as.matrix(alpha* mu_1)
  b_2_hat = as.matrix(alpha*(sigma_1_2 + mu_1^2))

  #Calculate the effect's contribution to KL objective (third term)
  lbf_model = maxlpo + log(weighted_sum_w)
  KL = -lbf_model -
    ((1/(2*sigma_2))*SER_ERSS_ss(dXtX, Xty, yty, b_hat, b_2_hat))


  output_list = list(
    alpha = I(list(alpha)),
    b_hat = I(list(drop(b_hat))),
    b_2_hat = I(list(drop(b_2_hat))),
    KL = KL
  )

  return(output_list)
}

#' Bayesian Factor with Summary Statistics (log-space)
#'
#' @param Xty p-vector of marginal associations X^T y.
#' @param dXtX p-vector = diag(X^T X).
#' @param sigma_2 Residual variance (scalar).
#' @param sigma_0_2 Prior effect variance (length p vector for the effect).
#' @param mu_0 Prior effect mean (length p vector for the effect).
#'
#' @return Numeric p-vector of log Bayes factors for inclusion vs null.
#' @export
#'
#' @examples
BF_ss <- function(Xty, dXtX, sigma_2, sigma_0_2, mu_0){

  beta_hat = Xty / dXtX
  shat2 = sigma_2 / dXtX

  log_BF = as.vector(dnorm(beta_hat,mu_0,sqrt(sigma_0_2 + shat2),log = TRUE) -
                       dnorm(beta_hat,0,sqrt(shat2),log = TRUE))

  return(log_BF)

}
