#' Safe wrapper for prior updates (summary stats, internal)
#'
#' @description Calls `update_priors_ss` with error handling; returns original
#'   `priors` on error.
#'
#' @param priors One-row data.frame/list of prior parameters for an effect.
#' @param Xty p-vector X^T y (standardized or on original scale as provided).
#' @param dXtX p-vector = diag(X^T X).
#' @param sigma_2 Residual variance (scalar).
#' @param prior_update_method Update strategy passed to `update_priors_ss`.
#' @param var_y Phenotypic variance `Var(y)` on the original (unstandardized)
#'   response scale; used to rescale `sigma_0_2` (given as a proportion of
#'   `Var(y)`), initialize the residual variance, and form
#'   `shat2 = var_y / diag(X'X)` in the summary-statistic path.
#'
#' @return Updated priors (same structure as input) or original on error.
#' @keywords internal
safe_update_priors_ss <- function(priors, Xty, dXtX, sigma_2, prior_update_method,
                                   var_y = NULL) {
  tryCatch({
    priors <- update_priors_ss(priors, Xty, dXtX, sigma_2, prior_update_method,
                               var_y = var_y)
    return(priors)
  }, error = function(e) {
    message("Error in update_priors: ", e)
    return(priors)
  })
}


#' Empirical Bayes updates of prior parameters (summary stats)
#'
#' @param priors One-row data.frame/list with fields `mu_0`, `sigma_0_2`,
#'   `prior_inclusion_weights`, `mu_0_scale_factor` for a single effect.
#' @param Xty p-vector X^T y (standardized scale for internal use).
#' @param dXtX p-vector = diag(X^T X).
#' @param sigma_2 Residual variance (scalar).
#' @param prior_update_method One of `"none"`, `"var"`, or `"clamped_scale_var"`.
#' @param var_y Phenotypic variance `Var(y)` on the original (unstandardized)
#'   response scale; used to rescale `sigma_0_2` (given as a proportion of
#'   `Var(y)`), initialize the residual variance, and form
#'   `shat2 = var_y / diag(X'X)` in the summary-statistic path.
#'
#' @return Updated one-row priors in the same structure as input.
#' @export
#'
#' @examples
update_priors_ss <- function(priors, Xty, dXtX, sigma_2, prior_update_method,
                             var_y = NULL){

  #Unpack priors
  list2env(lapply(priors, `[[`, 1), envir = environment())
  p = length(mu_0)

  beta_hat = Xty / dXtX
  shat2 = sigma_2 / dXtX
  if (is.null(var_y)) var_y <- mean(shat2 * dXtX)  # fallback estimate if not supplied

  init_sigma_0_2 = log(max(c((beta_hat-mu_0)^2 - shat2,0.0001))) #not very good guess for updating both, but optim isn't super sensitive to it either
  init_mu_0 = beta_hat[which.max(abs(beta_hat))]
  init_mu_0_scale_factor = 0

  mean_var = c(init_mu_0, init_sigma_0_2)
  scale_var = c(init_mu_0_scale_factor, init_sigma_0_2)

  .optim <- NULL

  if (prior_update_method == "var"){
    proposed_var = optim(
      par = init_sigma_0_2,
      fn = neg_log_lik_var_ss,
      mu_0 = mu_0,
      Xty = Xty,
      dXtX = dXtX,
      sigma_2 = sigma_2,
      prior_inclusion_weights = prior_inclusion_weights,
      method = "L-BFGS-B",
      lower = -30,
      upper = exp(100),
      control = list('maxit' = 10)
    )$par
    proposed_priors = list(
      mu_0 = I(list(mu_0)),
      sigma_0_2 = I(list(rep(exp(proposed_var),p))),
      prior_inclusion_weights = I(list(prior_inclusion_weights)),
      mu_0_scale_factor = mu_0_scale_factor,
      mu_0_scale_tau2 = mu_0_scale_tau2,
      sigma_0_2_scale_factor = if (exists("sigma_0_2_scale_factor")) sigma_0_2_scale_factor else 1,
      annotation = I(list(annotation))
    )
  } else if (prior_update_method == "clamped_scale_var"){
    log_sigma_floor <- log(0.01 * var_y)
    proposed_csv <- optim(
      par = c(mu_0_scale_factor, max(init_sigma_0_2, log_sigma_floor)),
      fn = neg_log_lik_scale_var_ss,
      annotation = annotation,
      Xty = Xty,
      dXtX = dXtX,
      sigma_2 = sigma_2,
      prior_inclusion_weights = prior_inclusion_weights,
      method = "L-BFGS-B",
      lower = c(0, log_sigma_floor),
      upper = c(100, 15),
      control = list(maxit = 10)
    )$par
    proposed_priors = data.frame(
      mu_0 = I(list(proposed_csv[[1]] * annotation)),
      sigma_0_2 = I(list(rep(exp(proposed_csv[[2]]), p))),
      prior_inclusion_weights = I(list(prior_inclusion_weights)),
      mu_0_scale_factor = proposed_csv[[1]],
      mu_0_scale_tau2 = mu_0_scale_tau2,
      sigma_0_2_scale_factor = if (exists("sigma_0_2_scale_factor")) sigma_0_2_scale_factor else 1,
      annotation = I(list(annotation))
    )
  }

  proposal_neg_log_lik = neg_log_lik_ss(proposed_priors, Xty, dXtX, sigma_2)
  current_neg_log_lik = neg_log_lik_ss(priors, Xty, dXtX, sigma_2)

  if (current_neg_log_lik > proposal_neg_log_lik){ #Accept update
    updated_priors = proposed_priors

  } else { #Reject update
    updated_priors = priors
  }

  return(updated_priors)
}

#' Negative log-likelihood under SER with summary stats (internal)
#'
#' @param priors One-row priors for a single effect.
#' @param Xty p-vector X^T y (standardized scale).
#' @param dXtX p-vector = diag(X^T X).
#' @param sigma_2 Residual variance.
#'
#' @return Scalar negative log-likelihood value.
#' @keywords internal
neg_log_lik_ss <- function(priors, Xty, dXtX, sigma_2){

  #Unpack priors
  list2env(lapply(priors, `[[`, 1), envir = environment())

  # Priors are already scaled during initialization
  #Calculate posterior inclusion probabilities
  log_BF = BF_ss(Xty, dXtX, sigma_2, sigma_0_2, mu_0)

  lpo = log_BF + log(prior_inclusion_weights + sqrt(.Machine$double.xmin)) #log of posterior odds - MC edit
  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)

  #Calculate posterior inclusion probabilities
  lbf_model = maxlpo + log(weighted_sum_w)
  log_lik = lbf_model #closest analog in summary case, dropping + sum(dnorm(y,0,sqrt(residual_variance),log = TRUE))

  return(-log_lik)
}

#' Negative log-likelihood parameterization (summary stats): scale factor c and log-variance (internal)
#'
#' @param scale_var_params Length-2 vector: \eqn{[c, log(sigma_0_2)]}.
#' @param annotation Length-p annotation vector (prior mean = c * annotation).
#' @param Xty p-vector X^T y.
#' @param dXtX p-vector = diag(X^T X).
#' @param sigma_2 Residual variance.
#' @param prior_inclusion_weights Length-p prior inclusion weights.
#'
#' @return Scalar negative log-likelihood value.
#' @keywords internal
neg_log_lik_scale_var_ss <- function(scale_var_params, annotation, Xty, dXtX, sigma_2,
                                     prior_inclusion_weights){
  c_val <- scale_var_params[[1]]
  sigma_0_2 <- exp(scale_var_params[[2]])
  mu_0 <- c_val * annotation

  priors <- data.frame(
    mu_0 = I(list(mu_0)),
    sigma_0_2 = I(list(sigma_0_2)),
    prior_inclusion_weights = I(list(prior_inclusion_weights))
  )

  return(neg_log_lik_ss(priors, Xty, dXtX, sigma_2))
}

#' Negative log-likelihood parameterization (summary stats): variance only (internal)
#'
#' @param sigma_0_2 Log prior variance (scalar or length-p).
#' @param mu_0 Length-p prior mean.
#' @param Xty p-vector X^T y.
#' @param dXtX p-vector = diag(X^T X).
#' @param sigma_2 Residual variance.
#' @param prior_inclusion_weights Length-p prior inclusion weights.
#'
#' @return Scalar negative log-likelihood value.
#' @keywords internal
neg_log_lik_var_ss <- function(sigma_0_2, mu_0, Xty, dXtX, sigma_2, prior_inclusion_weights){
  #Unpack priors
  sigma_0_2 = exp(sigma_0_2)

  priors = data.frame(
    mu_0 = I(list(mu_0)),
    sigma_0_2 = I(list(sigma_0_2)),
    prior_inclusion_weights = I(list(prior_inclusion_weights))
  )

  return(neg_log_lik_ss(priors, Xty, dXtX, sigma_2))
}

###
#' Estimate scaling factor for prior mean via Brent optimization (SS)
#'
#' @description Summary-statistics wrapper that estimates scaling factor `c` and
#'   variance `tau^2` using the unconditional model on inputs derived from
#'   `Xty` and `dXtX`.
#'
#' @param Xty Numeric p-vector of `X^T y` on standardized scale.
#' @param dXtX Numeric p-vector of `diag(X^T X)` on standardized scale.
#' @param var_y Scalar variance of y on original scale used for shat2.
#' @param mu_0_unscaled Numeric p-vector of unscaled prior means `a_j`.
#' @param log_tau2_lower Lower bound for `log(tau^2)` in Brent optimization.
#' @param log_tau2_upper Upper bound for `log(tau^2)` in Brent optimization.
#'
#' @return List with elements:
#'   - `c`: estimated scaling factor.
#'   - `tau2`: estimated extra variance component.
#' @keywords internal
estimate_mu_0_scale_factor_ss <- function(Xty,
                                          dXtX,
                                          var_y,
                                          mu_0_unscaled,
                                          log_tau2_lower = -30,
                                          log_tau2_upper = 15){
  beta_hat <- as.vector(Xty / dXtX)
  shat2    <- as.vector(var_y / dXtX)
  estimate_mu_0_scale_factor(beta_hat, shat2, mu_0_unscaled,
                             log_tau2_lower = log_tau2_lower,
                             log_tau2_upper = log_tau2_upper)
}

###
#' Estimate variance scale k via Brent optimization (SS wrapper)
#'
#' @description Summary-statistics wrapper for `estimate_sigma0_scale_factor`.
#'
#' @param Xty Numeric p-vector of `X^T y` on standardized scale.
#' @param dXtX Numeric p-vector of `diag(X^T X)` on standardized scale.
#' @param var_y Scalar variance of y on original scale used for shat2.
#' @param mu0_eff Numeric p-vector of effective prior means used as centers.
#' @param v0_unscaled Numeric p-vector of provided per-SNP prior variances.
#' @param log_k_lower Lower bound for `log(k)` in Brent optimization.
#' @param log_k_upper Upper bound for `log(k)` in Brent optimization.
#'
#' @return List with element `k` (estimated scale factor >= 0).
#' @keywords internal
estimate_sigma0_scale_factor_ss <- function(Xty,
                                            dXtX,
                                            var_y,
                                            mu0_eff,
                                            v0_unscaled,
                                            log_k_lower = -30,
                                            log_k_upper = 15){
  beta_hat <- as.vector(Xty / dXtX)
  shat2    <- as.vector(var_y / dXtX)
  estimate_sigma0_scale_factor(beta_hat, shat2, mu0_eff, v0_unscaled,
                               log_k_lower = log_k_lower,
                               log_k_upper = log_k_upper)
}
