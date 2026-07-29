#' Safe wrapper for prior updates (internal)
#'
#' @description Calls `update_priors` with error handling; returns original
#'   `priors` on error.
#'
#' @param priors One-row data.frame/list of prior parameters for an effect.
#' @param X Design matrix (n x p).
#' @param partial_residuals Working residual vector of length n.
#' @param sigma_2 Residual variance (scalar).
#' @param prior_update_method Update strategy passed to `update_priors`.
#' @param c_nonneg Logical; if TRUE, constrain the estimated annotation scale
#'   factor `c` to be non-negative (`c >= 0`) during scale-based EB updates by
#'   setting the optimizer lower bound to 0 instead of -100. Default FALSE.
#'
#' @return Updated priors (same structure as input) or original on error.
#' @keywords internal
safe_update_priors <- function(priors, X, partial_residuals, sigma_2, prior_update_method,
                               c_nonneg = FALSE) {
  tryCatch({
    priors <- update_priors(priors, X, partial_residuals, sigma_2, prior_update_method,
                            c_nonneg = c_nonneg)
    return(priors)
  }, error = function(e) {
    message("Error in update_priors: ", e)
    return(priors)
  })
}

#' Apply annealing to a one-row priors view (internal)
#'
#' @param priors_row One-row data.frame of priors (for a single effect).
#' @param Ti Temperature scalar for this iteration.
#' @param target Anneal target: one of "both","prior","likelihood".
#'
#' @return One-row data.frame with sigma_0_2 scaled by Ti if target includes prior.
#' @keywords internal
apply_annealing_to_priors_row <- function(priors_row, Ti, target){
  if (target %in% c("both","prior")){
    s0 <- priors_row$sigma_0_2[[1]]
    priors_row$sigma_0_2[[1]] <- s0 * Ti
  }
  priors_row
}

#' Maybe update priors with EB given warmup (internal)
#'
#' @param priors_row One-row priors for current effect.
#' @param X Design matrix.
#' @param y_resid Partial residuals for current effect.
#' @param sigma2_used Residual variance used this iteration (possibly annealed).
#' @param settings Settings list with prior_update_method and warmup.
#' @param iter_i Current iteration index (1-based).
#'
#' @return One-row priors (updated if EB performed; otherwise original).
#' @keywords internal
maybe_update_priors_row <- function(priors_row, X, y_resid, sigma2_used, settings, iter_i){
  pum <- settings$prior_update_method
  if (iter_i <= settings$prior_update_warmup) pum <- "none"
  if (pum == "none") return(priors_row)
  safe_update_priors(priors_row, X, y_resid, sigma2_used, pum,
                     c_nonneg = isTRUE(settings$c_nonneg))
}


#' Empirical Bayes updates of prior parameters (full data)
#'
#' @param priors One-row data.frame/list with fields `mu_0`, `sigma_0_2`,
#'   `prior_inclusion_weights`, `mu_0_scale_factor` for a single effect.
#' @param X Design matrix (n x p) with standardization attributes.
#' @param y Working residual vector (length n) for the current effect.
#' @param sigma_2 Residual variance (scalar).
#' @param prior_update_method One of `"none"`, `"var"`, or `"clamped_scale_var"`.
#' @param c_nonneg Logical; if TRUE, constrain the estimated annotation scale
#'   factor `c` to be non-negative (`c >= 0`) during scale-based EB updates by
#'   setting the optimizer lower bound to 0 instead of -100. Default FALSE.
#'
#' @return Updated one-row priors in the same structure as input.
#' @importFrom stats optim var dnorm
#' @export
#'
#' @examples
update_priors <- function(priors, X, y, sigma_2, prior_update_method,
                          c_nonneg = FALSE){

  #Unpack priors
  list2env(lapply(priors, `[[`, 1), envir = environment())
  p = length(mu_0)

  beta_hat = compute_Xty(X, y) / attr(X, "d")
  shat2 = sigma_2 / attr(X, "d")
  var_y <- var(y)

  init_sigma_0_2 = log(max(c((beta_hat-mu_0)^2 - shat2,1))) #not very good guess for updating both, but optim isn't super sensitive to it either
  init_mu_0 = beta_hat[which.max(abs(beta_hat))]
  init_mu_0_scale_factor = 1

  mean_var = list(init_mu_0, init_sigma_0_2)
  scale_var = list(init_mu_0_scale_factor, init_sigma_0_2)

  .optim <- NULL

  if (prior_update_method == "var"){
    proposed_var = optim(
      par = init_sigma_0_2,
      fn = neg_log_lik_var,
      mu_0 = mu_0,
      X = X,
      y = y,
      sigma_2 = sigma_2,
      prior_inclusion_weights = prior_inclusion_weights,
      method = "Brent",
      lower = -30,
      upper = 15
    )$par
    proposed_priors = data.frame(
      mu_0 = I(list(mu_0)),
      sigma_0_2 = I(list(rep(exp(proposed_var),p))),
      prior_inclusion_weights = I(list(prior_inclusion_weights)),
      mu_0_scale_factor = mu_0_scale_factor,
      mu_0_scale_tau2 = mu_0_scale_tau2,
      sigma_0_2_scale_factor = sigma_0_2_scale_factor,
      annotation = I(list(annotation))
    )
  } else if (prior_update_method == "clamped_scale_var"){
    log_sigma_floor <- log(0.01 * var_y)
    proposed_csv <- optim(
      par = c(mu_0_scale_factor, max(init_sigma_0_2, log_sigma_floor)),
      fn = neg_log_lik_scale_var,
      annotation = annotation,
      X = X,
      y = y,
      sigma_2 = sigma_2,
      prior_inclusion_weights = prior_inclusion_weights,
      method = "L-BFGS-B",
      lower = c(0, log_sigma_floor),
      upper = c(100, 15),
      control = list(maxit = 20)
    )$par
    proposed_priors = data.frame(
      mu_0 = I(list(proposed_csv[[1]] * annotation)),
      sigma_0_2 = I(list(rep(exp(proposed_csv[[2]]), p))),
      prior_inclusion_weights = I(list(prior_inclusion_weights)),
      mu_0_scale_factor = proposed_csv[[1]],
      mu_0_scale_tau2 = mu_0_scale_tau2,
      sigma_0_2_scale_factor = sigma_0_2_scale_factor,
      annotation = I(list(annotation))
    )
  }

  proposal_neg_log_lik = neg_log_lik(proposed_priors, X, y, sigma_2)
  current_neg_log_lik = neg_log_lik(priors, X, y, sigma_2)

  # Check against null (sigma=0) likelihood — skipped for methods that explicitly
  # constrain or penalize sigma away from zero, as the null check would override them.
  skip_null_check <- prior_update_method %in%
    c("clamped_scale_var")

  if (!skip_null_check) {
    null_priors <- proposed_priors
    null_priors$sigma_0_2[[1]][] <- 0
    null_neg_log_lik <- neg_log_lik(null_priors, X, y, sigma_2)
    # If null is better, or negligibly worse, set variance to 0
    if (null_neg_log_lik <= proposal_neg_log_lik + 1e-6) {
      proposed_priors$sigma_0_2[[1]][] <- 0
      proposal_neg_log_lik <- null_neg_log_lik
    }
  }

  if (current_neg_log_lik > proposal_neg_log_lik){ #Accept update
    updated_priors = proposed_priors

  } else { #Reject update
    updated_priors = priors
  }

  return(updated_priors)
}

#' Negative log-likelihood under SER given priors (internal)
#'
#' @param priors One-row priors for a single effect.
#' @param X Design matrix (n x p).
#' @param y Working residual vector (length n).
#' @param sigma_2 Residual variance.
#'
#' @return Scalar negative log-likelihood value.
#' @keywords internal
neg_log_lik <- function(priors, X, y, sigma_2){

  #Unpack priors
  list2env(lapply(priors, `[[`, 1), envir = environment())

  # Priors are already scaled during initialization
  #Calculate posterior inclusion probabilities
  log_BF = BF(X, y, sigma_2, sigma_0_2, mu_0) #log of Bayesian Factor

  lpo = log_BF + log(prior_inclusion_weights + sqrt(.Machine$double.xmin)) #log of posterior odds - MC edit
  maxlpo = max(lpo)
  w_weighted = exp(lpo - maxlpo)
  weighted_sum_w = sum(w_weighted)

  #Calculate posterior inclusion probabilities
  lbf_model = maxlpo + log(weighted_sum_w)
  log_lik = lbf_model + sum(dnorm(y,0,sqrt(sigma_2),log = TRUE))

  return(-log_lik)
}

#' Negative log-likelihood parameterization: scale factor c and log-variance (internal)
#'
#' @description Parameterizes the prior mean as `c * annotation` and variance
#'   as `exp(log_sigma_0_2)`, for joint optimization of `c` and `sigma_0_2`.
#'
#' @param scale_var_params Length-2 vector: \eqn{[c, log(sigma_0_2)]}.
#' @param annotation Length-p annotation vector (prior mean = c * annotation).
#' @param X Design matrix (n x p).
#' @param y Working residual vector (length n).
#' @param sigma_2 Residual variance.
#' @param prior_inclusion_weights Length-p prior inclusion weights.
#'
#' @return Scalar negative log-likelihood value.
#' @keywords internal
neg_log_lik_scale_var <- function(scale_var_params, annotation, X, y, sigma_2,
                                  prior_inclusion_weights){
  c_val <- scale_var_params[[1]]
  sigma_0_2 <- exp(scale_var_params[[2]])
  mu_0 <- c_val * annotation

  priors <- data.frame(
    mu_0 = I(list(mu_0)),
    sigma_0_2 = I(list(sigma_0_2)),
    prior_inclusion_weights = I(list(prior_inclusion_weights))
  )

  return(neg_log_lik(priors, X, y, sigma_2))
}

#' Negative log-likelihood parameterization: variance only (internal)
#'
#' @param sigma_0_2 Log prior variance (scalar or length-p).
#' @param mu_0 Length-p prior mean.
#' @param X Design matrix (n x p).
#' @param y Working residual vector (length n).
#' @param sigma_2 Residual variance.
#' @param prior_inclusion_weights Length-p prior inclusion weights.
#'
#' @return Scalar negative log-likelihood value.
#' @keywords internal
neg_log_lik_var <- function(sigma_0_2, mu_0, X, y, sigma_2, prior_inclusion_weights){
  #Unpack priors
  sigma_0_2 = exp(sigma_0_2)

  priors = data.frame(
    mu_0 = I(list(mu_0)),
    sigma_0_2 = I(list(sigma_0_2)),
    prior_inclusion_weights = I(list(prior_inclusion_weights))
  )

  return(neg_log_lik(priors, X, y, sigma_2))
}

###
#' Estimate scaling factor for prior mean via Brent optimization (full data)
#'
#' @description Estimates a scaling factor `c` for the prior mean `mu_0` by
#'   maximizing the marginal log-likelihood of the unconditional regression
#'   model: `m_j ~ N(c * a_j, v_j + tau^2)`, where `m_j` are univariate
#'   regression coefficients (`beta_hat`), `a_j` are annotations (the unscaled
#'   prior means `mu_0_unscaled`), `v_j` are their variances (`shat2`), and
#'   `tau^2` is an extra variance component. The optimization is performed over
#'   `tau^2` (on the log scale).
#'
#' @param beta_hat Numeric vector of univariate regression coefficients `m_j`.
#' @param shat2 Numeric vector of corresponding variances `v_j`.
#' @param mu_0_unscaled Numeric vector of unscaled prior means `a_j`.
#' @param log_tau2_lower Lower bound for `log(tau^2)` in Brent optimization.
#' @param log_tau2_upper Upper bound for `log(tau^2)` in Brent optimization.
#'
#' @return List with elements:
#'   - `c`: estimated scaling factor.
#'   - `tau2`: estimated extra variance component.
#' @keywords internal
estimate_mu_0_scale_factor <- function(beta_hat,
                                       shat2,
                                       mu_0_unscaled,
                                       log_tau2_lower = -30,
                                       log_tau2_upper = 15) {

  # Early exits and guards
  if (length(beta_hat) != length(shat2) || length(beta_hat) != length(mu_0_unscaled)){
    stop("beta_hat, shat2, and mu_0_unscaled must have the same length.")
  }
  if (all(mu_0_unscaled == 0) || any(!is.finite(beta_hat)) || any(!is.finite(shat2))) {
    return(list(c = 1, tau2 = 0))
  }

  # Negative log-likelihood as a function of log(tau^2)
  neg_log_lik_tau2 <- function(log_tau2, beta_hat, shat2, mu_0_unscaled) {
    tau2 <- exp(log_tau2)
    denom <- shat2 + tau2
    denom[denom <= .Machine$double.eps] <- .Machine$double.eps
    w <- 1 / denom

    sum_m_a_w <- sum(beta_hat * mu_0_unscaled * w)
    sum_a2_w  <- sum(mu_0_unscaled^2 * w)
    c_hat <- if (sum_a2_w > 0) sum_m_a_w / sum_a2_w else 0

    residuals <- beta_hat - c_hat * mu_0_unscaled
    # Include normalization constant for clarity (does not affect argmin)
    nll <- 0.5 * sum(log(2 * pi * denom) + residuals^2 / denom)
    return(nll)
  }

  # Initialize at a sensible log(tau^2)
  init_log_tau2 <- log(max(stats::var(beta_hat, na.rm = TRUE) - mean(shat2, na.rm = TRUE), 1e-6))
  init_log_tau2[!is.finite(init_log_tau2)] <- 0

  opt_result <- optim(
    par = init_log_tau2,
    fn = neg_log_lik_tau2,
    beta_hat = beta_hat,
    shat2 = shat2,
    mu_0_unscaled = mu_0_unscaled,
    method = "Brent",
    lower = log_tau2_lower,
    upper = log_tau2_upper
  )

  optimal_tau2 <- exp(opt_result$par)
  denom <- shat2 + optimal_tau2
  denom[denom <= .Machine$double.eps] <- .Machine$double.eps
  w_opt <- 1 / denom
  sum_m_a_w_opt <- sum(beta_hat * mu_0_unscaled * w_opt)
  sum_a2_w_opt  <- sum(mu_0_unscaled^2 * w_opt)
  c_hat_opt <- if (sum_a2_w_opt > 0) sum_m_a_w_opt / sum_a2_w_opt else 0

  return(list(c = c_hat_opt, tau2 = optimal_tau2))
}
###
#' Estimate scaling factor for per-SNP prior variances via Brent optimization (full data)
#'
#' @description Estimates a global multiplicative factor `k` for the provided
#'   per-SNP prior variances `v0_j` by maximizing the unconditional marginal
#'   likelihood under: m_j ~ N(mu0_j, shat2_j + k * v0_j). Optimization is over
#'   log(k) with simple Brent search. Returns k on natural scale.
#'
#' @param beta_hat Numeric vector of univariate regression coefficients m_j.
#' @param shat2 Numeric vector of corresponding variances s_j^2.
#' @param mu0_eff Numeric vector of effective prior means (after any scaling).
#' @param v0_unscaled Numeric vector of provided per-SNP prior variances.
#' @param log_k_lower Lower bound for log(k).
#' @param log_k_upper Upper bound for log(k).
#'
#' @return List with element `k` (estimated scale factor >= 0).
#' @keywords internal
estimate_sigma0_scale_factor <- function(beta_hat,
                                         shat2,
                                         mu0_eff,
                                         v0_unscaled,
                                         log_k_lower = -30,
                                         log_k_upper = 15){
  stopifnot(length(beta_hat) == length(shat2),
            length(beta_hat) == length(mu0_eff),
            length(beta_hat) == length(v0_unscaled))

  if (all(!is.finite(v0_unscaled)) || all(v0_unscaled <= 0)){
    return(list(k = 1))
  }

  nll_logk <- function(logk, m, s2, mu, v0){
    k <- exp(logk)
    denom <- s2 + k * pmax(v0, 0)
    denom[denom <= .Machine$double.eps] <- .Machine$double.eps
    res <- m - mu
    0.5 * sum(log(2*pi*denom) + (res^2)/denom)
  }

  init_logk <- 0
  opt <- optim(par = init_logk,
               fn = nll_logk,
               m = beta_hat,
               s2 = shat2,
               mu = mu0_eff,
               v0 = v0_unscaled,
               method = "Brent",
               lower = log_k_lower,
               upper = log_k_upper)
  list(k = as.numeric(exp(opt$par)))
}
