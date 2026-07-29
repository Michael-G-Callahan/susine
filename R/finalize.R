#' Finalize Priors Structure (internal)
#'
#' @description Converts the per-effect priors data.frame (rows = effects)
#'   into a list of matrices/vectors suitable for output.
#'
#' @param priors Data.frame with columns `mu_0`, `sigma_0_2`,
#'   `prior_inclusion_weights`, `mu_0_scale_factor`, `mu_0_scale_tau2`,
#'   `sigma_0_2_scale_factor`
#'   (each row = one effect).
#'
#' @return List with fields:
#'   `mu_0` (L x p), `sigma_0_2` (L x p),
#'   `prior_inclusion_weights` (L x p), `mu_0_scale_factor` (length L),
#'   `mu_0_scale_tau2` (length L), and `sigma_0_2_scale_factor` (length L).
#' @keywords internal
finalize_priors = function(priors){
  #Convert from data frame into list of matrices
  mu_0 = do.call(rbind, priors$mu_0)
  sigma_0_2 = do.call(rbind, priors$sigma_0_2)
  prior_inclusion_weights = do.call(rbind, priors$prior_inclusion_weights)

  annotation = do.call(rbind, priors$annotation)

  priors <- list(
    mu_0 = mu_0,
    sigma_0_2 = sigma_0_2,
    prior_inclusion_weights = prior_inclusion_weights,
    mu_0_scale_factor = priors$mu_0_scale_factor,
    mu_0_scale_tau2 = priors$mu_0_scale_tau2,
    sigma_0_2_scale_factor = priors$sigma_0_2_scale_factor,
    annotation = annotation
  )

  return(priors)
}

#' Finalize Effect Fits Structure (internal)
#'
#' @description Converts the per-effect fits data.frame (rows = effects)
#'   into a list of matrices/vectors.
#'
#' @param effect_fits Data.frame with columns `alpha`, `b_hat`, `b_2_hat`,
#'   and `KL` (each row = one effect).
#'
#' @return List with fields: `alpha` (L x p), `b_hat` (L x p),
#'   `b_2_hat` (L x p), and `KL` (length L).
#' @keywords internal
finalize_effect_fits = function(effect_fits){
  #Convert from data frame into list of matrices
  alpha = do.call(rbind, effect_fits$alpha)
  b_hat = do.call(rbind, effect_fits$b_hat)
  b_2_hat = do.call(rbind, effect_fits$b_2_hat)

  effect_fits <- list(
    alpha = alpha,
    b_hat = b_hat,
    b_2_hat = b_2_hat,
    KL = effect_fits$KL
  )

  return(effect_fits)
}


#' Finalize Model Fit (full data, internal)
#'
#' @description Adds summaries (PIPs, coefficients, fitted values) to
#'   the `model_fit` list and trims NA entries from traces.
#'
#' @param model_fit List with `sigma_2`, `elbo`, and `full_residuals`.
#' @param X Design matrix (n x p) used for fitting (standardized).
#' @param effect_fits List with `alpha`, `b_hat`, `b_2_hat`, and `KL`.
#'
#' @return Updated `model_fit` list with fields:
#'   `PIPs`, `std_coef`, `coef`, and `fitted_y` added.
#' @keywords internal
finalize_model_fit = function(model_fit, X, effect_fits){
  #Gather model summary information
  PIPs = aggregate_pips(effect_fits$alpha)
  std_coef = colSums(effect_fits$b_hat)
  coef = std_coef / compute_colSds(X) + colMeans(X)
  fitted_y = compute_Xb(X, std_coef)

  model_fit$sigma_2 = model_fit$sigma_2[!is.na(model_fit$sigma_2)]
  model_fit$elbo = model_fit$elbo[!is.na(model_fit$elbo)]
  if (!is.null(model_fit$alpha_diff)) {
    model_fit$alpha_diff = model_fit$alpha_diff[!is.na(model_fit$alpha_diff)]
  }

  model_fit = append(model_fit,
                     list(
                        PIPs = PIPs,
                        std_coef = std_coef,
                        coef = coef,
                        fitted_y = fitted_y
                      )
  )
}

#' Finalize Model Fit (summary stats, internal)
#'
#' @description Adds summaries (PIPs, standardized coefficients, fitted X^T y)
#'   to the `model_fit` list and trims NA entries from traces.
#'
#' @param model_fit List with `sigma_2`, `elbo`, and `full_residuals`.
#' @param XtX p x p cross-product matrix X^T X (standardized scale).
#' @param effect_fits List with `alpha`, `b_hat`, `b_2_hat`, and `KL`.
#'
#' @return Updated `model_fit` list with fields:
#'   `PIPs`, `std_coef`, and `fitted_Xty` added.
#' @keywords internal
finalize_model_fit_ss = function(model_fit, XtX, effect_fits){
  #Gather model summary information
  PIPs = aggregate_pips(effect_fits$alpha)
  std_coef = colSums(effect_fits$b_hat)
  # coef = std_coef / compute_colSds(X) + colMeans(X)
  fitted_Xty = compute_Xb(XtX, std_coef)

  model_fit$sigma_2 = model_fit$sigma_2[!is.na(model_fit$sigma_2)]
  model_fit$elbo = model_fit$elbo[!is.na(model_fit$elbo)]
  if (!is.null(model_fit$alpha_diff)) {
    model_fit$alpha_diff = model_fit$alpha_diff[!is.na(model_fit$alpha_diff)]
  }

  model_fit = append(model_fit,
                     list(
                       PIPs = PIPs,
                       std_coef = std_coef,
                       # coef = coef,
                       fitted_Xty = fitted_Xty
                     )
  )
}

#' Aggregate Posterior Inclusion Probabilities (PIPs)
#'
#' @param alpha L x p matrix of per-effect inclusion probabilities.
#'
#' @return Length-p vector of marginal variable inclusion probabilities.
#' @keywords internal
aggregate_pips = function(alpha){

  if (!is.matrix(alpha)){
    alpha <- as.matrix(alpha)
  }

  include_effects = which(rowSums((alpha - 1/ncol(alpha))^2) > 1e-8)

  if (length(include_effects) == 0){
    return(rep(0, ncol(alpha)))
  }

  alpha = alpha[include_effects, , drop = FALSE]

  complement = 1 - alpha
  column_product = apply(complement, 2, prod)
  PIP = 1 - column_product

  return(PIP)
}
