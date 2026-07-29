#' Expected Residual Sum of Squares (summary stats)
#'
#' @param XtX p x p cross-product matrix X^T X (on standardized scale).
#' @param Xty p-vector X^T y (on standardized scale).
#' @param yty Scalar y^T y.
#' @param b_hat L x p matrix of posterior mean effects.
#' @param b_2_hat L x p matrix of posterior second moments \eqn{E[b^2]}.
#'
#' @return Numeric scalar ERSS = \eqn{E[||y - X sum_l b_l||^2]}.
#' @export
#'
#' @examples
ERSS_ss <- function(XtX, Xty, yty, b_hat, b_2_hat){

  beta_bar = colSums(b_hat)
  XB2 = sum((b_hat %*% XtX) * b_hat)

  ER_2 = (yty - 2*sum(beta_bar * Xty)
          + sum(beta_bar * (XtX %*% beta_bar))
          - XB2
          + sum(diag(XtX) * t(b_2_hat)))

  return(ER_2)
}

#' SER Posterior Expectation of ||y - X b_l||^2 (summary stats)
#'
#' @param dXtX p-vector = diag(X^T X) (standardized scale).
#' @param Xty p-vector X^T y (standardized scale).
#' @param yty Scalar y^T y (unused here but kept for symmetry).
#' @param b_hat Posterior mean vector for a single effect (length p).
#' @param b_2_hat Posterior second moment vector for a single effect (length p).
#'
#' @return Numeric scalar \eqn{E[||y - X b_l||^2]} using summary statistics.
#' @export
#'
#' @examples
SER_ERSS_ss <- function(dXtX, Xty, yty, b_hat, b_2_hat){
  ER_2 <- (- 2*sum(b_hat*Xty) + sum(dXtX * as.vector(b_2_hat)))
  return(ER_2)
}


#' ELBO Objective (summary stats)
#'
#' @param XtX p x p cross-product matrix X^T X (standardized scale).
#' @param Xty p-vector X^T y (standardized scale).
#' @param yty Scalar y^T y.
#' @param n Sample size.
#' @param b_hat L x p matrix of posterior means for all effects.
#' @param b_2_hat L x p matrix of posterior second moments for all effects.
#' @param sigma_2 Residual variance estimate (scalar).
#' @param KL Vector of KL contributions (length L).
#' @param erss Optional precomputed ERSS scalar.
#'
#' @return Numeric scalar ELBO.
#' @export
#'
#' @examples
get_objective_ss <- function(XtX, Xty, yty, n, b_hat, b_2_hat, sigma_2, KL, erss = NULL){
  if (is.null(erss)) {
    erss <- ERSS_ss(XtX, Xty, yty, b_hat, b_2_hat)
  }
  E_log_lik = -(n/2) * log(2*pi*sigma_2) - (1/(2*sigma_2)) * erss

  objective = E_log_lik - sum(KL)

  return(objective)
}

#' Update ELBO and Residual Variance (summary stats, internal)
#'
#' @param model_fit List with `sigma_2`, `elbo`, `full_residuals`.
#' @param i Iteration index (integer >= 1).
#' @param XtX p x p cross-product matrix X^T X (standardized scale).
#' @param Xty p-vector X^T y (standardized scale).
#' @param yty Scalar y^T y.
#' @param n Sample size.
#' @param effect_fits Data frame/list with fields `b_hat`, `b_2_hat`, `KL`.
#' @param res_var_lwr_bound Lower bound for residual variance.
#' @param estimate_residual_variance Logical; if FALSE, keep residual variance
#'   fixed at its current value.
#'
#' @return Updated `model_fit` list with `elbo[i+1]` and `sigma_2[i+1]` set.
#' @keywords internal
update_model_fit_ss = function(model_fit, i, XtX, Xty, yty, n, effect_fits,
                               res_var_lwr_bound,
                               estimate_residual_variance = TRUE){
  #Aggregate full model parameters
  b_hat_full = do.call(rbind, effect_fits$b_hat)
  b_2_hat_full = do.call(rbind, effect_fits$b_2_hat)
  erss_i <- ERSS_ss(XtX, Xty, yty, b_hat_full, b_2_hat_full)

  #Update objective and residual variance
  model_fit$elbo[i+1] = get_objective_ss(
    XtX, Xty, yty, n, b_hat_full, b_2_hat_full, model_fit$sigma_2[i], effect_fits$KL, erss = erss_i
  )
  model_fit$sigma_2[i+1] <- if (isTRUE(estimate_residual_variance)) {
    pmax(erss_i/n, res_var_lwr_bound)
  } else {
    model_fit$sigma_2[i]
  }

  return(model_fit)
}
