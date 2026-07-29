#' Expected Residual Sum of Squares (ERSS)
#'
#' @param X Design matrix (n x p), standardized via `initialize_X` with
#'   attributes `"d"` (rowSums of standardized X^2), `"scaled:center"`, and
#'   `"scaled:scale"` present.
#' @param y Response vector of length n.
#' @param b_hat L x p matrix of per-effect posterior mean coefficients.
#' @param b_2_hat L x p matrix of per-effect posterior second moments
#'   \eqn{E[b^2]} (variance + mean^2).
#'
#' @return Numeric scalar ERSS = \eqn{E[||y - X sum_l b_l||^2]}.
#' @export
#'
#' @examples
#' # Simulate a small standardized design and response
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- initialize_X(matrix(rnorm(n * p), nrow = n, ncol = p))
#' y <- as.numeric(X[, 1] - X[, 2] + rnorm(n))
#' L <- 5
#' b_hat <- matrix(rnorm(L * ncol(X)), nrow = L)
#' b_2_hat <- b_hat^2 + 1  # Example second moments
#' erss <- ERSS(X, y, b_hat, b_2_hat)
ERSS <- function(X, y, b_hat, b_2_hat){

  Xr_L = compute_MXt(b_hat,X)
  ER_2 <- (sum((y - colSums(Xr_L))^2)
           + sum(attr(X,"d") * t(b_2_hat))
           - sum(Xr_L^2))
  return(ER_2)
}

#' SER Posterior Expectation of ||y - X b_l||^2
#'
#' @param X Design matrix (n x p) with standardization attributes.
#' @param y Response vector of length n.
#' @param b_hat Posterior mean vector for a single effect (length p).
#' @param b_2_hat Posterior second moment vector for a single effect
#'   (length p).
#' @param xb_hat Optional cached length-n vector of `X %*% b_hat`.
#'
#' @return Numeric scalar \eqn{E[||y - X b_l||^2]} under SER posterior.
#' @export
#'
#' @examples
SER_ERSS <- function(X, y, b_hat, b_2_hat, xb_hat = NULL){
  if (is.null(xb_hat)) {
    xb_hat <- compute_Xb(X,b_hat)
  }
  ER_2 <- (sum(y*y)
           - 2*sum(y*xb_hat)
           + sum(attr(X,"d") * b_2_hat)
           )

  return(ER_2)
}


#' Evidence Lower Bound (ELBO) Objective
#'
#' @param X Design matrix (n x p) with standardization attributes.
#' @param y Response vector of length n.
#' @param b_hat L x p matrix of posterior means for all effects.
#' @param b_2_hat L x p matrix of posterior second moments for all effects.
#' @param sigma_2 Residual variance estimate (scalar).
#' @param KL Vector of KL contributions (length L) from each effect.
#' @param erss Optional precomputed ERSS scalar.
#'
#' @return Numeric scalar ELBO = \eqn{E[log p(y|b, sigma2)]} - sum_l KL_l.
#' @export
#'
#' @examples
get_objective <- function(X, y, b_hat, b_2_hat, sigma_2, KL, erss = NULL){
  n = nrow(X)

  if (is.null(erss)) {
    erss <- ERSS(X, y, b_hat, b_2_hat)
  }
  E_log_lik = -(n/2) * log(2*pi*sigma_2) - (1/(2*sigma_2)) * erss

  objective = E_log_lik - sum(KL)

  return(objective)
}

#' Update ELBO and Residual Variance (internal)
#'
#' @param model_fit List containing `sigma_2`, `elbo`, and `full_residuals`.
#' @param i Iteration index (integer >= 1).
#' @param X Design matrix (n x p) with standardization attributes.
#' @param y Response vector of length n.
#' @param effect_fits Data frame/list with fields `b_hat`, `b_2_hat`, `KL` for L effects.
#' @param res_var_lwr_bound Lower bound for residual variance (numeric > 0).
#' @param estimate_residual_variance Logical; if FALSE, keep residual variance
#'   fixed at its current value.
#'
#' @return Updated `model_fit` list with `elbo[i+1]` and `sigma_2[i+1]` set.
#' @keywords internal
update_model_fit = function(model_fit, i, X, y, effect_fits, res_var_lwr_bound,
                            estimate_residual_variance = TRUE){
  #Aggregate full model parameters
  b_hat_full = do.call(rbind, effect_fits$b_hat)
  b_2_hat_full = do.call(rbind, effect_fits$b_2_hat)
  erss_i <- ERSS(X, y, b_hat_full, b_2_hat_full)

  #Update objective and residual variance
  model_fit$elbo[i+1] = get_objective(
    X, y, b_hat_full, b_2_hat_full, model_fit$sigma_2[i], effect_fits$KL, erss = erss_i
  )
  model_fit$sigma_2[i+1] <- if (isTRUE(estimate_residual_variance)) {
    pmax(erss_i/nrow(X), res_var_lwr_bound)
  } else {
    model_fit$sigma_2[i]
  }

  return(model_fit)
}
