########################### - Data matrix X - #################################

#TODO -- move all SuSiNE function init lines to an init_setup function

#' Standardize X data
#'
#' @param X (n x p) Data matrix
#'
#' @return X_scaled -- a centered and scaled (n x p) matrix
#' @export
#'
#' @examples
standardize_X = function(X){
  #Center X
  X_col_means <- colMeans(X)
  X_centered <- sweep(X, 2, X_col_means, "-")

  #Scale X
  X_scaled <- apply(X_centered, 2, scale)

  return(X_scaled)
}

#' Initialize the X matrix
#'
#' @param X (n x p) Data matrix
#'
#' @return X - a standardized matrix with column means (cm), column std.
#' deviations (csd), rowSums(X^2) (d) saved as pre-computed attributes.
#' @export
#'
#' @examples
initialize_X = function(X){
  # Add scaling attributes needed for prior standardization
  out = compute_colstats(X,center = TRUE,scale = TRUE)
  orig_center = out$cm
  orig_scale = out$csd

  X = standardize_X(X)
  out = compute_colstats(X,center = TRUE,scale = TRUE)
  attr(X,"scaled:center") = out$cm
  attr(X,"scaled:scale") = out$csd
  attr(X,"original:center") = orig_center
  attr(X,"original:scale") = orig_scale
  attr(X,"d") = out$d

  return(X)
}

#' Initialize the XtX matrix (summary stats)
#'
#' @param XtX p x p cross-product matrix X^T X (original scale).
#' @param n Sample size.
#'
#' @return Standardized `XtX` with attributes `"d" = diag(XtX)` and
#'   `"scaled:scale" = column standard deviations of X.
#' @keywords internal
initialize_XtX = function(XtX, n){
  dXtX = diag(XtX)
  csd = sqrt(dXtX/(n-1))
  csd[csd == 0] = 1
  XtX = t((1/csd) * XtX) / csd

  attr(XtX,"d") = diag(XtX)
  attr(XtX,"scaled:scale") = csd
  return(XtX)
}

#' Regularize LD/correlation matrix for RSS workflows
#'
#' @param R p x p LD/correlation matrix.
#' @param method One of `"none"`, `"project"`, `"shrink"`, `"project_shrink"`.
#' @param shrinkage Shrinkage weight toward identity (used for shrink methods).
#' @param min_eig Minimum eigenvalue floor used in PSD projection.
#'
#' @return Regularized correlation matrix with unit diagonal.
#' @keywords internal
regularize_ld_matrix <- function(R,
                                 method = c("none", "project", "shrink", "project_shrink"),
                                 shrinkage = 0.01,
                                 min_eig = 1e-8) {
  method <- match.arg(method)
  if (method == "none") return(R)

  if (!is.matrix(R) || nrow(R) != ncol(R)) {
    stop("R must be a square matrix for LD regularization.")
  }
  if (!is.finite(shrinkage) || shrinkage < 0 || shrinkage > 1) {
    stop("shrinkage must be in [0,1].")
  }
  if (!is.finite(min_eig) || min_eig < 0) {
    stop("min_eig must be non-negative.")
  }

  out <- (R + t(R)) / 2
  p <- ncol(out)

  if (method %in% c("project", "project_shrink")) {
    ee <- eigen(out, symmetric = TRUE)
    vals <- pmax(ee$values, min_eig)
    out <- ee$vectors %*% (vals * t(ee$vectors))
    out <- (out + t(out)) / 2
  }
  if (method %in% c("shrink", "project_shrink")) {
    out <- (1 - shrinkage) * out + shrinkage * diag(p)
  }

  # Normalize back to a correlation scale.
  d <- sqrt(pmax(diag(out), min_eig))
  out <- sweep(sweep(out, 1, d, "/"), 2, d, "/")
  diag(out) <- 1
  (out + t(out)) / 2
}


########################### - Model settings - #################################

#' Initialize user settings (fixed for model run)
#'
#' @param L Number of single effects (components).
#' @param prior_update_method Strategy for empirical Bayes updates; one of
#'   `"none"`, `"var"`, or `"clamped_scale_var"`.
#' @param convergence_method Convergence criterion; one of `"elbo"` or `"alpha"`.
#' @param tol Convergence tolerance for the selected `convergence_method`.
#' @param max_iter Maximum number of coordinate ascent iterations.
#' @param verbose Logical; print ELBO progress.
#' @param residual_variance_lowerbound Lower bound for residual variance.
#' @param estimate_residual_variance Logical; if FALSE, keep residual variance
#'   fixed at its initial value.
#' @param init_random If TRUE, randomize initial `b_hat` per effect.
#' @param init_random_sd Standard deviation for randomized `b_hat`.
#' @param init_seed Optional integer seed.
#' @param anneal_start_T Starting temperature (>= 1); 1 disables annealing.
#' @param anneal_schedule_type Temperature schedule type: "geometric" or "linear".
#' @param anneal_target Which variances to anneal: "both", "likelihood", or "prior".
#' @param anneal_burn_in Integer; first K iterations to anneal; after this, T=1.
#' @param prior_update_warmup Integer; first K iterations skip EB updates.
#' @param c_nonneg Logical; if TRUE, constrain the estimated annotation scale
#'   factor `c` to be non-negative (`c >= 0`) during scale-based EB updates by
#'   setting the optimizer lower bound to 0 instead of -100. Default FALSE.
#'
#' @return List of immutable settings used during the run.
#' @export
#'
#' @examples
initialize_settings = function(L, prior_update_method, tol, max_iter, verbose, residual_variance_lowerbound,
                              convergence_method = c("elbo","alpha"),
                              estimate_residual_variance = TRUE,
                              init_random = FALSE, init_random_sd = 0.05, init_seed = NULL,
                              anneal_start_T = 1,
                              anneal_schedule_type = c("geometric","linear"),
                              anneal_target = c("both","likelihood","prior"),
                              anneal_burn_in = 0,
                              prior_update_warmup = 0,
                              c_nonneg = FALSE){
  convergence_method <- match.arg(convergence_method)
  anneal_schedule_type <- match.arg(anneal_schedule_type)
  anneal_target <- match.arg(anneal_target)
  temp_sched <- build_temperature_schedule(max_iter, anneal_start_T, anneal_schedule_type, anneal_burn_in)
  settings = list(
    tol = tol,
    convergence_method = convergence_method,
    max_iter = max_iter,
    L = L,
    verbose = verbose,
    prior_update_method = prior_update_method,
    residual_variance_lowerbound = residual_variance_lowerbound,
    estimate_residual_variance = isTRUE(estimate_residual_variance),
    init_random = init_random,
    init_random_sd = init_random_sd,
    init_seed = init_seed,
    anneal_start_T = anneal_start_T,
    anneal_schedule_type = anneal_schedule_type,
    anneal_target = anneal_target,
    anneal_burn_in = anneal_burn_in,
    prior_update_warmup = prior_update_warmup,
    c_nonneg = isTRUE(c_nonneg),
    temperature_schedule = temp_sched
  )
  return(settings)
}

#' Build annealing temperature schedule (internal)
#'
#' @param max_iter Integer, number of iterations.
#' @param start_T Starting temperature (>= 1).
#' @param schedule_type One of "geometric" or "linear".
#' @param burn_in Number of initial iterations to anneal; after this, T=1.
#'
#' @return Numeric vector of length `max_iter` with temperatures per iteration.
#' @keywords internal
build_temperature_schedule <- function(max_iter, start_T, schedule_type = c("geometric","linear"), burn_in = 0){
  schedule_type <- match.arg(schedule_type)
  if (max_iter <= 0) return(numeric(0))
  if (burn_in <= 0 || start_T <= 1) return(rep(1, max_iter))
  burn_in <- min(max_iter, burn_in)
  if (burn_in == 1) {
    Ti_burn <- start_T
  } else {
    w <- seq(0, 1, length.out = burn_in) # 0 at i=1, 1 at i=burn_in
    if (schedule_type == "geometric"){
      Ti_burn <- exp((1 - w) * log(start_T) + w * log(1))
    } else {
      Ti_burn <- (1 - w) * start_T + w * 1
    }
  }
  c(Ti_burn, rep(1, max_iter - burn_in))
}


########################### - Model priors - ###################################

#' Coerce scalar/length-L/length-p/LxP input to an L x p matrix (internal)
#'
#' @param x Input value: scalar, length-L, length-p, or an L x p matrix.
#' @param L Number of effects (rows of the returned matrix).
#' @param p Number of predictors (columns of the returned matrix).
#' @param name Label used in the error message when `x` has invalid dimensions.
#'
#' @return An L x p matrix formed by recycling `x` according to its length.
#' @keywords internal
coerce_to_L_by_p_matrix <- function(x, L, p, name = "param"){
  if (is.matrix(x) && all(dim(x) == c(L, p))){
    return(matrix(x, nrow = L, ncol = p))
  }
  if (length(x) == 1){
    return(matrix(x, nrow = L, ncol = p))
  }
  if (length(x) == p){
    return(matrix(x, nrow = L, ncol = p, byrow = TRUE))
  }
  if (length(x) == L){
    return(matrix(x, nrow = L, ncol = p, byrow = FALSE))
  }
  print(paste0('Warning: ', name, ' does not have acceptable dimensions. Stopping execution.'))
  stop()
}

#' Create prior distribution parameters structure
#'
#' @param mu_0 Prior effect means. Accepts scalar, length-L, length-p, or L x p
#'   matrix; returned as L x p.
#' @param sigma_0_2 Prior effect variances (as proportion of Var(y)). Accepts
#'   scalar, length-L, length-p, or L x p; returned as L x p scaled by `var_y`.
#' @param prior_inclusion_weights Prior inclusion probabilities. Accepts length-L,
#'   length-p, or L x p; returned as L x p with rows summing to ~1.
#' @param p Number of predictors (columns of X).
#' @param L Number of effects (components).
#' @param var_y Variance of y (used to scale `sigma_0_2`).
#' @param X Optional design matrix (standardized scale) used to compute
#'   `beta_hat` and `shat2` if `auto_scale_mu_0=TRUE`.
#' @param y Optional response vector used with `X` for autoscaling.
#' @param Xty Optional p-vector `X^T y` (standardized/original handled upstream)
#'   used for autoscaling in summary-statistics workflows.
#' @param dXtX Optional p-vector `diag(X^T X)` for autoscaling in SS workflows.
#' @param n Optional sample size (not strictly required for autoscaling).
#' @param auto_scale_mu_0 Logical; if TRUE, estimate and store a global
#'   `mu_0_scale_factor = c` and `mu_0_scale_tau2` using the unconditional model.
#' @param auto_scale_sigma_0_2 Logical; if TRUE, estimate and store a global
#'   multiplicative scale `sigma_0_2_scale_factor = k` for the provided per-SNP
#'   prior variances using the unconditional model.
#'
#' @return Data.frame with L rows and columns:
#'   `mu_0`, `sigma_0_2`, `prior_inclusion_weights` (each a length-p list-column),
#'   `mu_0_scale_factor` (length-L vector), and `mu_0_scale_tau2` (length-L vector).
#' @export
#'
#' @examples
initialize_priors = function(mu_0,
                             sigma_0_2,
                             prior_inclusion_weights,
                             p,
                             L,
                             var_y,
                             X = NULL,
                             y = NULL,
                             Xty = NULL,
                             dXtX = NULL,
                             n = NULL,
                             auto_scale_mu_0 = FALSE,
                             auto_scale_sigma_0_2 = FALSE){

  has_values <- function(x) !is.null(x) && length(x) > 0

  user_provided_mu_0 <- has_values(mu_0)
  user_provided_sigma_0_2 <- has_values(sigma_0_2)

  #Prior inclusion weights
  if (!has_values(prior_inclusion_weights)){
    prior_inclusion_weights <- matrix(1/p, nrow = L, ncol = p)
  } else {
    prior_inclusion_weights <- coerce_to_L_by_p_matrix(prior_inclusion_weights, L, p, name = 'prior_inclusion_weights')
  }

  #Prior effect variance (conditional on selection)
  if (user_provided_sigma_0_2) {
    sigma_0_2 <- coerce_to_L_by_p_matrix(sigma_0_2, L, p, name = 'sigma_0_2') * var_y
  } else {
    sigma_0_2 <- matrix(1/L * var_y, nrow = L, ncol = p)
  }

  #Prior scaling factors
  mu_0_scale_factor = rep(1, length.out = L)
  mu_0_scale_tau2 = rep(NA_real_, length.out = L)
  sigma_0_2_scale_factor = rep(1, length.out = L)

  #Prior mean effect mean (conditional on selection)
  if (user_provided_mu_0) {
    mu_0 = coerce_to_L_by_p_matrix(mu_0, L, p, name = 'mu_0')
  } else {
    mu_0 = matrix(0, nrow = L, ncol = p)
  }


  # Scale priors if X is provided and standardized
  if (!is.null(X) && !is.null(attr(X, "original:scale"))) {
    X_sd = attr(X, "original:scale")

    # Scale mu_0 to standardized X space
    if (user_provided_mu_0) {
      if (is.matrix(mu_0) && all(dim(mu_0) == c(L, p))) {
        mu_0 = sweep(mu_0, 2, X_sd, "*")
      } else { # Handles scalar, vector of length p, vector of length L
        mu_0 = mu_0 * X_sd
      }
    }

    # NOTE: sigma_0_2 is NOT scaled by X_sd^2 here. When specified as a
    # proportion of var(y), it is already on the standardized-effect scale
    # (matching susieR's scaled_prior_variance semantics). The mu_0 scaling
    # above IS correct because annotations are in original-X units and need
    # conversion to standardized space. sigma_0_2 does not.
  }

  # Optionally auto-scale prior mean annotations (compute c and tau^2)
  if (isTRUE(auto_scale_mu_0)){
    mu_0_unscaled <- as.vector(mu_0[1, ])
    if (any(mu_0_unscaled != 0)){
      beta_hat <- shat2 <- NULL
      if (!is.null(X) && !is.null(y)){
        d <- attr(X, "d")
        beta_hat <- compute_Xty(X, y) / d
        shat2 <- var_y / d
      } else if (!is.null(Xty) && !is.null(dXtX)){
        beta_hat <- as.vector(Xty / dXtX)
        shat2 <- as.vector(var_y / dXtX)
      }
      if (!is.null(beta_hat) && !is.null(shat2)){
        est <- estimate_mu_0_scale_factor(beta_hat, shat2, mu_0_unscaled)
        mu_0_scale_factor <- rep(est$c, length.out = L)
        mu_0_scale_tau2 <- rep(est$tau2, length.out = L)
      }
    }
  }

  # Optionally auto-scale prior variances (compute global k on log scale)
  if (isTRUE(auto_scale_sigma_0_2) && user_provided_sigma_0_2){
    v0_unscaled <- as.vector(sigma_0_2[1, ])
    if (any(v0_unscaled > 0)){
      beta_hat <- shat2 <- NULL
      if (!is.null(X) && !is.null(y)){
        d <- attr(X, "d")
        beta_hat <- compute_Xty(X, y) / d
        shat2 <- var_y / d
      } else if (!is.null(Xty) && !is.null(dXtX)){
        beta_hat <- as.vector(Xty / dXtX)
        shat2 <- as.vector(var_y / dXtX)
      }
      if (!is.null(beta_hat) && !is.null(shat2)){
        # Use current (possibly scaled) mu_0 as center
        mu0_eff <- as.vector(mu_0[1, ]) * mu_0_scale_factor[1]
        estk <- estimate_sigma0_scale_factor(beta_hat, shat2, mu0_eff, v0_unscaled)
        sigma_0_2_scale_factor <- rep(estk$k, length.out = L)
      }
    }
  }

  # Stash the annotation vector (mu_0 in standardized-X space, before c scaling).

  # This is used by the "clamped_scale_var" method which optimizes c

  # in mu_0 = c * a, where a is this annotation vector.
  annotation <- mu_0  # L x p matrix

  # Apply the global scale factors to the stored priors so downstream code
  # uses pre-scaled values without additional multiplications.
  if (length(mu_0_scale_factor) == L){
    mu_0 <- sweep(mu_0, 1, mu_0_scale_factor, "*")
  }
  if (length(sigma_0_2_scale_factor) == L){
    sigma_0_2 <- sweep(sigma_0_2, 1, sigma_0_2_scale_factor, "*")
  }

  #Create a data frame where every row represents an effect
  priors <- data.frame(
    mu_0               = I(asplit(mu_0, 1)),
    sigma_0_2        = I(asplit(sigma_0_2, 1)),
    prior_inclusion_weights = I(asplit(prior_inclusion_weights, 1)),
    mu_0_scale_factor = mu_0_scale_factor,
    mu_0_scale_tau2 = mu_0_scale_tau2,
    sigma_0_2_scale_factor = sigma_0_2_scale_factor,
    annotation = I(asplit(annotation, 1))
  )

  return(priors)
}

########################### - Model fit - ######################################

#' Initialize effect fit containers (internal)
#'
#' @param L Number of effects (components).
#' @param p Number of predictors.
#' @param settings Settings list (from `initialize_settings`) providing
#'   `init_random`, `init_random_sd`, `init_seed`.
#'
#' @return Data.frame with L rows and columns `alpha`, `b_hat`, `b_2_hat`, `KL`.
#' @importFrom stats rnorm
#' @keywords internal
initialize_effect_fits = function(L,
                                  p,
                                  settings,
                                  init_alpha = NULL,
                                  init_effect_fits = NULL){

  alpha <- matrix(1/p, nrow = L, ncol = p, byrow = TRUE)
  b_hat <- matrix(0, nrow = L, ncol = p)
  b_2_hat <- matrix(0, nrow = L, ncol = p)
  KL <- rep(NA_real_, length.out = L)

  if (!is.null(init_alpha) && !is.null(init_effect_fits)) {
    stop("init_alpha and init_effect_fits cannot both be supplied.")
  }

  if (!is.null(init_effect_fits)) {
    exact_init <- coerce_init_effect_fits(init_effect_fits, L, p)
    alpha <- exact_init$alpha
    b_hat <- exact_init$b_hat
    b_2_hat <- exact_init$b_2_hat
  }

  if (!is.null(init_alpha)){
    alpha <- coerce_to_L_by_p_matrix(init_alpha, L, p, name = 'init_alpha')
    row_sums <- rowSums(alpha)
    if (any(!is.finite(row_sums)) || any(abs(row_sums - 1) > 1e-6)) {
      stop("init_alpha rows must sum to 1.")
    }
    if (any(alpha < 0)) {
      stop("init_alpha must be nonnegative.")
    }
  }

  if (isTRUE(settings$init_random) &&
      is.null(init_alpha) &&
      is.null(init_effect_fits)){
    if (!is.null(settings$init_seed)) set.seed(as.integer(settings$init_seed))
    for (l in 1:L){
      rb <- rnorm(p, mean = 0, sd = settings$init_random_sd)
      b_hat[l, ] <- rb
      b_2_hat[l, ] <- rb^2
    }
  }

  effect_fits <- data.frame(
    alpha = I(asplit(alpha, 1)),
    b_hat = I(asplit(b_hat, 1)),
    b_2_hat = I(asplit(b_2_hat, 1)),
    KL = KL
  )

  return(effect_fits)
}

#' Validate exact effect-fit warm starts (internal)
#'
#' @param init_effect_fits List-like object with `alpha`, `b_hat`, and `b_2_hat`.
#' @param L Number of effects.
#' @param p Number of variables.
#'
#' @return List with validated L x p matrices `alpha`, `b_hat`, and `b_2_hat`.
#' @keywords internal
coerce_init_effect_fits <- function(init_effect_fits, L, p) {
  if (!is.list(init_effect_fits)) {
    stop("init_effect_fits must be a list with alpha, b_hat, and b_2_hat.")
  }

  required <- c("alpha", "b_hat", "b_2_hat")
  missing_fields <- setdiff(required, names(init_effect_fits))
  if (length(missing_fields)) {
    stop(
      "init_effect_fits is missing required field(s): ",
      paste(missing_fields, collapse = ", "),
      "."
    )
  }

  alpha <- coerce_to_L_by_p_matrix(init_effect_fits$alpha, L, p, name = "init_effect_fits$alpha")
  row_sums <- rowSums(alpha)
  if (any(!is.finite(row_sums)) || any(abs(row_sums - 1) > 1e-6)) {
    stop("init_effect_fits$alpha rows must sum to 1.")
  }
  if (any(alpha < 0)) {
    stop("init_effect_fits$alpha must be nonnegative.")
  }

  b_hat <- coerce_to_L_by_p_matrix(init_effect_fits$b_hat, L, p, name = "init_effect_fits$b_hat")
  if (any(!is.finite(b_hat))) {
    stop("init_effect_fits$b_hat must contain only finite values.")
  }

  b_2_hat <- coerce_to_L_by_p_matrix(init_effect_fits$b_2_hat, L, p, name = "init_effect_fits$b_2_hat")
  if (any(!is.finite(b_2_hat))) {
    stop("init_effect_fits$b_2_hat must contain only finite values.")
  }
  if (any(b_2_hat < 0)) {
    stop("init_effect_fits$b_2_hat must be nonnegative.")
  }

  list(
    alpha = alpha,
    b_hat = b_hat,
    b_2_hat = b_2_hat
  )
}

#' Warm-start effect moments from init_alpha (internal)
#'
#' @param effect_fits Data.frame from `initialize_effect_fits`.
#' @param init_alpha Initial alpha matrix (already validated/coerced upstream).
#' @param priors Prior parameters data.frame from `initialize_priors`.
#' @param X Standardized design matrix.
#' @param y Centered response.
#' @param sigma2_init Initial residual variance.
#'
#' @return Updated effect_fits with b_hat and b_2_hat initialized.
#' @keywords internal
warm_start_effect_fits <- function(effect_fits,
                                         init_alpha,
                                         priors,
                                         X,
                                         y,
                                         sigma2_init){
  if (is.null(init_alpha)) return(effect_fits)

  alpha_mat <- do.call(rbind, effect_fits$alpha)
  L <- nrow(alpha_mat)
  p <- ncol(alpha_mat)

  d <- attr(X, "d")
  xty <- compute_Xty(X, y)

  b_hat_init <- matrix(0, nrow = L, ncol = p)
  b_2_hat_init <- matrix(0, nrow = L, ncol = p)

  for (l in 1:L) {
    mu_0_l <- as.vector(priors$mu_0[[l]])
    sigma_0_2_l <- as.vector(priors$sigma_0_2[[l]])
    den <- sigma2_init + sigma_0_2_l * d
    mu_1 <- (sigma2_init * mu_0_l + sigma_0_2_l * xty) / den
    sigma_1_2 <- (sigma2_init * sigma_0_2_l) / den
    b_hat_init[l, ] <- alpha_mat[l, ] * mu_1
    b_2_hat_init[l, ] <- alpha_mat[l, ] * (sigma_1_2 + mu_1^2)
  }

  effect_fits$b_hat <- I(asplit(b_hat_init, 1))
  effect_fits$b_2_hat <- I(asplit(b_2_hat_init, 1))

  return(effect_fits)
}

#' Warm-start effect moments from init_alpha (summary stats, internal)
#'
#' @param effect_fits Data.frame from `initialize_effect_fits`.
#' @param init_alpha Initial alpha matrix (already validated/coerced upstream).
#' @param priors Prior parameters data.frame from `initialize_priors`.
#' @param Xty Standardized X^T y vector (length p).
#' @param dXtX diag(X^T X) vector (length p) on standardized scale.
#' @param sigma2_init Initial residual variance.
#'
#' @return Updated effect_fits with b_hat and b_2_hat initialized.
#' @keywords internal
warm_start_effect_fits_ss <- function(effect_fits,
                                      init_alpha,
                                      priors,
                                      Xty,
                                      dXtX,
                                      sigma2_init){
  if (is.null(init_alpha)) return(effect_fits)

  alpha_mat <- do.call(rbind, effect_fits$alpha)
  L <- nrow(alpha_mat)
  p <- ncol(alpha_mat)

  b_hat_init <- matrix(0, nrow = L, ncol = p)
  b_2_hat_init <- matrix(0, nrow = L, ncol = p)

  for (l in 1:L) {
    mu_0_l <- as.vector(priors$mu_0[[l]])
    sigma_0_2_l <- as.vector(priors$sigma_0_2[[l]])
    den <- sigma2_init + sigma_0_2_l * dXtX
    mu_1 <- (sigma2_init * mu_0_l + sigma_0_2_l * Xty) / den
    sigma_1_2 <- (sigma2_init * sigma_0_2_l) / den
    b_hat_init[l, ] <- alpha_mat[l, ] * mu_1
    b_2_hat_init[l, ] <- alpha_mat[l, ] * (sigma_1_2 + mu_1^2)
  }

  effect_fits$b_hat <- I(asplit(b_hat_init, 1))
  effect_fits$b_2_hat <- I(asplit(b_2_hat_init, 1))

  return(effect_fits)
}

#' Initialize model fit traces (internal)
#'
#' @param max_iter Maximum number of iterations.
#' @param var_y Variance of y (to initialize `sigma_2[1]`).
#' @param p Number of predictors (used to size residuals placeholder).
#'
#' @return List with fields `sigma_2`, `elbo`, `alpha_diff`, and
#'   `full_residuals`.
#' @keywords internal
initialize_model_fit = function(max_iter,
                                var_y,
                                p,
                                sigma2_init = NULL,
                                elbo_init = NULL,
                                residuals_init = NULL){
  sigma_2 <- rep(NA_real_, length.out = max_iter + 1)
  sigma_2[1] <- if (!is.null(sigma2_init)) sigma2_init else var_y

  elbo <- rep(NA_real_, length.out = max_iter + 1)
  elbo[1] <- if (!is.null(elbo_init)) elbo_init else -Inf
  alpha_diff <- rep(NA_real_, length.out = max_iter)

  full_residuals <- if (!is.null(residuals_init)) {
    residuals_init
  } else {
    rep(NA_real_, p)
  }

  model_fit <- list(
    sigma_2 = sigma_2,
    elbo = elbo,
    alpha_diff = alpha_diff,
    full_residuals = full_residuals
  )
  return(model_fit)
}
