#' SuSiNE: Single-Effect Regression with Empirical Bayes updates
#'
#' IBSS overview (algorithmic outline)
#' - Standardize X and center y; initialize priors, effect fits, and model fit.
#' - For each outer iteration i = 1..max_iter:
#'   * Form full residuals from current sum of single-effect fits.
#'   * For each effect l = 1..L:
#'       - Add back effect-l to build its partial residuals.
#'       - Optionally anneal likelihood/prior based on iteration temperature.
#'       - Optionally update priors via EB (after warmup only).
#'       - Run single-effect update (SER) and subtract from residuals.
#'   * Update residual variance and ELBO; check convergence (post burn-in).
#' - Return priors, effect-level summaries, and overall model diagnostics.
#'
#' @description Fits a SuSiE-like model with L single-effect components to
#'   standardized X and y using coordinate-ascent variational inference. Priors
#'   for each effect can be optionally updated via empirical Bayes.
#'
#' @param L Number of single effects (components).
#' @param X Design matrix (n x p). Centered and scaled internally.
#' @param y Response vector of length n.
#' @param mu_0 Prior effect means (scalar, length-L, length-p, or L x p).
#' @param sigma_0_2 Prior effect variances as proportion of Var(y)
#'   (scalar, length-L, length-p, or L x p).
#' @param prior_inclusion_weights Prior inclusion probabilities (length-L,
#'   length-p, or L x p). If NULL, defaults to uniform 1/p.
#' @param residual_variance_lowerbound Lower bound for residual variance.
#' @param prior_update_method Empirical Bayes update strategy; one of
#'   `"none"`, `"var"`, or `"clamped_scale_var"`.
#' @param verbose Logical; print iteration progress.
#' @param convergence_method Convergence criterion; one of `"elbo"` or
#'   `"alpha"` (maximum absolute alpha change).
#' @param tol Convergence tolerance for the selected `convergence_method`.
#' @param max_iter Maximum number of iterations.
#' @param anneal_start_T Starting temperature (>= 1); 1 disables annealing.
#' @param anneal_schedule_type Schedule type: `"geometric"` or `"linear"`.
#' @param anneal_target Which variances to anneal: `"both"`, `"likelihood"`, or
#'   `"prior"`.
#' @param anneal_burn_in Integer; number of initial iterations to anneal.
#' @param prior_update_warmup Integer; number of initial iterations to skip EB
#'   updates (treat as `"none"`).
#' @param c_nonneg Logical; if TRUE, constrain the estimated annotation scale
#'   factor `c` to be non-negative (`c >= 0`) during scale-based EB updates by
#'   setting the optimizer lower bound to 0 instead of -100. Default FALSE.
#' @param init_random If TRUE, randomize initial `b_hat` for each effect.
#' @param init_random_sd Standard deviation for randomized `b_hat` (standardized scale).
#' @param init_seed Optional integer seed for randomized initialization.
#' @param auto_scale_mu_0 Logical; if TRUE, fit a global scaling c and noise
#'   tau^2 for the prior mean annotations via Brent optimization under the
#'   unconditional model m_j ~ N(c a_j, v_j + tau^2), and apply c as
#'   `mu_0_scale_factor`.
#' @param auto_scale_sigma_0_2 Logical; if TRUE, fit a global multiplicative
#'   scale k for provided per-SNP prior variances via unconditional likelihood
#'   m_j ~ N(mu_0j, shat2_j + k v0_j), and apply k as `sigma_0_2_scale_factor`.
#' @param init_alpha Optional L x p (or conformable) matrix giving initial
#'   posterior inclusion probabilities for each effect. If provided, these
#'   alphas are used to initialize effect moments from the current priors and
#'   residual variance.
#' @param init_effect_fits Optional list with conformable `alpha`, `b_hat`, and
#'   `b_2_hat` matrices for exact warm-start initialization. Mutually exclusive
#'   with `init_alpha`.
#' @param init_sigma2 Optional scalar initial residual variance.
#' @param record_history Logical; if TRUE, retain per-iteration PIP history in
#'   `model_fit$pip_history`.
#' @param estimate_residual_variance Logical; if FALSE, keep residual variance
#'   fixed at the initial value (`init_sigma2`, or `var(y)` when `init_sigma2`
#'   is NULL).
#'
#' @return List with elements:
#' - `priors`: L x p prior parameters as matrices/vectors.
#' - `effect_fits`: L x p posterior summaries for each effect.
#' - `model_fit`: overall model diagnostics and summaries (ELBO, sigma_2, PIPs,
#'   coefficients, fitted values).
#' @importFrom stats var
#' @export
#'
#' @examples
susine <- function(
    L,
    X,
    y,
    mu_0=NULL,
    sigma_0_2=NULL,
    prior_inclusion_weights=NULL,
    residual_variance_lowerbound = var(drop(y))/1e4,
    prior_update_method = c("none", "var", "clamped_scale_var"),
    verbose = FALSE,
    convergence_method = c("elbo","alpha"),
    tol = 1e-5,
    max_iter = 100,
    # Annealing controls (disabled by default)
    anneal_start_T = 1,
    anneal_schedule_type = c("geometric","linear"),
    anneal_target = c("both","likelihood","prior"),
    anneal_burn_in = 0,
    prior_update_warmup = 0,
    # Floor c at 0 for scale-based EB methods (FALSE = allow negative c)
    c_nonneg = FALSE,
    # Randomized initialization (disabled by default)
    init_random = FALSE,
    init_random_sd = 0.05,
    init_seed = NULL,
    auto_scale_mu_0 = FALSE,
    auto_scale_sigma_0_2 = FALSE,
    init_alpha = NULL,
    init_effect_fits = NULL,
    init_sigma2 = NULL,
    record_history = FALSE,
    estimate_residual_variance = TRUE
    ) {

  #Initialize immutable user settings
  prior_update_method = match.arg(prior_update_method)
  convergence_method <- match.arg(convergence_method)
  settings = initialize_settings(L, prior_update_method, tol, max_iter, verbose, residual_variance_lowerbound,
                                 convergence_method = convergence_method,
                                 init_random = init_random, init_random_sd = init_random_sd, init_seed = init_seed,
                                 anneal_start_T = anneal_start_T,
                                 anneal_schedule_type = anneal_schedule_type, anneal_target = anneal_target,
                                 anneal_burn_in = anneal_burn_in,
                                 prior_update_warmup = prior_update_warmup,
                                 estimate_residual_variance = estimate_residual_variance,
                                 c_nonneg = c_nonneg)

  #Initialize data
  y = y - mean(y)
  X = initialize_X(X)
  p = dim(X)[2]

  #Priors - fitting optional
  var_y = drop(var(y))
  priors = initialize_priors(
    mu_0, sigma_0_2, prior_inclusion_weights,
    p, L, var_y,
    X = X, y = y,
    auto_scale_mu_0 = auto_scale_mu_0,
    auto_scale_sigma_0_2 = auto_scale_sigma_0_2
  )

  # Validate that annotation-aware methods have non-zero annotations
  annotation_methods <- c("clamped_scale_var")
  if (prior_update_method %in% annotation_methods) {
    anno_check <- do.call(rbind, priors$annotation)
    if (all(anno_check == 0)) {
      stop("prior_update_method '", prior_update_method,
           "' requires a non-zero mu_0 annotation vector.")
    }
  }

  #Fitted values for each effect, and for overall model
  effect_fits = initialize_effect_fits(
    L,
    p,
    settings,
    init_alpha = init_alpha,
    init_effect_fits = init_effect_fits
  )
  model_fit = initialize_model_fit(
    settings$max_iter,
    var_y,
    p,
    sigma2_init = init_sigma2
  )
  model_fit$sigma_2[1] = max(model_fit$sigma_2[1], settings$residual_variance_lowerbound)

  # Warm start effect moments from init_alpha (if provided)
  if (!is.null(init_alpha)) {
    effect_fits <- warm_start_effect_fits(
      effect_fits,
      init_alpha,
      priors,
      X,
      y,
      model_fit$sigma_2[1]
    )
  }

  # Populate initial residuals and ELBO if not already set
  b_hat_full_init = do.call(rbind, effect_fits$b_hat)
  b_2_hat_full_init = do.call(rbind, effect_fits$b_2_hat)
  if (!all(is.na(b_hat_full_init))) {
    model_fit$full_residuals = y - compute_Xb(X, colSums(b_hat_full_init))
  }
  if (is.na(model_fit$elbo[1]) || is.infinite(model_fit$elbo[1])) {
    KL_init = effect_fits$KL
    KL_init[is.na(KL_init)] = 0
    model_fit$elbo[1] = get_objective(
      X,
      y,
      b_hat_full_init,
      b_2_hat_full_init,
      model_fit$sigma_2[1],
      KL_init
    )
  }

  pip_history <- NULL
  if (isTRUE(record_history)) {
    pip_history <- vector("list", length = settings$max_iter + 1)
    alpha_mat_init <- do.call(rbind, effect_fits$alpha)
    pip_history[[1]] <- aggregate_pips(alpha_mat_init)
  }

  # --- IBSS main loop ---
  for (i in 1:settings$max_iter){
    alpha_prev <- do.call(rbind, effect_fits$alpha)

    # Temperature for this iteration (schedule precomputed at init)
    Ti <- settings$temperature_schedule[i]

    # Residuals are carried forward from the prior outer iteration.
    # Recompute only if unavailable.
    if (all(is.na(model_fit$full_residuals))) {
      b_sum <- colSums(do.call(rbind, effect_fits$b_hat))
      model_fit$full_residuals <- y - compute_Xb(X, b_sum)
    }

    for (l in 1:L){
      # Add-back residuals for effect l
      partial_residuals = model_fit$full_residuals + compute_Xb(X, effect_fits$b_hat[[l]])

      # Determine annealed variances
      sigma2_used <- if (settings$anneal_target %in% c("both","likelihood")) model_fit$sigma_2[i] * Ti else model_fit$sigma_2[i]

      # Per-effect priors view (optionally anneal prior variance).
      # We do not persist annealed variance into priors; only EB updates persist.
      pl <- apply_annealing_to_priors_row(priors[l,], Ti, settings$anneal_target)

      # Optional EB prior update (with warmup override)
      pl_updated <- maybe_update_priors_row(pl, X, partial_residuals, sigma2_used, settings, i)
      if (!identical(pl_updated, pl)) priors[l,] <- pl_updated
      pl <- pl_updated

      # Single-effect update under annealed settings
      effect_fits[l,] = SER(
        X,
        partial_residuals,
        sigma2_used,
        pl
      )
      model_fit$full_residuals = partial_residuals - compute_Xb(X, effect_fits$b_hat[[l]])
    }

    # Update ELBO and residual variance estimate
    model_fit = update_model_fit(
      model_fit, i, X, y, effect_fits, settings$residual_variance_lowerbound,
      estimate_residual_variance = settings$estimate_residual_variance
    )
    alpha_curr <- do.call(rbind, effect_fits$alpha)
    model_fit$alpha_diff[i] <- max(abs(alpha_curr - alpha_prev))

    if (isTRUE(record_history)) {
      pip_history[[i + 1]] <- aggregate_pips(do.call(rbind, effect_fits$alpha))
    }

    if (settings$verbose){
      print(paste0('Iteration', i, ' ELBO: ', model_fit$elbo[i+1]))
    }

    conv_delta <- if (settings$convergence_method == "alpha") {
      model_fit$alpha_diff[i]
    } else {
      model_fit$elbo[i+1] - model_fit$elbo[i]
    }

    # Convergence: require at least anneal_burn_in iterations when annealing.
    if (is.finite(conv_delta) &&
        (conv_delta < settings$tol) &&
        (i >= settings$anneal_burn_in)){
      print(paste("Converged in",i,"iterations."))
      break
    }
  }

  if (i == settings$max_iter){
    print("Failed to converge!")
  }

  pip_history_trim <- NULL
  if (isTRUE(record_history) && !is.null(pip_history)) {
    iter_completed <- min(i, settings$max_iter)
    pip_history_trim <- pip_history[seq_len(iter_completed + 1)]
  }

  #Clean up output formats
  priors = finalize_priors(priors)
  effect_fits = finalize_effect_fits(effect_fits)
  model_fit = finalize_model_fit(model_fit, X, effect_fits)
  if (!is.null(pip_history_trim)) {
    model_fit$pip_history = pip_history_trim
  }

  output_list = list(
    priors=priors,
    effect_fits=effect_fits,
    model_fit=model_fit,
    settings=settings
  )
  return(output_list)
}
