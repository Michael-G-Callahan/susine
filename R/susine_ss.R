
#' SuSiNE with Summary Statistics
#'
#' @description Fits SuSiNE using summary statistics (XtX, Xty, yty, n) instead
#'   of raw X and y. Priors can be updated by empirical Bayes.
#'
#' @param L Number of single effects (components).
#' @param XtX p x p cross-product matrix X^T X (original scale; standardized internally).
#' @param Xty p-vector X^T y (original scale; standardized internally).
#' @param yty Scalar y^T y.
#' @param n Sample size.
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
#' @param init_random If TRUE, randomize initial `b_hat` for each effect.
#' @param init_random_sd Standard deviation for randomized `b_hat` (standardized scale).
#' @param init_seed Optional integer seed for randomized initialization.
#' @param init_alpha Optional L x p (or conformable) matrix giving initial
#'   posterior inclusion probabilities for each effect. If provided, these
#'   alphas are used to initialize effect moments from the current priors and
#'   residual variance.
#' @param init_effect_fits Optional list with conformable `alpha`, `b_hat`, and
#'   `b_2_hat` matrices for exact warm-start initialization. Mutually exclusive
#'   with `init_alpha`.
#' @param auto_scale_mu_0 Logical; if TRUE, fit a global scaling c and noise
#'   tau^2 for the prior mean annotations under the unconditional model and
#'   apply c as `mu_0_scale_factor`.
#' @param auto_scale_sigma_0_2 Logical; if TRUE, fit a global multiplicative
#'   scale k for provided per-SNP prior variances via unconditional likelihood
#'   m_j ~ N(mu_0j, shat2_j + k v0_j), and apply k as `sigma_0_2_scale_factor`.
#' @param estimate_residual_variance Logical; if FALSE, keep residual variance
#'   fixed at its initial standardized-scale value.
#' @param verbose_beta_hat Logical; if TRUE, print per-iteration, per-effect
#'   summaries of the current beta-hat distribution implied by the partial
#'   residuals (`Xty / dXtX`) before each SER update.
#' @param record_beta_hat_diagnostics Logical; if TRUE, store the same
#'   beta-hat summaries in `model_fit$beta_hat_diagnostics`.
#'
#' @return List with `priors`, `effect_fits`, and `model_fit` as in `susine` but
#' for the summary-statistics setting.
#' @importFrom graphics hist
#' @importFrom stats var
#' @export
#'
#' @examples
susine_ss <- function(
    L,
    XtX,
    Xty,
    yty,
    n,
    mu_0=NULL,
    sigma_0_2=NULL,
    prior_inclusion_weights=NULL,
    residual_variance_lowerbound = 1e-4,
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
    # Randomized initialization (disabled by default)
    init_random = FALSE,
    init_random_sd = 0.05,
    init_seed = NULL,
    init_alpha = NULL,
    init_effect_fits = NULL,
    auto_scale_mu_0 = FALSE,
    auto_scale_sigma_0_2 = FALSE,
    estimate_residual_variance = TRUE,
    verbose_beta_hat = FALSE,
    record_beta_hat_diagnostics = FALSE
) {

  summarise_beta_hat_distribution <- function(x){
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) {
      return(data.frame(
        mean = NA_real_,
        median = NA_real_,
        mode = NA_real_,
        q25 = NA_real_,
        q75 = NA_real_,
        q90 = NA_real_,
        q95 = NA_real_,
        variance = NA_real_
      ))
    }
    if (length(unique(x)) == 1) {
      mode_est <- x[[1]]
    } else {
      breaks <- pretty(range(x), n = 40)
      hist_obj <- hist(x, breaks = breaks, plot = FALSE)
      max_bin <- which.max(hist_obj$counts)
      mode_est <- mean(hist_obj$breaks[c(max_bin, max_bin + 1)])
    }
    qq <- stats::quantile(x, probs = c(0.25, 0.75, 0.9, 0.95), na.rm = TRUE, names = FALSE)
    data.frame(
      mean = mean(x),
      median = stats::median(x),
      mode = mode_est,
      q25 = qq[[1]],
      q75 = qq[[2]],
      q90 = qq[[3]],
      q95 = qq[[4]],
      variance = stats::var(x)
    )
  }

  format_beta_hat_diagnostic <- function(iteration, effect_index, diag_row){
    sprintf(
      paste0(
        "[beta_hat] iter=%d effect=%d ",
        "mean=%.6f median=%.6f mode=%.6f ",
        "q25=%.6f q75=%.6f q90=%.6f q95=%.6f var=%.6f"
      ),
      as.integer(iteration),
      as.integer(effect_index),
      diag_row$mean,
      diag_row$median,
      diag_row$mode,
      diag_row$q25,
      diag_row$q75,
      diag_row$q90,
      diag_row$q95,
      diag_row$variance
    )
  }

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
                                 estimate_residual_variance = estimate_residual_variance)

  #Initialize data - TODO
  XtX = initialize_XtX(XtX, n)
  Xty = as.vector(Xty / attr(XtX,"scaled:scale"))
  dXtX = diag(XtX)
  p = dim(XtX)[1]

  #Priors - fitting optional
  var_y = yty/(n-1)
  priors = initialize_priors(
    mu_0, sigma_0_2, prior_inclusion_weights,
    p, L, var_y,
    Xty = Xty, dXtX = dXtX, n = n,
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

  beta_hat_diag_buffer <- NULL
  if (isTRUE(verbose_beta_hat) || isTRUE(record_beta_hat_diagnostics)) {
    beta_hat_diag_buffer <- list()
  }

  #Fitted values for each effect, and for overall model
  effect_fits = initialize_effect_fits(
    L,
    p,
    settings,
    init_alpha = init_alpha,
    init_effect_fits = init_effect_fits
  )
  model_fit = initialize_model_fit(settings$max_iter, var_y, p)

  # Warm start effect moments from init_alpha (if provided)
  if (!is.null(init_alpha)) {
    effect_fits <- warm_start_effect_fits_ss(
      effect_fits,
      init_alpha,
      priors,
      Xty,
      dXtX,
      model_fit$sigma_2[1]
    )
  }

  #SuSiNE
  for (i in 1:settings$max_iter){
    alpha_prev <- do.call(rbind, effect_fits$alpha)

    # Temperature for this iteration
    Ti <- settings$temperature_schedule[i]

    # Residuals are carried across outer iterations; only initialize if needed.
    if (all(is.na(model_fit$full_residuals))) {
      b_hat_full <- colSums(do.call(rbind, effect_fits$b_hat))
      model_fit$full_residuals <- Xty - XtX %*% b_hat_full
    }

    for (l in 1:L){
      partial_residuals = as.vector(model_fit$full_residuals + XtX %*% effect_fits$b_hat[[l]])
      beta_hat_current <- as.vector(partial_residuals / dXtX)
      if (isTRUE(verbose_beta_hat) || isTRUE(record_beta_hat_diagnostics)) {
        beta_hat_diag <- summarise_beta_hat_distribution(beta_hat_current)
        beta_hat_diag$iteration <- as.integer(i)
        beta_hat_diag$effect <- as.integer(l)
        if (!is.null(beta_hat_diag_buffer)) {
          beta_hat_diag_buffer[[length(beta_hat_diag_buffer) + 1L]] <- beta_hat_diag
        }
        if (isTRUE(verbose_beta_hat) && i == 1L) {
          message(format_beta_hat_diagnostic(i, l, beta_hat_diag))
        }
      }

      # Determine annealed variances
      sigma2_used <- if (settings$anneal_target %in% c("both","likelihood")) model_fit$sigma_2[i] * Ti else model_fit$sigma_2[i]

      # Prepare per-effect priors view (optionally anneal prior variance)
      priors_iter <- priors[l,]
      if (settings$anneal_target %in% c("both","prior")){
        s0 <- priors_iter$sigma_0_2[[1]]
        priors_iter$sigma_0_2[[1]] <- s0 * Ti
      }

      # Optionally perform EB prior update (with warmup override)
      pum <- settings$prior_update_method
      if (i <= settings$prior_update_warmup) pum <- "none"
      if (pum != "none"){
        priors_updated <- safe_update_priors_ss(
          priors_iter,
          partial_residuals,
          dXtX,
          sigma2_used,
          pum,
          var_y    = var_y
        )
        priors[l,] <- priors_updated
        priors_for_SER <- priors_updated
      } else {
        priors_for_SER <- priors_iter
      }

      effect_fits[l,] = SER_ss(
        partial_residuals,
        dXtX,
        yty,
        sigma2_used,
        priors_for_SER,
        n
      )

      model_fit$full_residuals = as.vector(partial_residuals - XtX %*% effect_fits$b_hat[[l]])
    }

    #Update ELBO and residual variance estimate
    model_fit = update_model_fit_ss(
      model_fit, i, XtX, Xty, yty, n, effect_fits, settings$residual_variance_lowerbound,
      estimate_residual_variance = settings$estimate_residual_variance
    )
    alpha_curr <- do.call(rbind, effect_fits$alpha)
    model_fit$alpha_diff[i] <- max(abs(alpha_curr - alpha_prev))

    if (settings$verbose){
      print(paste0('Iteration', i, ': ', model_fit$elbo[i+1]))
    }

    conv_delta <- if (settings$convergence_method == "alpha") {
      model_fit$alpha_diff[i]
    } else {
      model_fit$elbo[i+1] - model_fit$elbo[i]
    }

    # Convergence: require at least anneal_burn_in iterations when annealing
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

  if (!is.null(beta_hat_diag_buffer)) {
    beta_hat_diag_df <- if (length(beta_hat_diag_buffer)) {
      do.call(rbind, beta_hat_diag_buffer)
    } else {
      data.frame(
        mean = numeric(),
        median = numeric(),
        mode = numeric(),
        q25 = numeric(),
        q75 = numeric(),
        q90 = numeric(),
        q95 = numeric(),
        variance = numeric(),
        iteration = integer(),
        effect = integer()
      )
    }
    if (isTRUE(record_beta_hat_diagnostics)) {
      model_fit$beta_hat_diagnostics <- beta_hat_diag_df
    }
    if (isTRUE(verbose_beta_hat) && nrow(beta_hat_diag_df) > 0) {
      final_iter <- max(beta_hat_diag_df$iteration, na.rm = TRUE)
      if (is.finite(final_iter) && final_iter > 1) {
        final_rows <- beta_hat_diag_df[beta_hat_diag_df$iteration == final_iter, , drop = FALSE]
        apply(final_rows, 1, function(row) {
          row_df <- as.data.frame(as.list(row), stringsAsFactors = FALSE)
          numeric_cols <- c("mean","median","mode","q25","q75","q90","q95","variance","iteration","effect")
          row_df[numeric_cols] <- lapply(row_df[numeric_cols], as.numeric)
          message(format_beta_hat_diagnostic(row_df$iteration, row_df$effect, row_df))
        })
      }
    }
  }

  #Clean up output formats
  priors = finalize_priors(priors)
  effect_fits = finalize_effect_fits(effect_fits)
  model_fit = finalize_model_fit_ss(model_fit, XtX, effect_fits)

  output_list = list(
    priors=priors,
    effect_fits=effect_fits,
    model_fit=model_fit,
    beta_hat_diagnostics=if (!is.null(model_fit$beta_hat_diagnostics)) {
      model_fit$beta_hat_diagnostics
    } else {
      NULL
    }
  )
  return(output_list)
}
