#' SuSiNE with z-scores and LD (RSS interface)
#'
#' @description Convenience wrapper to run `susine_ss` from GWAS z-scores and
#'   LD matrix R. Handles scaling for original or standardized effects.
#'
#' @param L Number of single effects (components).
#' @param z Length-p vector of GWAS z-scores.
#' @param R p x p LD correlation matrix.
#' @param n Sample size used to compute z.
#' @param var_y Optional variance of y on original scale.
#' @param s_hat Optional length-p vector of s.e.(betahat) on original scale.
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
#' @param ld_regularization LD regularization strategy for `R`; one of
#'   `"none"`, `"project"`, `"shrink"`, or `"project_shrink"`.
#' @param ld_shrinkage Shrinkage weight toward identity when using shrinkage.
#' @param ld_min_eig Eigenvalue floor used in PSD projection.
#' @param estimate_residual_variance Logical; if FALSE, keep residual variance
#'   fixed at the standardized RSS default (`sigma_2 = 1` when `var_y` is not
#'   supplied).
#' @param verbose_beta_hat Logical; if TRUE, print per-iteration, per-effect
#'   beta-hat distribution summaries from the summary-stat update path.
#' @param record_beta_hat_diagnostics Logical; if TRUE, store those summaries in
#'   `model_fit$beta_hat_diagnostics`.
#'
#' @return Output from `susine_ss` on constructed `XtX`, `Xty`, `yty`.
#' @export


susine_rss <- function(
    L,
    z,
    R,
    n,
    var_y,
    s_hat,
    mu_0=0,
    sigma_0_2=0,
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
    ld_regularization = c("none", "project", "shrink", "project_shrink"),
    ld_shrinkage = 0.01,
    ld_min_eig = 1e-8,
    estimate_residual_variance = TRUE,
    verbose_beta_hat = FALSE,
    record_beta_hat_diagnostics = FALSE
) {
    convergence_method <- match.arg(convergence_method)
    ld_regularization <- match.arg(ld_regularization)

    R <- regularize_ld_matrix(
      R = R,
      method = ld_regularization,
      shrinkage = ld_shrinkage,
      min_eig = ld_min_eig
    )

    #Get PVE adjusted z scores
    adj = (n-1)/(z^2 + n - 2)
    z   = sqrt(adj) * z

  if (!missing(s_hat) & !missing(var_y)) {

    # var_y, s_hat are provided, so the effects are on the
    # *original scale*.
    XtXdiag = var_y * adj/(s_hat^2)
    XtX = t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
    XtX = (XtX + t(XtX))/2
    Xty = z * sqrt(adj) * var_y / s_hat
  } else {

    # The effects are on the *standardized* X, y scale.
    XtX = (n-1)*R
    Xty = sqrt(n-1)*z
    var_y = 1
  }

  susine_output <- susine_ss(
    L,
    XtX,
    Xty,
    yty=(n-1)*var_y,
    n,
    mu_0,
    sigma_0_2,
    prior_inclusion_weights,
    residual_variance_lowerbound,
    prior_update_method,
    verbose,
    convergence_method,
    tol,
    max_iter,
    anneal_start_T,
    anneal_schedule_type,
    anneal_target,
    anneal_burn_in,
    prior_update_warmup,
    init_random,
    init_random_sd,
    init_seed,
    init_alpha,
    init_effect_fits,
    auto_scale_mu_0,
    estimate_residual_variance = estimate_residual_variance,
    verbose_beta_hat = verbose_beta_hat,
    record_beta_hat_diagnostics = record_beta_hat_diagnostics
  )

  return(susine_output)

}
