#' Aggregate SuSiNE fits across multiple runs
#'
#' @description Combines multiple `susine` (or `susine_ss`) fit objects into a
#'   single lightweight aggregate suitable for downstream metrics like
#'   `report_metrics()`. Aggregation is performed on the posterior inclusion
#'   probabilities (PIPs) and basic diagnostics.
#'
#' @param fits A list of fit objects returned by `susine()` or `susine_ss()`.
#' @param weighting One of `"equal"` (simple average) or `"elbo"` (evidence
#'   weighted by exp(ELBO - max(ELBO))).
#' @param temperature Optional temperature (> 0) to soften/harden ELBO weights
#'   when `weighting = "elbo"`. Higher values make weights more uniform.
#'
#' @return A list with a `model_fit` entry that mirrors the shape expected by
#'   `report_metrics()`, containing at least:
#'   - `PIPs`: aggregated PIPs (length-p numeric vector)
#'   - `elbo`: length-1 numeric vector (aggregated ELBO)
#'   - `sigma_2`: length-1 numeric vector (aggregated residual variance)
#'   Additionally includes `std_coef` if available from inputs.
#' @importFrom utils tail
#' @export
aggregate_susine_fits <- function(fits,
                                  weighting = c("equal","elbo"),
                                  temperature = 1){
  stopifnot(length(fits) >= 1)
  weighting <- match.arg(weighting)
  if (!is.finite(temperature) || temperature <= 0) temperature <- 1

  # Extract per-run quantities safely
  get_pips <- function(f){
    if (!is.null(f$model_fit$PIPs)) return(as.numeric(f$model_fit$PIPs))
    stop("Each fit must contain model_fit$PIPs.")
  }
  get_elbo_last <- function(f){
    x <- f$model_fit$elbo
    if (is.null(x)) return(NA_real_)
    as.numeric(tail(x, 1))
  }
  get_sigma2_last <- function(f){
    x <- f$model_fit$sigma_2
    if (is.null(x)) return(NA_real_)
    as.numeric(tail(x, 1))
  }
  get_std_coef <- function(f){
    if (!is.null(f$model_fit$std_coef)) return(as.numeric(f$model_fit$std_coef))
    NULL
  }

  pips_list <- lapply(fits, get_pips)
  p <- length(pips_list[[1]])
  if (!all(vapply(pips_list, length, integer(1)) == p)){
    stop("All runs must have PIPs of the same length (same p).")
  }

  elbos <- vapply(fits, get_elbo_last, numeric(1))
  sigma2s <- vapply(fits, get_sigma2_last, numeric(1))

  # Build weights
  if (weighting == "equal" || any(!is.finite(elbos))){
    w <- rep(1/length(fits), length(fits))
  } else {
    # Evidence weights with temperature softening
    z <- (elbos - max(elbos, na.rm = TRUE)) / temperature
    z[!is.finite(z)] <- -Inf
    w <- exp(z)
    sw <- sum(w)
    if (sw <= .Machine$double.eps){
      w <- rep(1/length(fits), length(fits))
    } else {
      w <- w / sw
    }
  }

  # Aggregate PIPs and diagnostics
  pip_agg <- Reduce(`+`, Map(function(pi_i, wi) wi * pi_i, pips_list, as.list(w)))
  elbo_agg <- sum(w * ifelse(is.finite(elbos), elbos, 0))
  sigma2_agg <- sum(w * ifelse(is.finite(sigma2s), sigma2s, 0))

  # Optional: aggregate standardized coefficients if present on all runs
  std_list <- lapply(fits, get_std_coef)
  if (all(vapply(std_list, function(x) !is.null(x) && length(x) == p, logical(1)))){
    std_coef_agg <- Reduce(`+`, Map(function(bi, wi) wi * bi, std_list, as.list(w)))
  } else {
    std_coef_agg <- NULL
  }

  model_fit <- list(
    PIPs = pip_agg,
    elbo = c(elbo_agg),
    sigma_2 = c(sigma2_agg)
  )
  if (!is.null(std_coef_agg)) model_fit$std_coef <- std_coef_agg

  # Return a minimal fit-like object so report_metrics() works out of the box
  list(model_fit = model_fit)
}

