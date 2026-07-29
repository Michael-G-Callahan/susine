test_that("alpha convergence mode records and enforces alpha-delta stopping", {
  set.seed(11)
  n <- 160
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(5, 17)] <- c(1.0, -0.8)
  y <- drop(X %*% beta + rnorm(n, sd = 1))

  fit_alpha <- susine(
    L = 5,
    X = X,
    y = y,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    convergence_method = "alpha",
    tol = 1e-3,
    max_iter = 100
  )

  expect_true(length(fit_alpha$model_fit$alpha_diff) >= 1)
  expect_true(all(is.finite(fit_alpha$model_fit$alpha_diff)))
  expect_true(tail(fit_alpha$model_fit$alpha_diff, 1) < 1e-3 + 1e-8)
})

test_that("non-uniform prior inclusion weights propagate correctly", {
  set.seed(12)
  n <- 180
  p <- 40
  X <- susine:::initialize_X(matrix(rnorm(n * p), n, p))
  y <- rep(0, n)

  prior_weights <- exp(seq_len(p) / p)
  prior_weights <- prior_weights / sum(prior_weights)
  uniform_weights <- rep(1 / p, p)

  weighted_priors <- data.frame(
    mu_0 = I(list(rep(0, p))),
    sigma_0_2 = I(list(rep(0.2, p))),
    prior_inclusion_weights = I(list(prior_weights))
  )
  uniform_priors <- data.frame(
    mu_0 = I(list(rep(0, p))),
    sigma_0_2 = I(list(rep(0.2, p))),
    prior_inclusion_weights = I(list(uniform_weights))
  )

  ser_weighted <- SER(X = X, y = y, sigma_2 = 1, priors = weighted_priors)
  ser_uniform <- SER(X = X, y = y, sigma_2 = 1, priors = uniform_priors)

  alpha_weighted <- as.numeric(ser_weighted$alpha[[1]])
  alpha_uniform <- as.numeric(ser_uniform$alpha[[1]])
  lhs <- log(alpha_weighted) - log(alpha_uniform)
  rhs <- log(prior_weights) - log(uniform_weights)

  expect_gt(cor(lhs, rhs), 0.999999)
  expect_equal(lhs - mean(lhs), rhs - mean(rhs), tolerance = 1e-6)
})

test_that("RSS LD regularization yields a valid PSD correlation matrix", {
  set.seed(13)
  p <- 20
  A <- matrix(rnorm(p * p), p, p)
  R <- cov2cor(crossprod(A))
  R[1, 2] <- 1.2
  R[2, 1] <- 1.2
  R <- (R + t(R)) / 2

  R_reg <- susine:::regularize_ld_matrix(
    R = R,
    method = "project_shrink",
    shrinkage = 0.05,
    min_eig = 1e-8
  )
  evals <- eigen(R_reg, symmetric = TRUE, only.values = TRUE)$values

  expect_equal(diag(R_reg), rep(1, p), tolerance = 1e-8)
  expect_true(min(evals) >= -1e-8)
})

test_that("vanilla baseline matches susieR when EB prior updates are disabled", {
  skip_if_not_installed("susieR")

  set.seed(14)
  n <- 220
  p <- 80
  L <- 8
  X <- matrix(rnorm(n * p), n, p)
  for (j in 2:p) {
    X[, j] <- 0.6 * X[, j - 1] + sqrt(1 - 0.6^2) * X[, j]
  }
  beta <- rep(0, p)
  beta[c(9, 35, 61)] <- c(1.0, -0.8, 0.7)
  y <- drop(X %*% beta + rnorm(n, sd = 1.2))

  susie_formals <- names(formals(susieR::susie))
  susie_args <- list(
    X = X,
    y = y,
    L = L,
    scaled_prior_variance = 0.2,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = TRUE,
    standardize = TRUE,
    intercept = TRUE,
    tol = 1e-5,
    max_iter = 100
  )
  if ("convergence_method" %in% susie_formals) {
    susie_args$convergence_method <- "elbo"
  }
  if ("refine" %in% susie_formals) {
    susie_args$refine <- FALSE
  }
  fit_susieR <- do.call(susieR::susie, susie_args)

  fit_susine <- susine(
    L = L,
    X = X,
    y = y,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    auto_scale_mu_0 = FALSE,
    auto_scale_sigma_0_2 = FALSE,
    tol = 1e-5,
    max_iter = 100
  )

  pip_diff <- abs(fit_susine$model_fit$PIPs - fit_susieR$pip)
  pip_cor <- cor(fit_susine$model_fit$PIPs, fit_susieR$pip)
  sigma_rel <- abs(tail(fit_susine$model_fit$sigma_2, 1) - fit_susieR$sigma2) / fit_susieR$sigma2

  expect_lt(max(pip_diff), 0.03)
  expect_gt(pip_cor, 0.99)
  expect_lt(sigma_rel, 0.1)
})
