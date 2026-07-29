# test-prior-update-methods.R

make_test_data <- function(seed = 99) {
  set.seed(seed)
  n <- 100
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  for (j in 2:p) {
    X[, j] <- 0.4 * X[, j - 1] + sqrt(1 - 0.4^2) * X[, j]
  }
  beta <- rep(0, p)
  beta[c(5, 17)] <- c(0.8, -0.6)
  annotation <- rnorm(p, mean = 0, sd = 0.5)
  annotation[5] <- 0.7
  annotation[17] <- -0.5
  y <- drop(X %*% beta + rnorm(n, sd = 1))
  list(X = X, y = y, beta = beta, annotation = annotation, n = n, p = p)
}

test_that("'none' method leaves priors unchanged", {
  d <- make_test_data()
  fit <- susine(
    L = 3,
    X = d$X,
    y = d$y,
    mu_0 = d$annotation,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    max_iter = 10,
    tol = 1e-5
  )

  for (l in 1:3) {
    s0 <- fit$priors$sigma_0_2[l, ]
    expect_true(all(abs(s0 - s0[1]) < 1e-12))
  }
  expect_true(all(is.finite(fit$model_fit$PIPs)))
})

test_that("'var' updates sigma_0_2 and leaves mu_0 fixed", {
  d <- make_test_data()
  fit_var <- susine(
    L = 3,
    X = d$X,
    y = d$y,
    mu_0 = d$annotation,
    sigma_0_2 = 0.2,
    prior_update_method = "var",
    max_iter = 50,
    tol = 1e-5
  )
  fit_none <- susine(
    L = 3,
    X = d$X,
    y = d$y,
    mu_0 = d$annotation,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    max_iter = 50,
    tol = 1e-5
  )

  expect_equal(fit_var$priors$mu_0, fit_none$priors$mu_0, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fit_var$priors$sigma_0_2, fit_none$priors$sigma_0_2, tolerance = 1e-8)))
})

test_that("'clamped_scale_var' enforces c >= 0 and a variance floor", {
  d <- make_test_data()
  sigma_init <- 0.2
  fit <- susine(
    L = 3,
    X = d$X,
    y = d$y,
    mu_0 = d$annotation,
    sigma_0_2 = sigma_init,
    prior_update_method = "clamped_scale_var",
    max_iter = 100,
    tol = 1e-4
  )

  expect_true(all(is.finite(fit$model_fit$PIPs)))
  expect_true(all(fit$priors$mu_0_scale_factor >= 0))
  expect_gt(min(fit$priors$sigma_0_2), 0)

  # A single EB update must floor c at 0 and sigma_0_2 at 1% of var(y),
  # enforced at the optimizer level regardless of c_nonneg.
  X_std <- initialize_X(d$X)
  y_cent <- d$y - mean(d$y)
  priors <- initialize_priors(
    mu_0 = d$annotation,
    sigma_0_2 = sigma_init,
    prior_inclusion_weights = NULL,
    p = d$p,
    L = 1,
    var_y = var(y_cent),
    X = X_std,
    y = y_cent
  )
  floor_val <- 0.01 * var(y_cent)
  updated <- update_priors(
    priors = priors[1, , drop = FALSE],
    X = X_std,
    y = y_cent,
    sigma_2 = var(y_cent),
    prior_update_method = "clamped_scale_var"
  )
  expect_gte(min(updated$sigma_0_2[[1]]), floor_val * 0.999)
  expect_gte(updated$mu_0_scale_factor[[1]], 0)
})

test_that("annotation-aware methods error when mu_0 is all zeros", {
  d <- make_test_data()

  expect_error(
    susine(
      L = 3,
      X = d$X,
      y = d$y,
      mu_0 = 0,
      sigma_0_2 = 0.2,
      prior_update_method = "clamped_scale_var",
      max_iter = 10
    ),
    "non-zero mu_0"
  )
})

test_that("retired prior_update_method names fail fast", {
  d <- make_test_data()
  retired <- c(
    "mean",
    "mu_var",
    "scale",
    "scale_var",
    "scale_var_outer",
    "penalized_scale_var",
    "alternating_scale_var"
  )

  for (method in retired) {
    expect_error(
      susine(
        L = 3,
        X = d$X,
        y = d$y,
        mu_0 = d$annotation,
        sigma_0_2 = 0.2,
        prior_update_method = method
      ),
      "should be one of"
    )
    expect_error(
      susine_ss(
        L = 1,
        XtX = crossprod(d$X),
        Xty = drop(crossprod(d$X, d$y)),
        yty = sum(d$y^2),
        n = d$n,
        mu_0 = d$annotation,
        sigma_0_2 = 0.2,
        prior_update_method = method
      ),
      "should be one of"
    )
  }

  # The retained annotation-aware method must NOT be rejected.
  expect_error(
    susine(
      L = 3,
      X = d$X,
      y = d$y,
      mu_0 = d$annotation,
      sigma_0_2 = 0.2,
      prior_update_method = "clamped_scale_var",
      max_iter = 5
    ),
    NA
  )
})

test_that("vanilla baseline still matches susieR with 'none'", {
  skip_if_not_installed("susieR")

  set.seed(14)
  n <- 220
  p <- 80
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p)
  beta[c(11, 52)] <- c(0.7, -0.5)
  y <- drop(X %*% beta + rnorm(n))

  fit_susine <- susine(
    L = 5,
    X = X,
    y = y,
    sigma_0_2 = 0.2,
    prior_update_method = "none",
    max_iter = 100,
    tol = 1e-6
  )
  fit_susie <- susieR::susie(
    X = X,
    y = y,
    L = 5,
    scaled_prior_variance = 0.2,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = TRUE,
    max_iter = 100,
    tol = 1e-6
  )

  expect_equal(fit_susine$model_fit$PIPs, fit_susie$pip, tolerance = 5e-3)
})
