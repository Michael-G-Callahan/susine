test_that("susine() fits stored fixture data with supported methods", {
  load("test_X_data.Rdata")
  load("test_susine_data.Rdata")

  fit_none <- susine(
    L = L,
    X = X,
    y = y,
    mu_0 = mu_0,
    sigma_0_2 = sigma_0_2,
    prior_inclusion_weights = prior_inclusion_weights,
    prior_update_method = "none"
  )
  fit_var <- susine(
    L = L,
    X = X,
    y = y,
    mu_0 = mu_0,
    sigma_0_2 = sigma_0_2,
    prior_inclusion_weights = prior_inclusion_weights,
    prior_update_method = "var"
  )

  expect_equal(round(fit_none$effect_fits$alpha, 5), round(susine.out$effect_fits$alpha, 5))
  expect_equal(round(fit_none$effect_fits$b_hat, 5), round(susine.out$effect_fits$b_hat, 5))
  expect_equal(round(fit_none$effect_fits$b_2_hat, 5), round(susine.out$effect_fits$b_2_hat, 5))
  expect_equal(round(fit_none$effect_fits$KL, 5), round(susine.out$effect_fits$KL, 5))
  expect_equal(round(fit_none$model_fit$sigma_2, 5), round(susine.out$model_fit$sigma_2, 5))

  expect_true(all(is.finite(fit_var$model_fit$PIPs)))
  expect_true(all(is.finite(fit_var$model_fit$elbo)))
})
