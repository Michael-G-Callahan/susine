test_that("ERSS() works", {
  load('test_X_data.Rdata')
  load('test_ERSS_data.Rdata')
  er2.mc =  ERSS(X, y, b_hat, b_2_hat)
  er2<- sum((y - rowSums(X %*% t(b_hat)))^2) +
    sum(X^2 %*% t(b_2_hat) - (X %*% t(b_hat))^2)
  expect_equal(round(er2,5), round(er2.mc,5))
})

test_that("SER_ERSS() works", {
  load('test_X_data.Rdata')
  load('test_ERSS_data.Rdata')
  b_hat_l <- b_hat[1, ]
  b_2_hat_l <- b_2_hat[1, ]
  ser.er2.mc =  SER_ERSS(X, y, b_hat_l, b_2_hat_l)
  ser.er2 = (sum(y^2) - 2*sum(y * (X %*% b_hat_l))
             + sum(X^2 %*% b_2_hat_l))
  expect_equal(round(ser.er2,5), round(ser.er2.mc,5))
})
