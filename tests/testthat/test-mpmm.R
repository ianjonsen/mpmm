## mpmm output is a list with class `fG_mpm`
##  tests expect that mpmm are 13-element lists (if optimiser does not crash)
##  that have S3 class mpmm
data(ellie.ice.short)

## REML on
test_that("mpmm returns mpmm list w 13 elements", {
  fit <- mpmm(~ ice + (1 | id), data = ellie.ice.short, control = mpmm_control(verbose = 0, REML = TRUE))
  expect_s3_class(fit, "mpmm")
  expect_equal(length(fit), 13)
})

#REML off
test_that("mpmm returns mpmm list w 13 elements", {
  fit <- mpmm(~ ice + (1 | id), data = ellie.ice.short, control = mpmm_control(verbose = 0, REML = FALSE))
  expect_s3_class(fit, "mpmm")
  expect_equal(length(fit), 13)
})

#optim, L-BFGS-B
test_that("mpmm returns mpmm list w 13 elements", {
  fit <- mpmm(~ ice + (1 | id), data = ellie.ice.short,
              control = mpmm_control(verbose = 0,
                                     optim = "optim",
                                     method = "L-BFGS-B"))
  expect_s3_class(fit, "mpmm")
  expect_equal(length(fit), 13)
})

#profile on
test_that("mpmm returns mpmm list w 13 elements", {
  fit <- mpmm(~ ice + (1 | id), data = ellie.ice.short, control = mpmm_control(verbose = 0, profile = TRUE))
  expect_s3_class(fit, "mpmm")
  expect_equal(length(fit), 13)
})

#fix rho
test_that("mpmm returns mpmm list w 13 elements", {
  fit <- mpmm(~ ice + (1 | id), data = ellie.ice.short,
              map = list(rho = factor(NA)),
              control = mpmm_control(verbose = 0, REML = TRUE))
  expect_s3_class(fit, "mpmm")
  expect_equal(length(fit), 13)
})
