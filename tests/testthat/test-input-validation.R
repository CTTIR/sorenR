test_that("fit, transform and diagnostics reject unusable arguments", {
  x <- cbind(a = sin(1:40), b = cos(1:40))
  expect_error(soren_fit(x[1, , drop = FALSE]), "two training observations")
  expect_error(soren_fit(x, control = list(seed = 1)), "complete soren_control")
  expect_error(soren_fit(x, n_components = 3L), "n_components")
  dup <- x; colnames(dup) <- c("a", "a")
  expect_error(soren_fit(dup), "nonempty and unique")
  model <- soren_fit(x, control = soren_control(iterations = 2L))
  expect_error(soren_transform(model, x, diagnostics = NA), "TRUE or FALSE")
  expect_error(soren_transform(model, x * 1e305), "Numeric range")
  expect_error(soren_embedding(model, type = "history"))
  expect_error(soren_load(c("a.rds", "b.rds")), "local path")
})

test_that("print and summary report the validated model", {
  x <- cbind(a = sin(1:40), b = cos(1:40))
  model <- soren_fit(x, control = soren_control(iterations = 2L))
  expect_output(out <- withVisible(print(model)),
                "Experimental SOREN model: [0-9]+ prototypes, 2 dimensions, update 0")
  expect_false(out$visible); expect_identical(out$value, model)
  expect_identical(summary(model), soren_diagnostics(model))
  broken <- model; broken$bandwidth[1] <- -1
  expect_error(print(broken), "Invalid or unsupported")
})
