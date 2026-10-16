test_that("plot labels can be customized without changing the model", {
  x <- cbind(a = sin(seq_len(20)), b = cos(seq_len(20)))
  model <- soren_fit(x, control = soren_control(iterations = 1L))
  before <- serialize(model, NULL)
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  device <- grDevices::dev.cur()
  on.exit({ grDevices::dev.off(device); unlink(path) }, add = TRUE)

  default <- withVisible(plot(model))
  expect_null(default$value)
  expect_false(default$visible)
  custom <- withVisible(plot(model, xlab = "Custom horizontal",
                             ylab = "Custom vertical", main = "Custom title",
                             col = "blue"))
  expect_null(custom$value)
  expect_false(custom$visible)
  expect_identical(serialize(model, NULL), before)

  one_dimension <- soren_fit(x, n_components = 1L,
                            control = soren_control(iterations = 1L))
  expect_error(plot(one_dimension), "at least two output dimensions")
})
