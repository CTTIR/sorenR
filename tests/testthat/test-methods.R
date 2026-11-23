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

test_that("the plot uses equal axis scales unless asked otherwise", {
  x <- cbind(a = 3 * sin(seq_len(30)), b = cos(seq_len(30)))
  model <- soren_fit(x, control = soren_control(iterations = 2L))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path, width = 7, height = 4)
  device <- grDevices::dev.cur()
  on.exit({ grDevices::dev.off(device); unlink(path) }, add = TRUE)
  units_per_inch <- function() {
    usr <- graphics::par("usr"); pin <- graphics::par("pin")
    c(diff(usr[1:2]) / pin[1], diff(usr[3:4]) / pin[2])
  }
  plot(model)
  scale <- units_per_inch()
  expect_equal(scale[1], scale[2], tolerance = 1e-6)
  plot(model, asp = NA)
  free <- units_per_inch()
  expect_gt(abs(free[1] / free[2] - 1), 0.05)
})
