test_that("diagnostics return the documented fields", {
  x <- cbind(a = sin(1:250), b = cos(1:250))
  model <- soren_fit(x, control = soren_control(iterations = 2L))
  d <- soren_diagnostics(model)
  expect_identical(names(d), c("prototypes", "prototypes_added", "retained_slots",
    "retained_unique_ids", "seen", "reservoir_population", "edges", "stress",
    "reference", "certificate", "calibration", "last_update", "support"))
  expect_identical(names(d$certificate),
                   c("mode", "radius", "layout_scale", "origin_shift_max", "holds", "updates"))
  expect_identical(names(d$calibration),
                   c("n", "resolution", "alpha", "level", "calibrated"))
  expect_null(d$last_update); expect_null(d$support)
  # 50 held-out rows: 1 / 51 exceeds novelty_alpha = 0.01 and is the level used.
  expect_identical(d$calibration$n, 50L)
  expect_identical(d$calibration$level, 1 / 51)
})

test_that("frozen batches count as seen but not as reservoir population", {
  x <- cbind(a = sin(1:60), b = cos(1:60))
  model <- soren_fit(x, control = soren_control(iterations = 2L, stability_mode = "frozen"))
  up <- soren_update(model, x[1:7, , drop = FALSE])
  d <- soren_diagnostics(up)
  expect_identical(d$last_update$status, "frozen")
  expect_identical(d$seen, 67L); expect_identical(d$reservoir_population, 60)
})
