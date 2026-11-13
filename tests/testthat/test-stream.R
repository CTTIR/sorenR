stream_data <- function(n, shift = 0, seed = 1) {
  set.seed(seed); centers <- rbind(c(0, 0, 0, 0), c(6, 0, 0, 0), c(0, 6, 0, 0))
  x <- centers[sample(3, n, TRUE), ] + matrix(rnorm(n * 4), n); x[, 3] <- x[, 3] + shift
  colnames(x) <- paste0("f", 1:4); x
}
novel_data <- function(n, seed = 2) {
  set.seed(seed); x <- matrix(rnorm(n * 4, sd = .3), n); x[, 4] <- x[, 4] + 15
  colnames(x) <- paste0("f", 1:4); x
}
stream_control <- function(...) soren_control(max_prototypes = 40L, initial_prototypes = 20L,
                                             memory_size = 60L, iterations = 30L, ...)

test_that("batches inside calibrated support leave the map exactly unchanged", {
  x <- stream_data(400); fit <- soren_fit(x, control = stream_control())
  expect_true(soren_diagnostics(fit)$calibration$calibrated)
  model <- fit
  for (i in 1:5) model <- soren_update(model, stream_data(60, seed = 10 + i))
  statuses <- vapply(model$history, function(r) r$status, "")
  expect_true(all(statuses %in% c("unchanged", "accepted")))
  if (all(statuses == "unchanged")) {
    expect_identical(soren_transform(model, x), soren_transform(fit, x))
    expect_identical(model$coordinates, fit$coordinates)
  }
  expect_identical(soren_transform(model, x, reference = "origin"), soren_transform(fit, x))
  expect_equal(model$n_seen, 700); expect_lte(nrow(model$memory$x), 60L)
  expect_gt(model$memory$seen, 400)
})

test_that("a new population is appended and every old prototype stays within the radius", {
  x <- stream_data(400); fit <- soren_fit(x, control = stream_control())
  model <- soren_update(fit, rbind(stream_data(60, seed = 5), novel_data(30)))
  record <- soren_diagnostics(model)$last_update; cert <- soren_diagnostics(model)$certificate
  expect_identical(record$status, "accepted"); expect_gte(record$growth, 1)
  m0 <- nrow(fit$prototypes)
  expect_identical(model$prototypes[seq_len(m0), ], fit$prototypes)
  expect_identical(model$bandwidth[seq_len(m0)], fit$bandwidth)
  expect_identical(model$origin[seq_len(m0), ], fit$coordinates)
  shift <- sqrt(rowSums((model$coordinates[seq_len(m0), ] - fit$coordinates)^2))
  expect_true(all(shift <= cert$radius + 1e-9)); expect_true(cert$holds)
  expect_equal(cert$radius, 0.05 * fit$layout_scale)
  # The original map is reproduced exactly from the updated model.
  expect_identical(soren_transform(model, x, reference = "origin"), soren_transform(fit, x))
  # The stated bound holds for every observation, including unseen ones.
  unseen <- rbind(stream_data(200, seed = 99), novel_data(50, seed = 98), stream_data(50, shift = 4, seed = 97))
  s <- soren_transform(model, unseen, diagnostics = TRUE)
  origin <- soren_transform(model, unseen, reference = "origin")
  z <- sorenR:::.preprocess(unseen, model$center, model$scale)
  new_part <- model; keep <- model$birth > 0L
  w <- sorenR:::.weights(z, model$prototypes, model$bandwidth)
  local <- vapply(1:2, function(k) rowSums(w[, keep, drop = FALSE] * (matrix(model$coordinates[keep, k],
    nrow(z), sum(keep), byrow = TRUE) + (z %*% t(model$linear[[k]][keep, , drop = FALSE]) -
    matrix(rowSums(model$linear[[k]][keep, , drop = FALSE] * model$prototypes[keep, , drop = FALSE]),
           nrow(z), sum(keep), byrow = TRUE)))), numeric(nrow(z)))
  omega <- s$support$new_weight
  gap <- sqrt(rowSums((local - omega * origin)^2))
  bound <- (1 - omega) * cert$origin_shift_max + gap
  expect_true(all(s$support$origin_displacement <= bound + 1e-8))
  inside <- omega < 1e-6
  expect_gt(sum(inside), 100); expect_true(all(s$support$origin_displacement[inside] <= cert$radius + 1e-6))
  # New observations are covered afterwards; old support keeps negligible new weight.
  again <- soren_transform(model, novel_data(40, seed = 7), diagnostics = TRUE)$support
  expect_lt(mean(again$novel), 0.5); expect_gt(mean(again$new_weight), 0.9)
  expect_lte(record$anchor_new_weight_max, 0.05 + 1e-12)
})

test_that("frozen mode only scores, and the lifetime radius survives repeated growth", {
  x <- stream_data(400); frozen <- soren_fit(x, control = stream_control(stability_mode = "frozen"))
  up <- soren_update(frozen, novel_data(30))
  expect_identical(up$history[[1]]$status, "frozen"); expect_gt(up$history[[1]]$flagged, 0)
  expect_identical(up$prototypes, frozen$prototypes); expect_identical(up$memory, frozen$memory)
  model <- fit <- soren_fit(x, control = stream_control(stability_radius = 0.02))
  for (i in 1:6) {
    batch <- novel_data(30, seed = i); batch[, 1] <- batch[, 1] + 8 * i
    model <- soren_update(model, batch)
    cert <- soren_diagnostics(model)$certificate
    expect_true(cert$holds); expect_lte(cert$origin_shift_max, cert$radius + 1e-9)
  }
  expect_gt(nrow(model$prototypes), nrow(fit$prototypes)); expect_lte(nrow(model$prototypes), 40L)
  expect_identical(soren_transform(model, x, reference = "origin"), soren_transform(fit, x))
  expect_lte(length(model$history), 32L)
  free <- soren_update(soren_fit(x, control = stream_control(stability_mode = "unrestricted")), novel_data(30))
  expect_identical(soren_diagnostics(free)$certificate$radius, Inf)
})

test_that("too few calibration rows disable flagging instead of guessing", {
  x <- stream_data(40); fit <- soren_fit(x, control = stream_control())
  d <- soren_diagnostics(fit)$calibration
  expect_false(d$calibrated); expect_identical(d$level, 0)
  up <- soren_update(fit, novel_data(30))
  expect_identical(up$history[[1]]$status, "unchanged"); expect_identical(up$history[[1]]$flagged, 0L)
  none <- soren_fit(stream_data(400), control = stream_control(calibration_size = 0L))
  expect_equal(nrow(none$calibration$x), 0L)
  expect_false(any(soren_transform(none, novel_data(5), diagnostics = TRUE)$support$novel))
})
