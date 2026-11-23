safety_base <- function(n = 400, seed = 1) {
  set.seed(seed); centers <- rbind(c(0, 0, 0, 0), c(6, 0, 0, 0), c(0, 6, 0, 0))
  x <- centers[sample(3, n, TRUE), ] + matrix(rnorm(n * 4), n)
  colnames(x) <- paste0("f", 1:4); x
}
safety_novel <- function(n = 30, seed = 2) {
  set.seed(seed); x <- matrix(rnorm(n * 4, sd = .3), n); x[, 4] <- x[, 4] + 15
  colnames(x) <- paste0("f", 1:4); x
}
safety_control <- function(...) soren_control(max_prototypes = 40L, initial_prototypes = 20L,
                                             memory_size = 60L, iterations = 30L, ...)
frozen_fields <- c("prototypes", "bandwidth", "coordinates", "origin", "birth",
                   "prototype_ids", "next_prototype_id", "graph", "linear", "calibration")

test_that("a candidate that fails verification is rolled back in full", {
  x <- safety_base(); fit <- soren_fit(x, control = safety_control())
  batch <- safety_novel()
  expect_identical(soren_update(fit, batch)$history[[1]]$status, "accepted")
  real <- .certificate; calls <- 0L
  local_mocked_bindings(.certificate = function(candidate, object, radius) {
    calls <<- calls + 1L; out <- real(candidate, object, radius)
    if (calls == 1L) out$holds <- FALSE
    out
  })
  up <- soren_update(fit, batch)
  record <- up$history[[1]]
  expect_identical(calls, 2L)
  expect_identical(record$status, "rollback")
  expect_identical(record$growth, 0L); expect_gt(record$flagged, 0)
  expect_equal(record$layout_shift_max, 0); expect_equal(record$anchor_displacement_max, 0)
  for (field in c(frozen_fields, "memory")) expect_identical(up[[field]], fit[[field]], info = field)
  expect_equal(up$n_seen, fit$n_seen + nrow(batch)); expect_identical(up$update_id, 1L)
  expect_identical(soren_transform(up, x), soren_transform(fit, x))
  expect_identical(soren_diagnostics(up)$last_update$status, "rollback")
})

test_that("proposals that would intrude on existing support are rejected", {
  x <- safety_base(); fit <- soren_fit(x, control = safety_control())
  local_mocked_bindings(.new_weight = function(object, x, m) rep(1, max(1L, nrow(x))))
  up <- soren_update(fit, safety_novel())
  record <- up$history[[1]]
  expect_identical(record$status, "unchanged")
  expect_identical(record$growth, 0L); expect_gte(record$growth_rejected, 1L)
  for (field in frozen_fields) expect_identical(up[[field]], fit[[field]], info = field)
})

test_that("isolated novel observations do not create prototypes", {
  x <- safety_base(); fit <- soren_fit(x, control = safety_control())
  scattered <- rbind(c(0, 0, 0, 15), c(0, 0, 15, 0), c(-15, 0, 0, 0))
  colnames(scattered) <- colnames(x)
  up <- soren_update(fit, scattered)
  record <- up$history[[1]]
  expect_identical(record$flagged, 3L)
  expect_identical(record$growth, 0L); expect_identical(record$growth_rejected, 0L)
  expect_identical(record$status, "unchanged")
  expect_identical(up$prototypes, fit$prototypes)
})

test_that("explicit row identifiers are retained and must stay distinct", {
  x <- cbind(a = sin(1:60), b = cos(1:60))
  fit <- soren_fit(x, control = soren_control(iterations = 2L),
                   row_ids = sprintf("obs-%02d", 1:60))
  expect_true(all(rownames(soren_embedding(fit)) %in% sprintf("obs-%02d", 1:60)))
  expect_error(soren_update(fit, x[1:2, ], row_ids = c("obs-01", "new")), "overlap")
  up <- soren_update(fit, x[1:2, ], row_ids = c("new-1", "new-2"))
  expect_identical(up$history[[1]]$status, "unchanged")
  expect_true(all(c("new-1", "new-2") %in% up$memory$ids))
  expect_error(soren_update(fit, x[1:2, ], row_ids = "one"), "row_ids")
  expect_error(soren_fit(x, row_ids = rep(c("a", "b"), 30)), "row_ids")
})
