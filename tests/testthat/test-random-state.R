restore_random_state <- function(kind, had_seed, seed) {
  suppressWarnings(do.call(RNGkind, as.list(kind)))
  if (had_seed) {
    assign(".Random.seed", seed, envir = .GlobalEnv) # nolint: object_name_linter.
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

test_that("a legacy sample kind is restored without repeated warnings", {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(restore_random_state(old_kind, had_seed, old_seed), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
  set.seed(11)
  seed <- get(".Random.seed", envir = .GlobalEnv)
  kind <- RNGkind()
  x <- cbind(a = sin(1:60), b = cos(1:60))
  control <- soren_control(iterations = 2L)
  expect_no_warning(model <- soren_fit(x, control = control))
  expect_no_warning(soren_update(model, x[1:5, , drop = FALSE]))
  expect_identical(RNGkind(), kind)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), seed)
  suppressWarnings(RNGkind(sample.kind = "Rejection"))
  set.seed(12)
  again <- soren_fit(x, control = control)
  fields <- setdiff(names(model), "fit")
  expect_identical(again[fields], model[fields])
})

test_that("a missing random seed stays missing after fitting and updating", {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(restore_random_state(old_kind, had_seed, old_seed), add = TRUE)
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  x <- cbind(a = sin(1:60), b = cos(1:60))
  model <- soren_fit(x, control = soren_control(iterations = 2L))
  soren_update(model, x[1:5, , drop = FALSE])
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
