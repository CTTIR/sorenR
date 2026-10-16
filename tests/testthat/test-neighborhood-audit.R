audit_fixture <- function() {
  x <- cbind(a = sin(1:20), b = cos(1:20), c = (1:20) / 10)
  soren_fit(x, control = soren_control(max_prototypes = 4L, initial_prototypes = 4L,
    memory_size = 8L, iterations = 2L))
}

# Independent finite enumeration, not the production fractional-count formula.
audit_choices <- function(d, k) {
  threshold <- sort(d)[k]
  closer <- which(d < threshold); boundary <- which(d == threshold)
  choices <- utils::combn(seq_along(boundary), k - length(closer), simplify = FALSE)
  lapply(choices, function(x) c(closer, boundary[x]))
}

audit_enumeration <- function(dx, dy, k) {
  x <- audit_choices(dx, k); y <- audit_choices(dy, k)
  mean(unlist(lapply(x, function(a) vapply(y, function(b)
    length(intersect(a, b)) / k, numeric(1)))))
}

test_that("fractional tie overlap agrees with independent exhaustive choices", {
  cases <- list(
    list(c(1, 2, 3), c(1, 2, 3), 1, 1),
    list(c(1, 2, 3), c(3, 2, 1), 1, 0),
    list(c(1, 2, 2), c(1, 2, 2), 2, .75),
    list(c(1, 2, 2), c(1, 2, 3), 2, .75),
    list(rep(1, 4), rep(2, 4), 2, .5),
    list(1:5, rep(2, 5), 2, .4),
    list(1:5, 5:1, 5, 1),
    list(1, 3, 1, 1),
    list(c(0, 0, 1, 1, 1, 1, 2), c(0, 1, 0, 1, 1, 2, 1), 4, .625))
  for (case in cases) {
    dx <- matrix(case[[1]], 1); dy <- matrix(case[[2]], 1); k <- case[[3]]
    counts <- sorenR:::.audit_counts(dx, dy, sort(dx)[k], sort(dy)[k])
    got <- sorenR:::.audit_overlap(counts, k)$recall
    expect_equal(got, case[[4]])
    expect_equal(got, audit_enumeration(as.vector(dx), as.vector(dy), k))
    for (b in c(1, 2, 3, 20)) {
      tx <- ty <- matrix(Inf, 1, k)
      blocks <- split(seq_along(dx), ceiling(seq_along(dx) / b))
      for (idx in blocks) {
        tx <- sorenR:::.audit_topk(tx, dx[, idx, drop = FALSE], k)
        ty <- sorenR:::.audit_topk(ty, dy[, idx, drop = FALSE], k)
      }
      streamed <- matrix(0, 1, 8)
      for (idx in blocks) streamed <- streamed + sorenR:::.audit_counts(
        dx[, idx, drop = FALSE], dy[, idx, drop = FALSE], tx[, k], ty[, k])
      expect_identical(streamed, counts)
      expect_identical(sorenR:::.audit_overlap(streamed, k)$recall, got)
    }
  }
  mixed <- sorenR:::.audit_counts(matrix(c(0, 0, 1, 1, 1, 1, 2), 1),
    matrix(c(0, 1, 0, 1, 1, 2, 1), 1), 1, 1)
  expect_equal(as.vector(mixed), c(2, 4, 2, 4, 1, 1, 1, 2))
})

test_that("audit is exact across blocks and agrees with a full-bank oracle", {
  model <- audit_fixture()
  q <- cbind(a = c(0, .2, -.2), b = c(.5, -.5, .5), c = c(.3, .6, .3))
  r <- rbind(q, q[1, , drop = FALSE], q[2, , drop = FALSE], -q)
  exclude <- c(1L, 0L, 3L); k <- 4L
  expected <- soren_neighborhood_audit(model, q, r, k, exclude, 1L)
  for (b in c(2L, 3L, 20L))
    expect_identical(soren_neighborhood_audit(model, q, r, k, exclude, b), expected)
  zq <- sweep(sweep(q, 2, model$center), 2, model$scale, "/")
  zr <- sweep(sweep(r, 2, model$center), 2, model$scale, "/")
  yq <- do.call(rbind, lapply(seq_len(nrow(q)), function(i)
    soren_transform(model, q[i, , drop = FALSE])))
  yr <- do.call(rbind, lapply(seq_len(nrow(r)), function(i)
    soren_transform(model, r[i, , drop = FALSE])))
  expect_identical(unname(sorenR:::.audit_block(model, q)$embedding), unname(yq))
  for (i in seq_len(nrow(q))) {
    dx <- vapply(seq_len(nrow(r)), function(j) sum((zq[i, ] - zr[j, ])^2), numeric(1))
    dy <- vapply(seq_len(nrow(r)), function(j) sum((yq[i, ] - yr[j, ])^2), numeric(1))
    if (exclude[i]) { dx[exclude[i]] <- Inf; dy[exclude[i]] <- Inf }
    expect_equal(expected$per_query$recall_at_k[i], audit_enumeration(dx, dy, k))
  }
  perm <- c(4L, 8L, 2L, 5L, 1L, 3L, 7L, 6L)
  reindex <- ifelse(exclude == 0L, 0L, match(exclude, perm))
  reordered <- soren_neighborhood_audit(model, q, r[perm, ], k, reindex, 2L)
  cols <- setdiff(names(expected$per_query), "excluded_reference")
  expect_identical(reordered$per_query[cols], expected$per_query[cols])
  qp <- c(3, 1, 2)
  reordered_q <- soren_neighborhood_audit(model, q[qp, ], r, k, exclude[qp], 2L)
  expect_identical(reordered_q$per_query$recall_at_k, expected$per_query$recall_at_k[qp])
})

test_that("collapse, duplicates and explicit exclusions have honest denominators", {
  model <- audit_fixture(); model$coordinates[,] <- 0
  q <- matrix(0, 2, 3, dimnames = list(c("same", "same"), model$schema$names))
  r <- q[rep(1, 5), , drop = FALSE]
  a <- soren_neighborhood_audit(model, q, r, 2L, c(0L, 2L), 2L)
  expect_equal(a$per_query$recall_at_k, c(2/5, 2/4))
  expect_equal(a$per_query$eligible_reference_count, c(5, 4))
  expect_equal(a$per_query$input_boundary_size, c(5, 4))
  expect_equal(a$per_query$embedding_boundary_weight, c(2/5, 2/4))
  expect_identical(a$per_query$query_index, 1:2)
  expect_equal(a$summary$input_tied_query_fraction, 1)
  expect_false(a$calibrated)
  self <- soren_neighborhood_audit(model, q, q, 1L, seq_len(nrow(q)))
  expect_equal(self$per_query$recall_at_k, c(1, 1))
  expect_equal(self$per_query$input_boundary_size, c(1, 1))
})

test_that("empty and invalid audit arguments respect the existing schema", {
  model <- audit_fixture(); x <- model$memory$x
  colnames(x) <- model$schema$names
  empty <- soren_neighborhood_audit(model, x[FALSE, ], x, 2L, integer())
  expect_identical(nrow(empty$per_query), 0L)
  expect_identical(empty$per_query$query_index, integer())
  expect_identical(empty$per_query$recall_at_k, numeric())
  expect_identical(empty$summary$mean_recall_at_k, NA_real_)
  expect_identical(empty$summary$input_tied_query_fraction, NA_real_)
  expect_error(soren_neighborhood_audit(model, x[FALSE, ], x, nrow(x) + 1), "eligible")
  expect_error(soren_neighborhood_audit(model, x, x[FALSE, ], 1), "finite")
  for (k in list(0, -1, 1.5, NA, Inf, 1+1i, TRUE, c(1, 2)))
    expect_error(soren_neighborhood_audit(model, x, x, k))
  for (b in list(0, 1.5, NA, Inf, 1+1i))
    expect_error(soren_neighborhood_audit(model, x, x, 1, chunk_size = b))
  for (ex in list(1, rep(NA_real_, nrow(x)), rep(.5, nrow(x)), rep(-1, nrow(x)),
                  rep(Inf, nrow(x)), rep(nrow(x) + 1, nrow(x)), matrix(0, nrow(x), 1)))
    expect_error(soren_neighborhood_audit(model, x, x, 1, ex), "exclude_reference")
  expect_error(soren_neighborhood_audit(model, x, x, nrow(x), seq_len(nrow(x))), "eligible")
  expect_error(soren_neighborhood_audit(model, x[, 3:1], x, 1), "schema")
  expect_error(soren_neighborhood_audit(model, x, unname(x), 1), "schema")
  expect_error(soren_neighborhood_audit(model, as.data.frame(x), x, 1), "dense numeric")
  bad <- x; bad[1, 1] <- Inf
  expect_error(soren_neighborhood_audit(model, bad, x, 1), "finite")
  xi <- x; storage.mode(xi) <- "integer"
  expect_equal(soren_neighborhood_audit(model, xi, xi, 1)$per_query$recall_at_k,
    soren_neighborhood_audit(model, xi * 1, xi * 1, 1)$per_query$recall_at_k)
})

test_that("audit preserves model, inputs, RNG and process settings on success and error", {
  model <- audit_fixture(); x <- model$memory$x; colnames(x) <- model$schema$names
  model_before <- serialize(model, NULL); input_before <- serialize(x, NULL)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(815); seed <- .Random.seed; kind <- RNGkind(); opts <- options(); wd <- getwd()
  a <- soren_neighborhood_audit(model, x, x, 2)
  expect_error(soren_neighborhood_audit(model, x, x, nrow(x) + 1))
  expect_identical(.Random.seed, seed); expect_identical(RNGkind(), kind)
  expect_identical(options(), opts); expect_identical(getwd(), wd)
  rm(".Random.seed", envir = .GlobalEnv)
  soren_neighborhood_audit(model, x, x, 2)
  expect_error(soren_neighborhood_audit(model, x, x, 0))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), kind)
  expect_identical(serialize(model, NULL), model_before)
  expect_identical(serialize(x, NULL), input_before)
  expect_identical(names(a), c("per_query", "summary", "reference", "tie_rule", "calibrated"))
  expect_true(all(vapply(a$per_query, is.numeric, logical(1))))
  expect_true(all(lengths(a$per_query) == nrow(x)))
})

test_that("audit uses bounded blocks", {
  model <- audit_fixture(); x <- model$memory$x; colnames(x) <- model$schema$names
  block <- sorenR:::.audit_block; distances <- sorenR:::.audit_distances
  topk <- sorenR:::.audit_topk; seen <- list()
  local_mocked_bindings(
    .audit_block = function(object, x) {
      expect_lte(nrow(x), 2); seen[[length(seen) + 1L]] <<- nrow(x); block(object, x)
    },
    .audit_distances = function(query, reference, excluded, reference_indices) {
      expect_lte(nrow(query), 2); expect_lte(nrow(reference), 2)
      distances(query, reference, excluded, reference_indices)
    },
    .audit_topk = function(previous, distances, k) {
      expect_lte(nrow(previous), 2); expect_identical(ncol(previous), as.integer(k))
      expect_lte(ncol(distances), 2); topk(previous, distances, k)
    }, .package = "sorenR")
  soren_neighborhood_audit(model, x[1:5, ], x, 3, chunk_size = 2)
  expect_true(length(seen) > 2)
})

test_that("audit detects intermediate range failures without norm cancellation", {
  model <- audit_fixture(); x <- model$memory$x; colnames(x) <- model$schema$names
  huge <- x[1, , drop = FALSE]; huge[,] <- 1e308
  expect_error(soren_neighborhood_audit(model, huge, x, 1, chunk_size = 2), "Numeric range")
  a <- matrix(c(1e100, -1e100), 1); b <- a + 1e85
  expect_equal(sorenR:::.dist2(a, b)[1, 1], sum((a - b)^2))
})
