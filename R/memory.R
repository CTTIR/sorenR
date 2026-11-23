# Half of the slots are a uniform reservoir sample, filled with Algorithm R as
# given by Vitter (1985), ACM Trans. Math. Softw. 11(1), 37-57,
# doi:10.1145/3147.3165, Section 2, pp. 39-40; the other half cover the data by
# farthest-point selection.
.memory <- function(old, x, ids, size, seed) {
  cap <- floor(size / 2)
  if (is.null(old)) old <- list(reservoir = x[FALSE, , drop = FALSE], reservoir_ids = character(),
                               coverage = x[FALSE, , drop = FALSE], coverage_ids = character(), seen = 0)
  reservoir <- old$reservoir; rid <- old$reservoir_ids; seen <- old$seen
  result <- .local_seed(seed, function() {
    for (i in seq_len(nrow(x))) {
      seen <- seen + 1
      if (nrow(reservoir) < cap) { reservoir <- rbind(reservoir, x[i, ]); rid <- c(rid, ids[i]) }
      else {
        j <- floor(stats::runif(1) * seen) + 1
        if (j <= cap) { reservoir[j, ] <- x[i, ]; rid[j] <- ids[i] }
      }
    }
    list(reservoir = reservoir, reservoir_ids = rid, seen = seen)
  })
  pool <- rbind(old$coverage, x); pool_ids <- c(old$coverage_ids, ids)
  pick <- .cover_indices(pool, size - cap)
  result$coverage <- pool[pick, , drop = FALSE]; result$coverage_ids <- pool_ids[pick]
  result$x <- rbind(result$reservoir, result$coverage)
  result$ids <- c(result$reservoir_ids, result$coverage_ids)
  # Reservoir and coverage may share an observation; each retained slot is an anchor.
  result
}
