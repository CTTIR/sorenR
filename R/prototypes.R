.cover_indices <- function(x, size, first = 1L) {
  size <- min(size, nrow(x))
  selected <- first
  nearest <- as.vector(.dist2(x, x[first, , drop = FALSE]))
  while (length(selected) < size && max(nearest) > 0) {
    j <- which.max(nearest)
    selected <- c(selected, j)
    nearest <- pmin(nearest, as.vector(.dist2(x, x[j, , drop = FALSE])))
  }
  selected
}

.move_prototypes <- function(c, x, rate, clip) {
  assignment <- max.col(-.dist2(x, c), ties.method = "first")
  for (j in seq_len(nrow(c))) {
    z <- x[assignment == j, , drop = FALSE]
    if (nrow(z)) {
      r <- sweep(z, 2L, c[j, ], "-")
      lengths <- sqrt(rowSums(r^2))
      r <- r * pmin(1, clip / pmax(lengths, .Machine$double.eps))
      c[j, ] <- c[j, ] + rate * colMeans(r)
    }
  }
  c
}

.bandwidth <- function(c) {
  if (nrow(c) < 2L) return(1)
  d <- .dist2(c, c); diag(d) <- Inf
  nearest <- sqrt(apply(d, 1L, min)); nearest <- nearest[nearest > 0 & is.finite(nearest)]
  if (!length(nearest)) 1 else max(stats::median(nearest), 1e-6)
}
