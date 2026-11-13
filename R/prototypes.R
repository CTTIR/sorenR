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
      r <- z - rep(c[j, ], each = nrow(z))
      lengths <- sqrt(rowSums(r^2))
      r <- r * pmin(1, clip / pmax(lengths, .Machine$double.eps))
      c[j, ] <- c[j, ] + rate * colMeans(r)
    }
  }
  c
}

# One bandwidth per prototype: the distance to its second nearest prototype.
.bandwidths <- function(c, rank = 2L) {
  m <- nrow(c)
  if (m < 2L) return(rep(1, m))
  d <- .dist2(c, c); diag(d) <- Inf
  h <- sqrt(apply(d, 1L, function(r) sort(r)[min(rank, m - 1L)]))
  ok <- is.finite(h) & h > 0
  h[!ok] <- if (any(ok)) stats::median(h[ok]) else 1
  pmax(h, 1e-6)
}
