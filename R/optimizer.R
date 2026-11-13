# Weighted stress with weights 1 / target^2, normalised so that the value is the
# mean squared relative distance error over prototype pairs.
.stress_setup <- function(target) {
  w <- 1 / pmax(target, 1e-12)^2; w[target <= 0] <- 0; diag(w) <- 0
  lower <- lower.tri(target); total <- sum(w[lower] * target[lower]^2)
  list(target = target, w = w, target_lower = target[lower], w_lower = w[lower],
       total = if (total > 0) total else 1)
}

# Pairs are taken in the order of stats::dist(), the lower triangle by column.
.stress <- function(y, s)
  sum(s$w_lower * (as.vector(stats::dist(y)) - s$target_lower)^2) / s$total

.classical <- function(target, d) {
  m <- nrow(target); y <- matrix(0, m, d)
  if (m < 2L) return(y)
  b <- -0.5 * target^2
  b <- b - rep(colMeans(b), each = m); b <- b - rowMeans(b)
  e <- eigen((b + t(b)) / 2, symmetric = TRUE)
  use <- seq_len(min(d, m - 1L)); lambda <- pmax(e$values[use], 0)
  v <- e$vectors[, use, drop = FALSE]
  # Eigenvector signs are arbitrary: fix them so the layout is reproducible.
  for (k in use) { i <- which.max(abs(v[, k])); if (v[i, k] < 0) v[, k] <- -v[, k] }
  y[, use] <- v * rep(sqrt(lambda), each = m)
  y + matrix(sin(seq_along(y)) * 1e-6 * max(target), m)
}

.project <- function(y, center, radius, rows) {
  if (!length(rows) || !is.finite(radius)) return(y)
  delta <- y[rows, , drop = FALSE] - center[rows, , drop = FALSE]
  len <- sqrt(rowSums(delta^2)); over <- len > radius
  if (any(over)) y[rows[over], ] <- center[rows[over], , drop = FALSE] +
    delta[over, , drop = FALSE] * (radius / len[over])
  y
}

# Stress majorization. Rows in `rows` stay inside a ball of `radius` around
# `center`; a step is a convex combination with a projected Guttman transform
# and is taken only if stress does not increase, so every iterate is feasible.
.layout <- function(y, target, iterations, center = NULL, radius = Inf, rows = integer()) {
  m <- nrow(y); names_in <- dimnames(y); s <- .stress_setup(target)
  if (m < 2L) return(list(y = y, initial = 0, stress = 0, iterations = 0L))
  v <- -s$w; diag(v) <- rowSums(s$w)
  vinv <- solve(v + 1 / m) - 1 / m
  current <- .stress(y, s); initial <- current; used <- 0L
  for (iter in seq_len(iterations)) {
    d <- sqrt(.dist2(y, y)); ratio <- s$w * s$target / d; ratio[d == 0] <- 0
    b <- -ratio; diag(b) <- rowSums(ratio)
    proposal <- vinv %*% (b %*% y)
    # Stress ignores translation: keep the constrained rows (or all rows) centered
    # where they are, so that no part of the radius is spent on a rigid shift.
    fix <- if (length(rows)) rows else seq_len(m)
    proposal <- proposal + rep(colMeans(y[fix, , drop = FALSE]) -
                               colMeans(proposal[fix, , drop = FALSE]), each = m)
    proposal <- .project(proposal, center, radius, rows)
    step <- 1; accepted <- FALSE
    for (trial in seq_len(20L)) {
      candidate <- y + step * (proposal - y); value <- .stress(candidate, s)
      if (is.finite(value) && value <= current) { accepted <- TRUE; break }
      step <- step / 2
    }
    if (!accepted) break
    gain <- current - value; y <- candidate; current <- value; used <- iter
    if (gain <= 1e-9 * max(current, 1e-12)) break
  }
  dimnames(y) <- names_in
  list(y = y, initial = initial, stress = current, iterations = used)
}
