# Exact unordered-pair objective. Smoothing is part of the loss, not a hidden clip.
.shape <- function(y, graph, repulsion, gradient = TRUE) {
  m <- nrow(y); grad <- matrix(0, m, ncol(y)); loss <- 0
  pairs <- which(upper.tri(graph), arr.ind = TRUE)
  if (!nrow(pairs)) return(list(loss = loss, gradient = grad))
  w <- graph[pairs]; pos <- sum(w); neg <- sum(w == 0)
  for (r in seq_len(nrow(pairs))) {
    i <- pairs[r, 1L]; j <- pairs[r, 2L]; delta <- y[i, ] - y[j, ]; s <- sum(delta^2)
    if (w[r] > 0) {
      a <- w[r] / pos; loss <- loss + a * log1p(s); coeff <- a / (1 + s)
    } else {
      a <- repulsion / neg
      loss <- loss + a * (log1p(s) - log(s + 1e-4))
      coeff <- a * (1 / (1 + s) - 1 / (s + 1e-4))
    }
    if (gradient) {
      g <- 2 * coeff * delta; grad[i, ] <- grad[i, ] + g; grad[j, ] <- grad[j, ] - g
    }
  }
  list(loss = loss, gradient = grad)
}
