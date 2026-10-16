.graph <- function(c, h, k) {
  m <- nrow(c); edges <- matrix(0, m, m)
  if (m < 2L) return(edges)
  d <- .dist2(c, c); diag(d) <- Inf
  neighbors <- matrix(FALSE, m, m)
  for (i in seq_len(m)) neighbors[i, order(d[i, ], seq_len(m))[seq_len(min(k, m - 1L))]] <- TRUE
  mutual <- neighbors & t(neighbors)
  edges[mutual] <- exp(-d[mutual] / h^2)
  edges
}
