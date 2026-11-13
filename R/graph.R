# Symmetric nearest-neighbor graph with Euclidean edge lengths (zero = no edge).
# Components are joined by their shortest connecting edge so that every pair
# of prototypes has a finite path length.
.graph <- function(c, k) {
  m <- nrow(c); edges <- matrix(0, m, m)
  if (m < 2L) return(edges)
  d <- sqrt(.dist2(c, c)); adj <- matrix(FALSE, m, m)
  for (i in seq_len(m)) {
    o <- order(d[i, ], seq_len(m)); o <- o[o != i]
    adj[i, o[seq_len(min(k, m - 1L))]] <- TRUE
  }
  adj <- adj | t(adj)
  repeat {
    comp <- .component(adj, 1L)
    if (all(comp)) break
    cross <- d; cross[!comp, ] <- Inf; cross[, comp] <- Inf
    at <- which(cross == min(cross), arr.ind = TRUE)[1L, ]
    adj[at[1L], at[2L]] <- adj[at[2L], at[1L]] <- TRUE
  }
  edges[adj] <- pmax(d[adj], 1e-12)
  edges
}

.component <- function(adj, start) {
  seen <- rep(FALSE, nrow(adj)); seen[start] <- TRUE; frontier <- seen
  while (any(frontier)) {
    reach <- colSums(adj[frontier, , drop = FALSE]) > 0
    frontier <- reach & !seen; seen <- seen | frontier
  }
  seen
}

# Layout target distances between prototypes.
.target <- function(c, graph, geodesic) {
  m <- nrow(c)
  if (m < 2L || !geodesic) return(sqrt(.dist2(c, c)))
  d <- graph; d[graph == 0] <- Inf; diag(d) <- 0
  for (k in seq_len(m)) {
    via <- d[, k] + rep(d[k, ], each = m); shorter <- via < d
    d[shorter] <- via[shorter]
  }
  if (any(!is.finite(d))) stop("Numeric range exceeded in path lengths.")
  d
}

# Frozen local linear maps: for each prototype the minimum-norm ridge solution
# sending differences to graph neighbors onto their layout differences.
.linear_maps <- function(c, y, graph, rows, old = NULL) {
  m <- nrow(c); p <- ncol(c); d <- ncol(y)
  out <- lapply(seq_len(d), function(k) {
    a <- matrix(0, m, p); if (!is.null(old)) a[seq_len(nrow(old[[k]])), ] <- old[[k]]; a
  })
  for (j in rows) {
    nb <- which(graph[j, ] > 0)
    if (!length(nb)) next
    dc <- c[nb, , drop = FALSE] - rep(c[j, ], each = length(nb))
    dy <- y[nb, , drop = FALSE] - rep(y[j, ], each = length(nb))
    g <- tcrossprod(dc); lambda <- 1e-3 * mean(diag(g)) + 1e-12
    a <- crossprod(dc, solve(g + diag(lambda, length(nb)), dy))
    for (k in seq_len(d)) out[[k]][j, ] <- a[, k]
  }
  out
}
