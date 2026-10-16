.optimize <- function(y, graph, control) {
  before <- .shape(y, graph, control$repulsion)
  initial <- before$loss; used <- 0L
  for (iter in seq_len(control$iterations)) {
    if (max(abs(before$gradient)) < 1e-8) break
    step <- 1; accepted <- FALSE
    for (trial in seq_len(20L)) {
      proposal <- y - step * before$gradient
      after <- .shape(proposal, graph, control$repulsion)
      if (is.finite(after$loss) && after$loss <= before$loss) {
        accepted <- TRUE; break
      }
      step <- step / 2
    }
    if (!accepted) break
    y <- proposal; before <- after; used <- iter
  }
  list(y = y, initial_loss = initial, loss = before$loss, iterations = used)
}
