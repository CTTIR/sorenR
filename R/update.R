#' Transactionally update a representation
#' @param object A fitted model, which is not modified in place.
#' @param newdata Dense numeric matrix with the original feature schema.
#' @param row_ids Optional identifiers, unique within this batch and distinct from
#'   currently retained identifiers. Unretained historical IDs cannot be checked.
#' @return A new model. Inspect the last history record for accepted, rollback or
#'   frozen status, actual displacement and budget.
#' @details Feasibility is checked on old retained anchors after every proposal.
#'   Growth proposals that violate the budget are rolled back in full. Same-size
#'   proposals are backtracked. A finite-anchor mean budget does not constrain
#'   every point, held-out data, or cumulative displacement over many updates.
#' @export
soren_update <- function(object, newdata, row_ids = NULL) {
  start <- proc.time()[[3L]]; .validate_model(object)
  x <- .input(newdata, object$schema)
  ids <- .ids(row_ids, nrow(x), object$n_seen)
  if (any(ids %in% object$memory$ids)) stop("row_ids overlap retained observations; supply distinct explicit row_ids.")
  z <- .preprocess(x, object$center, object$scale)
  ctl <- object$control; evidence <- .controller(object, z)
  old_target <- .map(object, object$memory$x); candidate <- object
  if (ctl$stability_mode != "frozen") {
  candidate$prototypes <- .move_prototypes(object$prototypes, z, ctl$learning_rate, ctl$clip)
  grew <- FALSE
  if (evidence$streak >= 2L && nrow(object$prototypes) < ctl$max_prototypes) {
    far <- which.max(evidence$residual)
    support <- sum(.dist2(z, z[far, , drop = FALSE]) <= object$bandwidth^2)
    if (evidence$residual[far] > ctl$growth_threshold * object$bandwidth && support >= ctl$growth_min_support) {
      candidate$prototypes <- rbind(candidate$prototypes, z[far, ])
      candidate$coordinates <- rbind(candidate$coordinates, .map(object, z[far, , drop = FALSE]))
      candidate$prototype_ids <- c(candidate$prototype_ids, object$next_prototype_id)
      candidate$next_prototype_id <- object$next_prototype_id + 1L
      grew <- TRUE
    }
  }
  candidate$bandwidth <- .bandwidth(candidate$prototypes)
  candidate$graph <- .graph(candidate$prototypes, candidate$bandwidth, ctl$neighbors)
  opt <- .optimize(candidate$coordinates, candidate$graph, ctl)
  candidate$coordinates <- opt$y
  } else {
    grew <- FALSE
    opt <- list(iterations = 0L)
  }
  displacement <- .displacement(candidate, object$memory$x, old_target)
  proposal_displacement <- displacement$mean
  budget <- evidence$budget; step <- 1; status <- "accepted"
  if (ctl$stability_mode == "frozen") {
    candidate <- object; status <- "frozen"; step <- 0
  } else if (displacement$mean > budget + ctl$tolerance) {
    full <- candidate; feasible <- FALSE
    if (!grew) for (attempt in seq_len(20L)) {
      step <- 2^(-attempt)
      candidate$prototypes <- object$prototypes + step * (full$prototypes - object$prototypes)
      candidate$coordinates <- object$coordinates + step * (full$coordinates - object$coordinates)
      candidate$bandwidth <- object$bandwidth + step * (full$bandwidth - object$bandwidth)
      candidate$graph <- .graph(candidate$prototypes, candidate$bandwidth, ctl$neighbors)
      displacement <- .displacement(candidate, object$memory$x, old_target)
      if (displacement$mean <= budget + ctl$tolerance) { feasible <- TRUE; break }
    }
    if (!feasible) { candidate <- object; status <- "rollback"; step <- 0 }
  }
  displacement <- .displacement(candidate, object$memory$x, old_target)
  if (status == "accepted") {
    local_seed <- as.integer((as.double(ctl$seed) + object$update_id + 1) %% .Machine$integer.max)
    candidate$memory <- .memory(object$memory, z, ids, ctl$memory_size, local_seed)
  }
  candidate$n_seen <- object$n_seen + nrow(x)
  candidate$update_id <- object$update_id + 1L
  candidate$controller <- evidence[c("ewma", "streak", "calibrated")]
  record <- list(update_id = candidate$update_id, status = status, budget = budget,
    displacement = displacement$mean, maximum_displacement = displacement$maximum,
    constraint_residual = if (is.finite(budget)) displacement$mean - budget else NA_real_,
    proposed_displacement = proposal_displacement, step = step,
    growth_proposed = grew, growth_accepted = grew && status == "accepted",
    support_score = evidence$score, calibrated = FALSE, iterations = opt$iterations,
    seconds = proc.time()[[3L]] - start)
  candidate$history <- utils::tail(c(object$history, list(record)), ctl$history_size)
  .validate_model(candidate)
  candidate
}
