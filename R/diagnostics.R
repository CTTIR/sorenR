#' Inspect a fitted representation
#' @param object A fitted model.
#' @param newdata Optional observations for support diagnostics.
#' @return A list of graph counts, retained memory, raw finite-anchor displacement,
#'   support information and prototype distance distortion. Scores are uncalibrated.
#' @details Prototype stress compares all pairwise standardized input distances
#'   with raw output distances; it is not a biological ground truth or an uncertainty
#'   probability. Input sensitivity and refit stability require separate experiments.
#' @export
soren_diagnostics <- function(object, newdata = NULL) {
  .validate_model(object)
  dx <- as.numeric(stats::dist(object$prototypes))
  dy <- as.numeric(stats::dist(object$coordinates))
  stress <- if (length(dx) && sum(dx^2) > 0) sqrt(sum((dx - dy)^2) / sum(dx^2)) else NA_real_
  list(prototypes = nrow(object$prototypes), retained_slots = nrow(object$memory$x),
    retained_unique_ids = length(unique(object$memory$ids)), seen = object$n_seen,
    reservoir_population = object$memory$seen,
    edges = sum(object$graph > 0) / 2, isolated_prototypes = sum(rowSums(object$graph) == 0),
    prototype_distance_stress = stress, reference = "frozen standardized Euclidean input",
    calibrated = FALSE, controller = object$controller,
    last_update = if (length(object$history)) utils::tail(object$history, 1L)[[1L]] else NULL,
    support = if (!is.null(newdata)) soren_transform(object, newdata, TRUE)$support else NULL)
}

#' Retrieve the coordinates of retained items
#' @param object A fitted model.
#' @param type Either replay slots or prototypes.
#' @return A coordinate matrix with identifiers in row names. Replay IDs may occur
#'   twice when reservoir and coverage slots retain the same observation.
#' @export
soren_embedding <- function(object, type = c("replay", "prototypes")) {
  .validate_model(object); type <- match.arg(type)
  if (type == "prototypes") {
    out <- object$coordinates; rownames(out) <- paste0("prototype-", object$prototype_ids)
  } else {
    out <- .map(object, object$memory$x); rownames(out) <- object$memory$ids
  }
  out
}
