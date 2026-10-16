#' Configure an experimental SOREN model
#'
#' All budgets refer to standardized input space or raw output coordinates.
#' Adaptive support scores are heuristics, not calibrated tests of drift.
#' @param max_prototypes Maximum codebook size (at least two).
#' @param initial_prototypes Initial codebook size, at most the maximum.
#' @param memory_size Maximum retained observations (at least two).
#' @param neighbors Maximum mutual nearest neighbors per prototype.
#' @param iterations Maximum layout descent iterations per fit or update.
#' @param learning_rate Prototype step in (0,1].
#' @param stability_mode One of fixed, adaptive, frozen, unrestricted.
#' @param stability_budget Mean squared displacement budget on old replay anchors.
#' @param tolerance Absolute numerical feasibility tolerance.
#' @param repulsion Nonedge loss multiplier, nonnegative.
#' @param clip Maximum individual standardized residual norm.
#' @param growth_threshold Growth residual threshold in bandwidth units.
#' @param growth_min_support Minimum observations supporting a growth proposal.
#' @param history_size Maximum retained update diagnostic records.
#' @param chunk_size Maximum rows per transform allocation.
#' @param seed Local random seed; caller random state is restored.
#' @return A validated list of class \code{soren_control}.
#' @export
soren_control <- function(max_prototypes = 32L, initial_prototypes = 8L,
                          memory_size = 128L, neighbors = 3L, iterations = 30L,
                          learning_rate = 0.2, stability_mode = "fixed",
                          stability_budget = 0.01, tolerance = 1e-10,
                          repulsion = 0.1, clip = 3, growth_threshold = 2,
                          growth_min_support = 3L, history_size = 32L,
                          chunk_size = 1024L, seed = 42L) {
  args <- as.list(environment())
  counts <- c("max_prototypes", "initial_prototypes", "memory_size", "neighbors",
              "iterations", "growth_min_support", "history_size", "chunk_size", "seed")
  for (nm in counts) {
    lower <- if (nm == "seed") 0 else if (nm %in% c("max_prototypes", "memory_size")) 2 else 1
    .scalar(args[[nm]], nm, lower, .Machine$integer.max, integer = TRUE)
  }
  if (initial_prototypes > max_prototypes) stop("initial_prototypes exceeds max_prototypes.")
  .scalar(learning_rate, "learning_rate", .Machine$double.eps, 1)
  for (nm in c("stability_budget", "tolerance", "repulsion")) .scalar(args[[nm]], nm, 0)
  for (nm in c("clip", "growth_threshold")) .scalar(args[[nm]], nm, .Machine$double.eps)
  if (!is.character(stability_mode) || length(stability_mode) != 1L || is.na(stability_mode) ||
      !stability_mode %in% c("fixed", "adaptive", "frozen", "unrestricted")) stop("Invalid stability_mode.")
  structure(args, class = "soren_control")
}

.scalar <- function(x, name, lower, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper ||
      (integer && x != floor(x))) stop(name, " must be a finite ",
        if (integer) "integer " else "number ", "in [", lower, ", ", upper, "].", call. = FALSE)
  invisible(x)
}

.local_seed <- function(seed, fun) {
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  fun()
}
