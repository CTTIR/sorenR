.controller <- function(object, x) {
  residual <- sqrt(apply(.dist2(x, object$prototypes), 1L, min))
  score <- stats::median(residual) / object$reference_residual
  if (!is.finite(score)) stop("Numeric range exceeded in support score.")
  ewma <- 0.8 * object$controller$ewma + 0.2 * score
  streak <- if (score > 1.5) min(1000000L, object$controller$streak + 1L) else 0L
  budget <- object$control$stability_budget
  if (object$control$stability_mode == "adaptive" && streak >= 2L)
    budget <- budget * (1 + min(4, max(0, ewma - 1)))
  if (!is.finite(budget)) stop("Numeric range exceeded in adaptive stability budget.")
  if (object$control$stability_mode == "unrestricted") budget <- Inf
  list(ewma = ewma, streak = streak, calibrated = FALSE, score = score,
       budget = budget, residual = residual)
}

.displacement <- function(candidate, anchors, target) {
  squared <- rowSums((.map(candidate, anchors) - target)^2)
  list(mean = mean(squared), maximum = max(squared))
}
