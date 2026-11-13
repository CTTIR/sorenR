#' Map observations with a fitted reference map
#' @param object A fitted \code{soren_model}.
#' @param newdata Dense numeric matrix with the exact training column schema.
#' @param diagnostics Return support diagnostics alongside coordinates?
#' @param reference \code{"current"} for the present map, or \code{"origin"} for
#'   the map as it was originally fitted.
#' @return A coordinate matrix, or a list with coordinates and a data frame of
#'   support diagnostics.
#' @details The map is deterministic and consumes no random numbers; no
#'   preprocessing is refitted. An updated model still contains everything
#'   needed to evaluate the originally fitted map, so \code{reference = "origin"}
#'   returns the coordinates the first model version would have returned.
#'
#'   Diagnostics hold, per observation: \code{score} (distance to the nearest
#'   prototype in bandwidth units), the conformal \code{p_value} of that score
#'   against the held-out calibration scores, \code{novel} (p-value at or below
#'   the calibration level), \code{new_weight} (share of map weight on prototypes
#'   added after the fit) and \code{origin_displacement} (distance between the
#'   current and the original coordinates). For every observation
#'   \code{origin_displacement} is at most
#'   \code{(1 - new_weight) * origin_shift_max} plus \code{new_weight} times the
#'   distance between the new-prototype part of the map and the original map.
#'   P-values are valid on average for observations exchangeable with the
#'   calibration rows; they are not posterior probabilities.
#' @examples
#' x <- cbind(a = sin(1:60), b = cos(1:60))
#' model <- soren_fit(x, control = soren_control(iterations = 5L))
#' newdata <- cbind(a = sin(61:65), b = cos(61:65))
#' soren_transform(model, newdata)
#' soren_transform(model, newdata, diagnostics = TRUE)$support
#' @export
soren_transform <- function(object, newdata, diagnostics = FALSE,
                            reference = c("current", "origin")) {
  .validate_model(object); reference <- match.arg(reference)
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics))
    stop("diagnostics must be TRUE or FALSE.")
  x <- .input(newdata, object$schema, allow_empty = TRUE)
  z <- .preprocess(x, object$center, object$scale)
  if (!diagnostics) return(.map(object, z, origin = reference == "origin"))
  support <- .support(object, z)
  list(coordinates = if (reference == "origin") support$origin else support$current,
       support = support$table)
}

.support <- function(object, z) {
  now <- .map(object, z, weights_new = TRUE); then <- .map(object, z, origin = TRUE)
  score <- .score(z, object$prototypes, object$bandwidth)
  p <- .pvalue(score, object$calibration$score)
  list(current = now$y, origin = then, table = data.frame(score = score, p_value = p,
    novel = p <= .level(object), new_weight = now$new_weight,
    origin_displacement = sqrt(rowSums((now$y - then)^2)), row.names = rownames(z)))
}
