#' Map observations with a frozen representation
#' @param object A fitted \code{soren_model}.
#' @param newdata Dense numeric matrix with the exact training column schema.
#' @param diagnostics Return support diagnostics alongside coordinates?
#' @return A coordinate matrix, or a list with coordinates and support diagnostics.
#' @details The normalized radial basis map is deterministic and consumes no random
#'   numbers. It lies in the convex hull of prototype coordinates. Far-field
#'   behavior is not a reliable extrapolation. No preprocessing is refitted.
#' @export
soren_transform <- function(object, newdata, diagnostics = FALSE) {
  .validate_model(object)
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics))
    stop("diagnostics must be TRUE or FALSE.")
  x <- .input(newdata, object$schema, allow_empty = TRUE)
  z <- .preprocess(x, object$center, object$scale); coordinates <- .map(object, z)
  if (diagnostics) return(list(coordinates = coordinates, support = .support(object, z)))
  coordinates
}

.support <- function(object, z) {
  residual <- numeric(nrow(z))
  if (nrow(z)) for (start in seq.int(1L, nrow(z), by = object$control$chunk_size)) {
    idx <- start:min(nrow(z), start + object$control$chunk_size - 1L)
    residual[idx] <- sqrt(apply(.dist2(z[idx, , drop = FALSE], object$prototypes), 1L, min))
  }
  data.frame(residual = residual, support_ratio = residual / object$reference_residual,
             outside_reference = residual > object$reference_residual,
             row.names = rownames(z))
}
