#' Construct a tiltdens object
#' @keywords internal
#' @noRd
new_tiltdens <- function(x, bw, weights, n, from, to, data_name, call,
                         method, extra = list(), kernel = NULL) {
  if (is.null(from)) from <- min(x) - 3 * bw
  if (is.null(to))   to   <- max(x) + 3 * bw
  grid <- seq(from, to, length.out = n)

  obj <- c(
    list(
      x = grid,
      y = as.vector(kernel_matrix(grid, x, bw, kernel) %*% weights),
      bw = bw,
      n = length(x),
      call = call,
      data.name = data_name,
      has.na = FALSE,
      weights = weights,
      data = x,
      method = method
    ),
    extra
  )
  structure(obj, class = c("tiltdens", "density"))
}

#' Print a tilted or sharpened density fit
#'
#' @param x A `"tiltdens"` object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.tiltdens <- function(x, digits = getOption("digits") - 2L, ...) {
  label <- switch(
    x$method,
    tilt_density     = "Tilted density estimate (Doosti & Hall 2016)",
    tilt_density_cv  = "Tilted density estimate, cross-validated (Doosti, Hall & Mateu 2018)",
    sharpen_density  = "Data-sharpened density estimate (Doosti & Hall 2016)",
    "Perturbed density estimate"
  )
  cat("\n", label, "\n\n", sep = "")
  cat("Call:      "); print(x$call)
  cat("Data:      ", x$data.name, " (", x$n, " obs.)\n", sep = "")
  cat("Bandwidth: ", format(x$bw, digits = digits),
      if (!is.null(x$kernel)) paste0("  (", x$kernel$name, " kernel)") else "",
      "\n", sep = "")

  if (!is.null(x$comparator)) {
    cat("Comparator:", x$comparator, "with bandwidth",
        format(x$comparator_bw, digits = digits), "\n")
  }
  cat("Blocks:    ", x$n_groups,
      if (identical(x$n_groups, x$n)) " (one per observation)" else "", "\n", sep = "")
  if (length(x$breaks)) {
    cat("Breaks:    at ranks ", paste(x$breaks, collapse = ", "),
        "  (x = ", paste(format(x$break_values, digits = 3), collapse = ", "), ")",
        ", chosen by \"", x$breaks_method, "\"\n", sep = "")
  }

  if (!is.null(x$distance2)) {
    cat("Distance:  ", format(x$distance2, digits = digits),
        " (squared L2 to the comparator)\n", sep = "")
  }
  if (!is.null(x$cv)) {
    cat("CV:        ", format(x$cv, digits = digits), "\n", sep = "")
  }
  cat("Minimum:   ", format(min(x$y), digits = digits),
      "  (a proper density cannot go below zero)\n", sep = "")
  cat("\n")
  invisible(x)
}

#' Evaluate a fitted density at new points
#'
#' @param object A `"tiltdens"` object.
#' @param newdata Numeric vector of points at which to evaluate. Defaults to the
#'   grid stored in `object`.
#' @param ... Ignored.
#'
#' @return A numeric vector of density values.
#'
#' @examples
#' set.seed(1)
#' fit <- tilt_density(rnorm(60), m = 3)
#' predict(fit, newdata = c(-1, 0, 1))
#'
#' @export
predict.tiltdens <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$y)
  newdata <- as.vector(newdata)
  as.vector(kernel_matrix(newdata, object$data, object$bw,
                          object$kernel) %*% object$weights)
}

#' Integrated squared error against a known density
#'
#' Approximates \eqn{\int (\hat f - f)^2} by the trapezoidal rule over the grid
#' of the fit. Useful for simulation studies where the truth is known.
#'
#' @param object A `"tiltdens"` or `"density"` object.
#' @param true_density A function of one numeric vector.
#'
#' @return A single number.
#'
#' @examples
#' set.seed(1)
#' fit <- tilt_density(rnorm(80), m = 3)
#' ise(fit, dnorm)
#'
#' @export
ise <- function(object, true_density) {
  grid <- object$x
  diff_sq <- (object$y - true_density(grid))^2
  sum(diff(grid) * (diff_sq[-1L] + diff_sq[-length(diff_sq)]) / 2)
}
