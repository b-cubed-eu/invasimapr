#' Standardize a site-by-species matrix by site (row-wise z) with optional robust SD
#'
#' @description
#' For each site (row), subtract the row-mean and divide by a row-scale.
#' The scale can be classical SD or a robust max(MAD, SD) with a small floor.
#'
#' @param X Numeric matrix (sites × species).
#' @param robust Logical. If TRUE, uses max(MAD, SD, 1e-8) per row; else classical SD with a safe 1 fallback.
#'
#' @return List with \code{z} (z-scored matrix), \code{mu} (row means), \code{sd} (row scales).
#' @examples
#' \dontrun{
#' std = standardise_by_site(C_js, robust = TRUE)
#' C_js_z = std$z
#' }
#' @importFrom matrixStats rowMeans2 rowSds
#' @export
standardise_by_site = function(X, robust = FALSE) {
  stopifnot(is.matrix(X))
  mu = matrixStats::rowMeans2(X, na.rm = TRUE)
  if (!robust) {
    sdv = matrixStats::rowSds(X, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv == 0] = 1
  } else {
    .robust_sd = function(x) {
      m = stats::mad(x, center = stats::median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE)
      s = stats::sd(x, na.rm = TRUE)
      max(m, s, 1e-8)
    }
    sdv = apply(X, 1, .robust_sd)
  }
  z = sweep(sweep(X, 1, mu, "-"), 1, sdv, "/")
  z[!is.finite(z)] = 0
  list(z = z, mu = mu, sd = sdv)
}
