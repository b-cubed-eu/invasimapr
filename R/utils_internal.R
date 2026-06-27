#' Internal utilities
#'
#' Small helper functions used inside the package. Not part of the public API.
#'
#' @name utils_internal
#' @keywords internal
#' @aliases new_invasimapr_fit .standardise_df .scale_like .row_z
NULL

new_invasimapr_fit = function(x = list()) {
  class(x) = unique(c("invasimapr_fit", class(x)))
  x
}

#' Internal infix helper
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' invasimapr_fit object
#'
#' A lightweight S3 container used throughout **invasimapr** to store assembled
#' inputs, fitted models, sensitivities, invasion fitness, and derived summaries.
#'
#' Objects of class `invasimapr_fit` are created by
#' \link{prepare_inputs} (via \link{new_invasimapr_fit}) and are progressively
#' enriched by downstream modelling and prediction functions.
#'
#' @details
#' The object typically contains the following components:
#' \describe{
#'   \item{inputs}{Core data structures returned by \link{assemble_matrices}.}
#'   \item{meta}{Basic metadata such as number of sites, residents, and traits.}
#'   \item{residents, traits, sensitivities, fitness, prob, summary}{Optional
#'     components added by downstream workflows.}
#' }
#'
#' @name invasimapr_fit
#' @rdname invasimapr_fit
#' @keywords internal
NULL


#' @title Compact print method for invasimapr_fit
#' @description Internal print method. Not part of the public API.
#' @param x An object of class `invasimapr_fit`.
#' @param ... Unused; included for compatibility with the `print` generic.
#' @return Invisibly returns `x`, the `invasimapr_fit` object. Called for its
#'   side effect of printing a compact summary to the console.
#' @keywords internal
#' @export
print.invasimapr_fit = function(x, ...) {
  cat("<invasimapr_fit>\n")
  stages = c("inputs","traits","crowding","residents","sensitivities",
             "invaders","fitness","prob","summary")
  present = stages[stages %in% names(x)]
  if (length(present)) cat(" stages:", paste(present, collapse = " -> "), "\n")
  if (!is.null(x$meta)) {
    cat(sprintf(" sites: %s | residents: %s | invaders: %s\n",
                x$meta$n_sites %||% NA_integer_,
                x$meta$n_residents %||% NA_integer_,
                x$meta$n_invaders %||% NA_integer_))
  }
  invisible(x)
}

# ------------------------------ small helpers ------------------------------

# column-wise z for data.frame (numeric columns only; factors/characters kept)
.standardise_df = function(df, ref_means = NULL, ref_sds = NULL) {
  stopifnot(is.data.frame(df))
  out = df
  num = vapply(df, is.numeric, logical(1))
  if (!any(num)) return(list(data = out, means = numeric(0), sds = numeric(0)))

  X = as.matrix(df[, num, drop = FALSE])
  if (is.null(ref_means) || is.null(ref_sds)) {
    m = colMeans(X, na.rm = TRUE)
    s = apply(X, 2, stats::sd, na.rm = TRUE)
  } else {
    m = ref_means[colnames(X)]
    s = ref_sds  [colnames(X)]
  }
  s[!is.finite(s) | s == 0] = 1
  Z = sweep(sweep(X, 2, m, "-"), 2, s, "/")

  out[, num] = Z
  list(data = out, means = m, sds = s)
}

# scale a df like a reference (means/sds taken from residents)
.scale_like = function(df, ref_means, ref_sds) {
  .standardise_df(df, ref_means, ref_sds)$data
}

# quick row-wise z (matrix)
.row_z = function(M, robust = FALSE) {
  std = standardise_by_site(as.matrix(M), robust = robust)
  std
}

