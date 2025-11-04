# ======================================================================
# 0) PREPARE INPUTS
# ======================================================================

#' @title Prepare inputs (assemble & align core tables)
#' @description Thin wrapper over \code{assemble_matrices()} that
#'   stores results in an \code{invasimapr_fit} object.
#' @inheritParams assemble_matrices
#' @return \code{invasimapr_fit}
#' @export
prepare_inputs = function(..., make_plots = FALSE) {
  # assemble_matrices
  # try({
  #   source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/assemble_matrices.R")
  # }, silent = TRUE)
  #
  # if (!exists("assemble_matrices")) {
  #   stop("assemble_matrices() not found.")
  # }
  #
  # # utils_internal
  # try({
  #   source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/utils_internal.R")
  # }, silent = TRUE)
  #
  # if (!exists("new_invasimapr_fit")) {
  #   stop("new_invasimapr_fit() not found.")
  # }

  core = assemble_matrices(..., make_plots = make_plots)

  fit = new_invasimapr_fit(list(
    inputs = core,
    meta = list(
      n_sites     = core$n_sites,
      n_residents = core$n_residents,
      n_traits    = core$n_traits
    )
  ))
  fit
}
