#' Prepare inputs (assemble and align core tables)
#'
#' @description
#' Thin wrapper around \link{assemble_matrices} that standardises and stores
#' assembled inputs in a \link{new_invasimapr_fit} object for downstream modelling.
#' Use this function to run the input-assembly step once and pass a single
#' structured container through subsequent pipelines.
#'
#' @inheritDotParams assemble_matrices
#' @param ... Arguments passed directly to \link{assemble_matrices}.
#' @param make_plots Logical; if `TRUE`, propagate `make_plots` to
#'   \link{assemble_matrices} so that diagnostic plots are produced during
#'   assembly. Defaults to `FALSE`.
#'
#' @details
#' \strong{What this does}
#' \enumerate{
#'   \item Calls \link{assemble_matrices} to build core site-by-species and
#'         site-by-trait matrices, perform consistency checks, and harmonise
#'         keys and factor levels.
#'   \item Wraps the returned list into a lightweight S3 container of class
#'         `invasimapr_fit`, with a dedicated `inputs` slot and a small `meta`
#'         block containing basic counts used by later steps.
#' }
#'
#' \strong{Object layout}
#' \preformatted{
#' invasimapr_fit
#'   inputs : list   (output from assemble_matrices)
#'   meta   : list   (n_sites, n_residents, n_traits)
#' }
#'
#' \strong{Notes}
#' \itemize{
#'   \item No mutation of the assembled inputs occurs here; the function only
#'         packages and annotates them.
#'   \item To inspect the assembled data, use `fit$inputs` on the returned object.
#' }
#'
#' @return
#' An object of class `invasimapr_fit` with components:
#' \describe{
#'   \item{inputs}{The full list returned by \link{assemble_matrices}.}
#'   \item{meta}{A list with elements `n_sites`, `n_residents`, and `n_traits`.}
#' }
#'
#' @seealso
#' \link{assemble_matrices},
#' \link{new_invasimapr_fit}
#'
#' @examples
#' \dontrun{
#' fit <- prepare_inputs(
#'   sites      = site_df,
#'   residents  = resident_df,
#'   invaders   = invader_df,
#'   traits     = trait_df,
#'   make_plots = TRUE
#' )
#'
#' str(fit, 1)
#' str(fit$inputs, 1)
#' fit$meta
#' }
#'
#' @export
prepare_inputs = function(..., make_plots = FALSE) {
  # -- Step 1: assemble the core matrices/tables -------------------------------
  # Delegates all validation, alignment, and key harmonisation to
  # assemble_matrices(). Any plotting is controlled via make_plots.
  core = assemble_matrices(..., make_plots = make_plots)

  # -- Step 2: wrap results in a lightweight S3 container ---------------------
  # The meta block carries minimal counts used by downstream summaries/guards.
  fit = new_invasimapr_fit(list(
    inputs = core,
    meta   = list(
      n_sites     = core$n_sites,
      n_residents = core$n_residents,
      n_traits    = core$n_traits
    )
  ))

  # -- Return: a ready-to-pass container for subsequent modelling -------------
  fit
}
