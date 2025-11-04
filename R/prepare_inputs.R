#' Prepare inputs (assemble & align core tables)
#'
#' Thin wrapper around [`assemble_matrices()`] that standardises and stores the
#' assembled inputs in an [`invasimapr_fit`] object for downstream modelling.
#' Use this to run the input-assembly step once and pass a single structured
#' object through subsequent pipelines.
#'
#' @inheritParams assemble_matrices
#' @param make_plots logical. If `TRUE`, propagate `make_plots` to
#'   [`assemble_matrices()`] so that diagnostic plots are produced during
#'   assembly. Defaults to `FALSE`.
#'
#' @details
#' **What this does**
#' 1. Calls [`assemble_matrices()`] to build core site × species/trait matrices,
#'    perform consistency checks, and harmonise keys and factor levels.
#' 2. Wraps the returned list (`core`) into a lightweight S3 container
#'    `invasimapr_fit` with a dedicated `inputs` slot and a small `meta` block
#'    (counts of sites, residents, traits) used by later steps.
#'
#' **Object layout**
#' ```
#' invasimapr_fit
#' ├─ inputs : <list>   # output from assemble_matrices()
#' └─ meta   : <list>   # n_sites, n_residents, n_traits
#' ```
#'
#' **Notes**
#' - No mutation of the `core` object occurs here; this function only packages
#'   and annotates it.
#' - If you need to inspect the assembled inputs, use `fit$inputs` after return.
#'
#' @return An object of class `invasimapr_fit` with components:
#' \describe{
#'   \item{`inputs`}{The full list returned by [`assemble_matrices()`].}
#'   \item{`meta`}{A list with `n_sites`, `n_residents`, `n_traits`.}
#' }
#'
#' @seealso
#' - Input assembly: \href{https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.html}{`assemble_matrices()`}
#' - Container constructor: [`new_invasimapr_fit()`]
#'
#' @examples
#' \dontrun{
#' # Minimal pattern; pass the same arguments you would to assemble_matrices()
#' fit <- prepare_inputs(
#'   sites      = site_df,
#'   residents  = resident_df,
#'   invaders   = invader_df,
#'   traits     = trait_df,
#'   make_plots = TRUE
#' )
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
