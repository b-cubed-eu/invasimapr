#' Summarise invasiveness / invasibility (tables, maps, rankings)
#'
#' Thin wrapper around
#' \href{https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.html}{`summarise_invasiveness_invasibility()`}
#' that pulls inputs from an [`invasimapr_fit`] container, forwards optional
#' arguments, and (optionally) overlays a boundary layer on the site map.
#' Designed as a drop-in reporting step after computing invasion fitness and/or
#' establishment probabilities.
#'
#' @section What it uses from `fit`:
#' - `fit$fitness$lambda_is` (optional): site × invader invasion fitness matrix.
#' - `fit$prob$p_is` (optional): site × invader establishment probability matrix.
#' - `fit$inputs$site_df` (optional): site metadata with `site`, `x`, `y`.
#' - `fit$inputs$comm_res` (optional): resident community matrix for context stats.
#' - `fit$invaders$traits_inv_raw` **or** `fit$invaders$traits_inv_glmm`
#'   (optional): invader traits used for invader-ranked summaries.
#'
#' At least one of `lambda_is` or `p_is` must be available; otherwise an error
#' is thrown.
#'
#' @param fit An object of class `invasimapr_fit` containing `fitness` and/or
#'   `prob` components produced by upstream steps (see **See also**).
#' @param boundary_sf Optional **sf** object to overlay on the site map
#'   (e.g., national/park boundary).
#' @param boundary_params Named list of aesthetics for `ggplot2::geom_sf()`.
#'   Defaults to `list(inherit.aes = FALSE, fill = NA, color = "black", size = 0.3)`.
#' @param traits_inv Optional override of the invader trait table used in
#'   summaries. If `NULL`, falls back to `fit$invaders$traits_inv_raw`, then
#'   `fit$invaders$traits_inv_glmm` when available.
#' @param ... Additional arguments forwarded to
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.html}{`summarise_invasiveness_invasibility()`}.
#'
#' @details
#' **Behaviour**
#' 1. Extracts `lambda_is` and/or `p_is` from `fit`; errors if neither is found.
#' 2. Chooses `site_df` from `fit$inputs$site_df` when present.
#' 3. Selects an effective trait table for invader-ranked summaries in the order:
#'    user-supplied `traits_inv` → `fit$invaders$traits_inv_raw` →
#'    `fit$invaders$traits_inv_glmm`.
#' 4. Calls
#'    \href{https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.html}{`summarise_invasiveness_invasibility()`}
#'    to construct tidy summaries and ggplots (site maps, rankings, heatmaps).
#' 5. If `boundary_sf` is supplied and plots are available, overlays the boundary
#'    with `geom_sf()` using `boundary_params`.
#'
#' **Output layout**
#' The returned `fit` gains a `summary` element mirroring the structure from
#' `summarise_invasiveness_invasibility()` (e.g., `tables`, `plots`, `args_used`).
#'
#' @return The input `fit` with an added/updated `fit$summary` list containing
#'   summary tables and plots. Invisibly returns `fit` to support piping.
#'
#' @seealso
#' - Fitness: \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.html}{`compute_invasion_fitness()`}
#' - Probabilities: \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.html}{`compute_establishment_probability()`}
#' - Summaries: \href{https://b-cubed-eu.github.io/invasimapr/reference/summarise_invasiveness_invasibility.html}{`summarise_invasiveness_invasibility()`}
#' - Inputs: \href{https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.html}{`assemble_matrices()`}
#'
#' @examples
#' \dontrun{
#' fit <- prepare_inputs(sites = site_df, residents = resident_df,
#'                       invaders = invader_df, traits = trait_df)
#' fit <- learn_sensitivities(fit)
#' fit <- predict_invaders(fit, traits_inv = invader_traits)
#'
#' # Optionally compute lambda and/or probability
#' fit$fitness <- compute_invasion_fitness_from_fit(fit)   # hypothetical helper
#' fit$prob    <- compute_establishment_probability_from_fit(fit)  # idem
#'
#' # Overlay a boundary (sf) on the site map
#' fit <- summarise_results(fit, boundary_sf = rsa_boundary)
#' fit$summary$plots$site_map
#' }
#'
#' @export
summarise_results = function(
    fit,
    boundary_sf = NULL,
    boundary_params = list(inherit.aes = FALSE, fill = NA, color = "black", size = 0.3),
    traits_inv = NULL,
    ...
){
  # --- Preconditions -----------------------------------------------------------
  `%||%` = function(a,b) if(!is.null(a)) a else b
  stopifnot(inherits(fit, "invasimapr_fit"))

  # Extract lambda and/or probability from the container (robust to try-errors)
  lambda  = if (!is.null(fit$fitness) && !"try-error" %in% class(fit$fitness))
    fit$fitness$lambda_is else NULL
  p_is    = if (!is.null(fit$prob)    && !"try-error" %in% class(fit$prob))
    fit$prob$p_is else NULL
  if (is.null(lambda) && is.null(p_is))
    stop("summarise_results(): need either fit$fitness$lambda_is or fit$prob$p_is.")

  # Optional site metadata (coordinates for mapping)
  site_df = fit$inputs$site_df %||% NULL

  # Choose effective traits for invader-ranked summaries:
  # user override -> raw traits saved at prediction time -> GLMM-scale traits
  traits_inv_eff = traits_inv %||%
    fit$invaders$traits_inv_raw %||%
    fit$invaders$traits_inv_glmm %||% NULL

  # Optional resident community matrix for contextual summaries
  comm_res = fit$inputs$comm_res %||% NULL

  # --- Delegate to the summary constructor ------------------------------------
  summ = summarise_invasiveness_invasibility(
    lambda_is  = lambda,
    p_is       = p_is,
    site_df    = site_df,
    traits_inv = traits_inv_eff,
    comm_res   = comm_res,
    ...
  )

  # --- Optional boundary overlay on site map ----------------------------------
  if (!is.null(boundary_sf) &&
      !is.null(summ$plots) &&
      !is.null(summ$plots$site_map) &&
      requireNamespace("ggplot2", quietly = TRUE)) {
    summ$plots$site_map = do.call(
      `+`,
      list(
        summ$plots$site_map,
        do.call(ggplot2::geom_sf, c(list(data = boundary_sf), boundary_params))
      )
    )
  }

  # Attach and return for chaining
  fit$summary = summ
  fit
}
