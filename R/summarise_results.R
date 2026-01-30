#' Summarise invasiveness and invasibility (tables, maps, rankings)
#'
#' @description
#' Thin wrapper around \link{summarise_invasiveness_invasibility} that extracts
#' inputs from a \link{new_invasimapr_fit} container, forwards optional arguments,
#' and optionally overlays a boundary layer on the site map. This function is
#' intended as a reporting step after computing invasion fitness and/or
#' establishment probabilities.
#'
#' @param fit An object produced by the invasimapr workflow, containing
#'   `fitness` and/or `prob` components created by upstream steps.
#' @param boundary_sf Optional object of class `sf` providing a boundary layer
#'   to overlay on the site map.
#' @param boundary_params Named list of aesthetics passed to
#'   `ggplot2::geom_sf()`. Defaults to
#'   `list(inherit.aes = FALSE, fill = NA, color = "black", size = 0.3)`.
#' @param traits_inv Optional override of the invader trait table used for
#'   invader-ranked summaries. If `NULL`, traits are taken from the container.
#' @param ... Additional arguments forwarded to
#'   \link{summarise_invasiveness_invasibility}.
#'
#' @details
#' \strong{Inputs used from the container}
#' \itemize{
#'   \item `fit$fitness$lambda_is`: site-by-invader invasion fitness matrix (optional).
#'   \item `fit$prob$p_is`: site-by-invader establishment probability matrix (optional).
#'   \item `fit$inputs$site_df`: site metadata with coordinates `x` and `y` (optional).
#'   \item `fit$inputs$comm_res`: resident community matrix for context summaries (optional).
#'   \item invader trait tables stored in the container, used for invader rankings.
#' }
#'
#' At least one of `lambda_is` or `p_is` must be present; otherwise an error
#' is thrown.
#'
#' \strong{Behaviour}
#' \enumerate{
#'   \item Extract invasion fitness and/or establishment probability matrices
#'         from the container.
#'   \item Select site metadata and an effective invader trait table when available.
#'   \item Call \link{summarise_invasiveness_invasibility} to construct tidy
#'         summary tables and plots (site maps, rankings, heatmaps).
#'   \item If a boundary layer is supplied and plots are available, overlay it
#'         on the site map using `geom_sf()`.
#' }
#'
#' \strong{Output layout}
#' The returned container gains a `summary` element mirroring the structure
#' returned by \link{summarise_invasiveness_invasibility}, typically containing
#' summary tables, plots, and the arguments used.
#'
#' @return
#' The input `fit` object with an added or updated `fit$summary` list.
#' Invisibly returns `fit` to support piping.
#'
#' @seealso
#' \link{compute_invasion_fitness},
#' \link{compute_establishment_probability},
#' \link{summarise_invasiveness_invasibility},
#' \link{assemble_matrices}
#'
#' @examples
#' \dontrun{
#' fit <- prepare_inputs(
#'   sites = site_df,
#'   residents = resident_df,
#'   invaders = invader_df,
#'   traits = trait_df
#' )
#'
#' fit <- learn_sensitivities(fit)
#' fit <- predict_invaders(fit, traits_inv = invader_traits)
#'
#' fit <- summarise_results(fit)
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
