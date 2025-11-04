# ======================================================================
# 6) SUMMARISE / REPORT  — DROP-IN: add boundary layer support
#   (keeps previous fixes: no lambda_long; read prob from $p_is)
# ======================================================================

#' @title Summarise invasiveness/invasibility
#' @description Thin wrapper over \code{summarise_invasiveness_invasibility()}.
#' @param fit object with \code{$fitness} and/or \code{$prob}.
#' @param boundary_sf An \strong{sf} object to overlay as a boundary on the site map.
#' @param boundary_params Named list of aesthetics for \code{ggplot2::geom_sf};
#'   defaults to \code{list(inherit.aes=FALSE, fill=NA, color="black", size=0.3)}.
#' @param ... forwarded to \code{summarise_invasiveness_invasibility()}.
#' @return updated \code{invasimapr_fit}
#' @export
summarise_results = function(
    fit,
    boundary_sf = NULL,
    boundary_params = list(inherit.aes = FALSE, fill = NA, color = "black", size = 0.3),
    traits_inv = NULL,   # =- allow user override
    ...
){
  # try({ source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/summarise_invasiveness_invasibility.R") }, silent = TRUE)
  # if (!exists("summarise_invasiveness_invasibility")) stop("summarise_invasiveness_invasibility() not found.")
  `%||%` = function(a,b) if(!is.null(a)) a else b
  stopifnot(inherits(fit, "invasimapr_fit"))

  lambda  = if (!is.null(fit$fitness) && !"try-error" %in% class(fit$fitness)) fit$fitness$lambda_is else NULL
  p_is    = if (!is.null(fit$prob)    && !"try-error" %in% class(fit$prob))    fit$prob$p_is         else NULL
  if (is.null(lambda) && is.null(p_is)) stop("summarise_results(): need either fit$fitness$lambda_is or fit$prob$p_is.")

  site_df = fit$inputs$site_df %||% NULL

  # Prefer RAW traits if we saved them; fall back to override or the standardized ones
  traits_inv_eff = traits_inv %||%
    fit$invaders$traits_inv_raw %||%
    fit$invaders$traits_inv_glmm %||% NULL

  comm_res = fit$inputs$comm_res %||% NULL

  summ = summarise_invasiveness_invasibility(
    lambda_is  = lambda,
    p_is       = p_is,
    site_df    = site_df,
    traits_inv = traits_inv_eff,  # =- pass effective traits
    comm_res   = comm_res,
    ...
  )

  # (boundary overlay unchanged)
  if (!is.null(boundary_sf) &&
      !is.null(summ$plots) &&
      !is.null(summ$plots$site_map) &&
      requireNamespace("ggplot2", quietly = TRUE)) {
    summ$plots$site_map = do.call(
      `+`,
      list(summ$plots$site_map,
           do.call(ggplot2::geom_sf, c(list(data = boundary_sf), boundary_params)))
    )
  }

  fit$summary = summ
  fit
}

