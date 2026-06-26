#' Compute invasion fitness (\eqn{\lambda}) and optional establishment probability (\eqn{P})
#'
#' Wraps
#' \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.html}{`compute_invasion_fitness()`}
#' to generate invader-by-site invasion fitness \eqn{\lambda_{is}} under model
#' Options **A-E**, and (optionally) maps \eqn{\lambda} to probabilities via
#' \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.html}{`compute_establishment_probability()`}.
#' Supports calibration of \eqn{\kappa} so that **resident** mean \eqn{\lambda}
#' is approximately zero, and can overlay an `sf` boundary on map-like plots.
#'
#' @param fit An `invasimapr_fit` from
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.html}{`predict_invaders()`},
#'   containing `invaders$r_is_z/C_is_z/S_is_z`, resident summaries, and
#'   learned sensitivities.
#' @param option Character, one of `c("A","B","C","D","E")` selecting the
#'   invasion-fitness specification:
#'   \describe{
#'     \item{A}{Baseline: \eqn{\lambda = r^{(z)} - \alpha C^{(z)} - \beta S^{(z)}}.}
#'     \item{B}{Global abiotic scaling \eqn{\theta_0 \cdot r^{(z)}}.}
#'     \item{C}{Trait-varying abiotic scaling \eqn{\theta_i \cdot r^{(z)}}.}
#'     \item{D}{Site-varying abiotic and crowding \eqn{\Gamma_{is}, \alpha_{is}}.}
#'     \item{E}{Signed saturation effect (facilitation allowed via signed \eqn{\beta_i}).}
#'   }
#'
#' @param calibrate_kappa Logical; if `TRUE`, set \eqn{\kappa} so mean **resident**
#'   \eqn{\lambda \approx 0} (scale alignment for communication/comparison).
#' @param prob_method (legacy) `NULL` or one of `c("probit","logit","hard")`;
#'   preserved for backward compatibility.
#' @param prob_args (legacy) List of arguments passed to
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.html}{`compute_establishment_probability()`}
#'   (e.g., `sigma`, `tau`, `predictive`, `sigma_mat`, `site_df`, `return_long`, `make_plots`).
#' @param method Alias of `prob_method` (preferred in user code).
#' @param prob_scale Alias of `prob_args` (preferred in user code).
#' @param boundary_sf Optional **sf** object to overlay on map-like probability plots.
#' @param boundary_params Named list for `ggplot2::geom_sf()` aesthetics; default
#'   `list(inherit.aes=FALSE, fill=NA, color="black", size=0.3)`.
#'
#' @details
#' **What this wrapper does**
#' 1. Chooses the appropriate sensitivity inputs for the requested `option`
#'    (e.g., `theta_i` for C, `Gamma_is`/`alpha_is` for D).
#' 2. Calls ... with site-standardised inputs \eqn{(r^{(z)}, C^{(z)}, S^{(z)})} held in `fit`.
#' 3. (Optional) Converts \eqn{\lambda} to probability \eqn{P} via probit, logit,
#'    or a hard threshold, forwarding `prob_*` settings.
#' 4. (Optional) If `boundary_sf` is supplied and plots are available, overlays
#'    the boundary on map-like plots using `geom_sf()`.
#'
#' **Calibration (`calibrate_kappa = TRUE`)**
#' Aligns invader-resident scales by shifting \eqn{\lambda} so that the **resident**
#' mean is ~0, using resident moments and trait-plane slopes stored in `fit`.
#'
#' @return The input `fit` with:
#' \describe{
#'   \item{`$fitness`}{List returned by `compute_invasion_fitness()` (including
#'         `lambda_is`, `lambda_long`, `option`, and any plots/maps).}
#'   \item{`$prob`}{If requested, list returned by `compute_establishment_probability()` with
#'         probability tables and plots (with boundary overlay applied when provided).}
#' }
#'
#' @seealso
#' - Fitness: \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_invasion_fitness.html}{`compute_invasion_fitness()`}
#' - Probabilities: \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_establishment_probability.html}{`compute_establishment_probability()`}
#' - Sensitivities: \href{https://b-cubed-eu.github.io/invasimapr/reference/learn_sensitivities.html}{`learn_sensitivities()`}
#' - Invader predictors: \href{https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.html}{`predict_invaders()`}
#'
#' @examples
#' \dontrun{
#' fit = fit |>
#'   predict_invaders(traits_inv) |>
#'   predict_establishment(option = "C", calibrate_kappa = TRUE,
#'                         method = "probit", prob_scale = list(sigma = 1))
#'
#' # Overlay a boundary on probability maps
#' fit = predict_establishment(fit, option = "D", method = "logit",
#'                             prob_scale = list(tau = 1),
#'                             boundary_sf = rsa_boundary)
#' }
#'
#' @export
predict_establishment = function(
    fit,
    option = c("A","B","C","D","E"),
    calibrate_kappa = FALSE,
    prob_method = c(NULL,"probit","logit","hard"),
    prob_args = list(),
    method = NULL,
    prob_scale = NULL,
    boundary_sf = NULL,
    boundary_params = list(inherit.aes = FALSE, fill = NA, color = "black", size = 0.3)
) {
  # --- Preconditions & argument normalisation ----------------------------------
  `%||%` = function(a, b) if (!is.null(a)) a else b
  stopifnot(inherits(fit, "invasimapr_fit"))

  option = match.arg(option)
  method_in = method %||% prob_method
  prob_method_norm = if (is.null(method_in)) NULL else
    match.arg(as.character(method_in)[1L], c("probit","logit","hard"))

  prob_args_norm = prob_scale %||% prob_args
  if (is.null(prob_args_norm)) prob_args_norm = list()
  if (!is.list(prob_args_norm)) stop("`prob_scale`/`prob_args` must be a list.")

  # --- Select site-/trait-varying sensitivities as required by `option` --------
  GI = switch(
    EXPR = as.character(option),
    "A" = NULL,
    "B" = NULL,
    "C" = fit$sensitivities$theta_i,
    "D" = if (!is.null(fit$sensitivities$site_gamma))
      fit$sensitivities$site_gamma$Gamma_is else NULL,
    "E" = NULL,
    stop("Unknown option: ", option)
  )

  AI = switch(
    EXPR = as.character(option),
               D = if (!is.null(fit$sensitivities$site_alpha)) fit$sensitivities$site_alpha$alpha_is else NULL,
               NULL)

  # --- Compute invasion fitness λ ----------------------------------------------
  fin = compute_invasion_fitness(
    r_is_z   = fit$invaders$r_is_z,
    C_is_z   = fit$invaders$C_is_z,
    S_is_z   = fit$invaders$S_is_z,
    option   = option,
    alpha_i  = fit$sensitivities$alpha_i,
    beta_i   = fit$sensitivities$beta_i,
    theta0   = fit$sensitivities$theta0,
    theta_i  = fit$sensitivities$theta_i,
    Gamma_is = GI,
    alpha_is = AI,
    beta_signed_i = fit$sensitivities$beta_signed_i,
    calibrate_kappa = calibrate_kappa,
    # resident-scale inputs used by calibration (and for completeness)
    r_js_z   = fit$residents$r_js_z,
    C_js_z   = fit$residents$C_js_z,
    S_js_z   = fit$residents$S_js_z,
    Q_res    = fit$traits$Q_res,
    a0 = fit$sensitivities$a0, a1 = fit$sensitivities$a1, a2 = fit$sensitivities$a2,
    b0 = fit$sensitivities$b0, b1 = fit$sensitivities$b1, b2 = fit$sensitivities$b2,
    site_df  = fit$inputs$site_df,
    return_long = TRUE
  )
  fit$fitness = fin

  # --- Optional: map λ → P via probit/logit/hard and overlay boundary ----------
  if (!is.null(prob_method_norm)) {
    fit$prob = do.call(
      what = compute_establishment_probability,
      args = c(list(
        lambda_is    = fin$lambda_is,
        method       = prob_method_norm,
        site_df      = fit$inputs$site_df,  # fill coords for maps / long table
        option_label = fin$option           # carry over the λ option label
      ), prob_args_norm)
    )

    # Add boundary overlay to site-level maps when available
    if (!is.null(boundary_sf) &&
        !is.null(fit$prob$plots) &&
        requireNamespace("ggplot2", quietly = TRUE)) {

      add_boundary = function(p) {
        if (is.null(p) || !"ggplot" %in% class(p)) return(p)
        do.call(`+`, list(
          p,
          do.call(ggplot2::geom_sf, c(list(data = boundary_sf), boundary_params))
        ))
      }
      pl = fit$prob$plots
      if (!is.null(pl$site_mean))  pl$site_mean  = add_boundary(pl$site_mean)
      if (!is.null(pl$expected_n)) pl$expected_n = add_boundary(pl$expected_n)
      if (!is.null(pl$band))       pl$band       = add_boundary(pl$band)
      fit$prob$plots = pl
    }
  } else {
    fit$prob = NULL
  }

  fit
}
