# ======================================================================
# 5) INVASION FITNESS λ and (optionally) establishment probability
#    DROP-IN: add boundary layer support for returned plots
# ======================================================================

#' @title Compute invasion fitness (λ) and optional establishment probability (P)
#' @param fit from \code{predict_invaders()}.
#' @param option character in \code{c("A","B","C","D","E")}.
#' @param calibrate_kappa logical; set κ so mean resident λ ≈ 0.
#' @param prob_method (legacy) NULL or one of \code{"probit","logit","hard"}.
#' @param prob_args (legacy) list of args for \code{compute_establishment_probability()}.
#' @param method (alias) same as \code{prob_method}; preferred in user code.
#' @param prob_scale (alias) same as \code{prob_args}; preferred in user code.
#' @param boundary_sf An \strong{sf} object to overlay as a boundary on map-like plots.
#' @param boundary_params Named list of aesthetics for \code{ggplot2::geom_sf};
#'   defaults to \code{list(inherit.aes=FALSE, fill=NA, color="black", size=0.3)}.
#' @return updated \code{invasimapr_fit}
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
  `%||%` = function(a, b) if (!is.null(a)) a else b
  # try({ source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/compute_invasion_fitness.R") }, silent = TRUE)
  # if (!exists("compute_invasion_fitness")) stop("compute_invasion_fitness() not found.")
  # try({ source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/compute_establishment_probability.R") }, silent = TRUE)
  # if (!exists("compute_establishment_probability")) stop("compute_establishment_probability() not found.")
  stopifnot(inherits(fit, "invasimapr_fit"))

  option = match.arg(option)
  method_in = method %||% prob_method
  prob_method_norm = if (is.null(method_in)) NULL else match.arg(as.character(method_in)[1L], c("probit","logit","hard"))
  prob_args_norm = prob_scale %||% prob_args; if (is.null(prob_args_norm)) prob_args_norm = list()
  if (!is.list(prob_args_norm)) stop("`prob_scale`/`prob_args` must be a list.")

  GI = switch(option,
               A = NULL,
               B = NULL,
               C = fit$sensitivities$theta_i,
               D = if (!is.null(fit$sensitivities$site_gamma)) fit$sensitivities$site_gamma$Gamma_is else NULL,
               E = NULL)
  AI = switch(option,
               D = if (!is.null(fit$sensitivities$site_alpha)) fit$sensitivities$site_alpha$alpha_is else NULL,
               NULL)

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

  if (!is.null(prob_method_norm)) {
    fit$prob = do.call(
      what = compute_establishment_probability,
      args = c(list(
        lambda_is   = fin$lambda_is,
        method      = prob_method_norm,
        # give site coords so plots/long-table are populated
        site_df     = fit$inputs$site_df,
        # label plots/long table with the option used to make λ
        option_label= fin$option
      ),
      prob_args_norm)
    )

    if (!is.null(boundary_sf) &&
        !is.null(fit$prob$plots) &&
        requireNamespace("ggplot2", quietly = TRUE)) {
      add_boundary = function(p) {
        if (is.null(p) || !"ggplot" %in% class(p)) return(p)
        do.call(`+`, list(p, do.call(ggplot2::geom_sf, c(list(data = boundary_sf), boundary_params))))
      }
      pl = fit$prob$plots
      if (!is.null(pl$site_mean))  pl$site_mean  = add_boundary(pl$site_mean)
      fit$prob$plots = pl
    }
  } else {
    fit$prob = NULL
  }

  fit
}

