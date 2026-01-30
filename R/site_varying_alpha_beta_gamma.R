#' Site-varying alpha and gamma with trait-dependent beta
#'
#' @description
#' Constructs site-by-invader crowding penalties and density-dependence scalars
#' by combining trait-derived slopes with site-level random effects.
#'
#' Site-level random slopes for \eqn{C_z} and \eqn{r_z} are added to the
#' trait-dependent systems to produce site-varying crowding penalties
#' \eqn{\alpha_{is}} and density-dependence scalars \eqn{\Gamma_{is}}.
#' Saturation effects \eqn{\beta_i} are trait-varying only; no site-varying
#' \eqn{\beta_{is}} is constructed because \eqn{S_z} is site-only.
#'
#' The function also exposes signed slopes and clamping diagnostics returned by
#' \code{derive_sensitivities()}.
#'
#' @param fit_coeffs Fitted residents-only GLMM (e.g. `glmmTMB`) used to extract
#'   fixed effects and site-level random slopes.
#' @param Q_inv Data frame of invader trait scores with columns `tr1` and `tr2`;
#'   row names must correspond to `inv_ids`.
#' @param sites Character vector of site identifiers (rows of the output matrices).
#' @param inv_ids Character vector of invader identifiers (columns of the output matrices).
#' @param lrt Logical; passed to \code{derive_sensitivities()} to control likelihood-ratio
#'   testing for trait effects.
#' @param quiet Logical; if `TRUE`, suppress warnings about missing or negligible
#'   site-level random effects.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{alpha_is}: matrix of site-by-invader crowding penalties
#'   \item \code{Gamma_is}: matrix of site-by-invader density-dependence scalars
#'   \item \code{beta_i}, \code{beta_signed_i}: invader-level saturation effects
#'   \item \code{alpha_signed_i}: signed crowding slopes prior to clamping
#'   \item \code{theta_i}, \code{gamma_i}: invader-level abiotic sensitivities
#'   \item \code{slope_C_i}: trait-derived crowding slopes
#'   \item \code{delta_C_s}, \code{delta_r_s}: site-level random effects
#'   \item \code{notes}: character vector of diagnostics
#'   \item \code{df}: tidy long-format data frame
#' }
#'
#' @export
site_varying_alpha_beta_gamma = function(fit_coeffs, Q_inv, sites, inv_ids,
                                          lrt = TRUE, quiet = FALSE) {

  stopifnot(is.data.frame(Q_inv), all(c("tr1","tr2") %in% names(Q_inv)))

  # 1) Trait-side systems
  sens = derive_sensitivities(fit_coeffs, Q_inv, inv_ids, lrt = lrt)

  # Recompute slope_C_i (signed) for combining with site deltas
  Qi = Q_inv[inv_ids, c("tr1","tr2"), drop = FALSE]
  slope_C_i = with(Qi, sens$a0 + sens$a1*tr1 + sens$a2*tr2)
  slope_C_i = stats::setNames(as.numeric(slope_C_i), inv_ids)

  # 2) Site random slopes
  re_site = try(ranef(fit_coeffs)$cond$site, silent = TRUE)
  notes = character()
  if (inherits(re_site, "try-error")) {
    delta_C_s = stats::setNames(rep(0, length(sites)), sites)
    delta_r_s = delta_C_s
    notes = c(notes, "No site-level random effects; using delta_C_s = 0 and delta_r_s = 0.")
  } else {
    cols = colnames(re_site)
    jC = grep("(^|:)[Cc](?:_[A-Za-z]+)?_z($|:)", cols, perl = TRUE)
    jR = grep("(^|:)[Rr](?:_[A-Za-z]+)?_z($|:)", cols, perl = TRUE)

    if (!length(jC)) {
      delta_C_s = stats::setNames(rep(0, nrow(re_site)), rownames(re_site))
      notes = c(notes, paste0("Random slope for C not found. RE cols: ",
                               paste(cols, collapse = ", "), ". Using delta_C_s = 0."))
    } else {
      delta_C_s = stats::setNames(as.numeric(re_site[, jC[1], drop = TRUE]), rownames(re_site))
      if (stats::sd(delta_C_s, na.rm = TRUE) < 1e-8)
        notes = c(notes, "delta_C_s ~ 0 variance; alpha_is similar alpha_i across sites.")
    }

    if (!length(jR)) {
      delta_r_s = stats::setNames(rep(0, nrow(re_site)), rownames(re_site))
      notes = c(notes, paste0("Random slope for r not found. RE cols: ",
                               paste(cols, collapse = ", "), ". Using delta_r_s = 0."))
    } else {
      delta_r_s = stats::setNames(as.numeric(re_site[, jR[1], drop = TRUE]), rownames(re_site))
      if (stats::sd(delta_r_s, na.rm = TRUE) < 1e-8)
        notes = c(notes, "delta_r_s ~ 0 variance; Gamma_is similar theta_i across sites.")
    }
  }

  # Align to requested sites
  delta_C_s = delta_C_s[sites]; delta_C_s[is.na(delta_C_s)] = 0
  delta_r_s = delta_r_s[sites]; delta_r_s[is.na(delta_r_s)] = 0

  # 3) Combine to site-varying fields
  SLOPE_C_is = outer(delta_C_s, as.numeric(slope_C_i), `+`)
  dimnames(SLOPE_C_is) = list(sites, inv_ids)

  alpha_is = -SLOPE_C_is; alpha_is[alpha_is < 0 | !is.finite(alpha_is)] = 0
  Gamma_is = outer(delta_r_s, as.numeric(sens$theta_i), `+`)
  dimnames(Gamma_is) = list(sites, inv_ids)

  # 4) Tidy long table (now includes alpha_signed_i for reference)
  df = data.frame(
    site            = rep(sites,  times = length(inv_ids)),
    invader         = rep(inv_ids, each  = length(sites)),
    alpha_is        = as.vector(alpha_is),
    SLOPE_C_is      = as.vector(SLOPE_C_is),
    slope_C_i       = rep(as.numeric(slope_C_i), each = length(sites)),
    alpha_signed_i  = rep(as.numeric(sens$alpha_signed_i), each = length(sites)),
    delta_C_s       = rep(as.numeric(delta_C_s), times = length(inv_ids)),
    Gamma_is        = as.vector(Gamma_is),
    theta_i         = rep(as.numeric(sens$theta_i), each = length(sites)),
    delta_r_s       = rep(as.numeric(delta_r_s), times = length(inv_ids)),
    beta_i          = rep(as.numeric(sens$beta_i), each = length(sites)),
    beta_signed_i   = rep(as.numeric(sens$beta_signed_i), each = length(sites)),
    row.names       = NULL, check.names = FALSE
  )

  if (!quiet && length(notes)) warning(paste(notes, collapse = " "))

  list(
    alpha_is       = alpha_is,
    Gamma_is       = Gamma_is,
    beta_i         = sens$beta_i,
    beta_signed_i  = sens$beta_signed_i,
    alpha_signed_i = sens$alpha_signed_i,
    theta_i        = sens$theta_i,
    gamma_i        = sens$gamma_i,   # for options where you want gamma_i (not Gamma_{is})
    slope_C_i      = slope_C_i,
    delta_C_s      = delta_C_s,
    delta_r_s      = delta_r_s,
    notes          = notes,
    df             = df
  )
}
