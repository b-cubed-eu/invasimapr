#' Site-varying α_{is} and Γ_{is}, plus β_i (and signed variants)
#'
#' @description
#' Adds site random-slope adjustments for C_z and r_z to produce α_{is} and Γ_{is}.
#' β is trait-varying only (no site-varying β_{is} because S_z is site-only).
#' Also surfaces signed slopes and clamping info from `derive_sensitivities()`.
#' @return list(alpha_is, Gamma_is, beta_i, beta_signed_i, alpha_signed_i,
#'              theta_i, gamma_i, slope_C_i, delta_C_s, delta_r_s, notes, df)
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
        notes = c(notes, "delta_C_s ~ 0 variance; α_is ≈ α_i across sites.")
    }

    if (!length(jR)) {
      delta_r_s = stats::setNames(rep(0, nrow(re_site)), rownames(re_site))
      notes = c(notes, paste0("Random slope for r not found. RE cols: ",
                               paste(cols, collapse = ", "), ". Using delta_r_s = 0."))
    } else {
      delta_r_s = stats::setNames(as.numeric(re_site[, jR[1], drop = TRUE]), rownames(re_site))
      if (stats::sd(delta_r_s, na.rm = TRUE) < 1e-8)
        notes = c(notes, "delta_r_s ~ 0 variance; Γ_is ≈ θ_i across sites.")
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
    gamma_i        = sens$gamma_i,   # for options where you want γ_i (not Γ_{is})
    slope_C_i      = slope_C_i,
    delta_C_s      = delta_C_s,
    delta_r_s      = delta_r_s,
    notes          = notes,
    df             = df
  )
}
