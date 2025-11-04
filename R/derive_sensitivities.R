#' Derive trait-varying sensitivities (alpha_i, beta_i) and abiotic slope (theta, gamma_i)
#'
#' @description
#' From an auxiliary GLMM on standardized predictors, extract fixed-effect
#' systems for `r_z`, `C_z`, and `S_z`, evaluate them at invader trait
#' coordinates, and (optionally) test whether the abiotic slope on `r_z`
#' should vary with traits (`tr1`, `tr2`). Biotic terms are enforced as
#' nonnegative penalties on invasion fitness.
#'
#' @details
#' **Biotic penalties (prior/constraint).** Invasion fitness is modeled as
#' \deqn{\lambda_{is} = \Gamma_{is} \, r^{(z)}_{is} \;-\; \alpha_i \, C^{(z)}_{is} \;-\; \beta_i \, S^{(z)}_{is}.}
#' We therefore set \eqn{\alpha_i = \max(0, -\text{slope}_{C,i})} and
#' \eqn{\beta_i = \max(0, -\text{slope}_{S,i})}, where \eqn{\text{slope}_{C,i}}
#' and \eqn{\text{slope}_{S,i}} are the fitted (possibly trait-varying) slopes
#' for `C_z` and `S_z`. Signed versions (before clamping) are returned for
#' diagnostics.
#'
#' **Testing trait-varying abiotic slope.** The argument `lrt` controls
#' whether and how to test the joint null \eqn{H_0: r_z:tr1 = r_z:tr2 = 0}.
#' - For backward compatibility, `lrt = TRUE` performs a **Wald** joint χ²
#'   test using the fixed-effects variance–covariance (no refit).
#' - `lrt = FALSE` performs no test.
#' - Alternatively, pass a character string: `"wald"`, `"lrt"`, or `"none"`.
#'   - `"lrt"` refits a reduced model (dropping both interactions) and runs a
#'     likelihood-ratio test (`anova(..., test = "Chisq")`).
#'
#' If the chosen test is significant at \eqn{\alpha = 0.05}, the function
#' uses the trait-varying `theta_i`; otherwise a constant `theta0` is used
#' for `gamma_i`.
#'
#' @param fit_coeffs A fitted `glmmTMB` model on `(r_z + C_z + S_z) × (tr1 + tr2)`.
#' @param Q_inv Data frame with invader trait coordinates `tr1`, `tr2`
#'   (rownames = invader IDs).
#' @param inv_ids Character vector of invader IDs specifying the output order.
#' @param lrt Logical **or** character. If logical:
#'   - `TRUE` → Wald test; `FALSE` → no test.
#'   If character, one of `"wald"`, `"lrt"`, or `"none"`. Default is equivalent
#'   to `"wald"` for `TRUE` and `"none"` for `FALSE`.
#'
#' @return A list with components:
#' \describe{
#'   \item{alpha_i, beta_i}{Nonnegative penalties for `C_z` and `S_z`.}
#'   \item{alpha_signed_i, beta_signed_i}{Signed (pre-clamp) slopes.}
#'   \item{theta0}{Intercept for `r_z` slope system.}
#'   \item{theta_i}{Trait-varying `r_z` slope evaluated at invader traits.}
#'   \item{gamma_i}{Either `theta_i` (if test is significant) or a constant
#'                 vector equal to `theta0` (if not).}
#'   \item{a0,a1,a2; b0,b1,b2}{Coefficient components for `C_z` and `S_z`
#'                 systems: main effect and interactions with `tr1`,`tr2`.}
#'   \item{df}{Per-invader tidy data frame with traits, signed slopes,
#'             clamping flags, and `theta_i`/`gamma_i`.}
#'   \item{clamp_summary}{Counts and proportions of clamped \eqn{\alpha} and \eqn{\beta}.}
#'   \item{use_theta}{Logical; `TRUE` if the chosen test was significant.}
#'   \item{lrt}{\strong{Back-compat alias} of `test_tab` (data frame with test result).}
#'   \item{test_tab}{Data frame with test statistic, df, and p-value (or `NULL`).}
#'   \item{test_type}{One of `"wald"`, `"lrt"`, or `"none"`.}
#'   \item{prior_note}{Reminder of the nonnegative-penalty prior.}
#' }
#'
#' @seealso [glmmTMB::glmmTMB()], [stats::anova()]
#' @export
derive_sensitivities = function(fit_coeffs, Q_inv, inv_ids, lrt = TRUE) {
  # ---- Robust input checks ----------------------------------------------------
  stopifnot(is.data.frame(Q_inv), all(c("tr1","tr2") %in% colnames(Q_inv)))
  stopifnot(is.character(inv_ids), length(inv_ids) >= 1L)
  if (length(setdiff(inv_ids, rownames(Q_inv)))) {
    miss_inv <- setdiff(inv_ids, rownames(Q_inv))
    stop("Q_inv is missing invader rows: ", paste(miss_inv, collapse = ", "))
  }

  # ---- Backward-compatible test selection ------------------------------------
  # lrt can be logical (legacy) or character {"wald","lrt","none"}.
  if (is.logical(lrt)) {
    test_type <- if (isTRUE(lrt)) "wald" else "none"
  } else if (is.character(lrt) && length(lrt) == 1L) {
    test_type <- match.arg(tolower(lrt), c("wald","lrt","none"))
  } else {
    stop("`lrt` must be logical or one of {'wald','lrt','none'}.")
  }

  # ---- Extract fixed effects and (if available) vcov --------------------------
  cf <- glmmTMB::fixef(fit_coeffs)$cond                   # named vector of fixed effects
  Vb <- try(stats::vcov(fit_coeffs)$cond, silent = TRUE)  # fixed-effect VCOV (may fail)

  # Helper: extract a main/interaction coefficient safely
  .get_cf <- function(cf, main, by = NULL, default = 0) {
    if (is.null(by)) {
      if (main %in% names(cf)) return(unname(cf[main]))
      return(default)
    }
    nm <- c(paste0(main, ":", by), paste0(by, ":", main))
    nm <- nm[nm %in% names(cf)]
    if (length(nm)) unname(cf[nm[1]]) else default
  }

  # Align Q_inv to requested invader order and keep only the trait columns
  Q_inv <- Q_inv[inv_ids, c("tr1","tr2"), drop = FALSE]

  # ---- Trait-varying slopes for biotic penalties (C_z and S_z) ---------------
  # Slopes are linear functions of traits: a0 + a1*tr1 + a2*tr2, etc.
  a0 <- .get_cf(cf, "C_z")
  a1 <- .get_cf(cf, "C_z", "tr1")
  a2 <- .get_cf(cf, "C_z", "tr2")

  b0 <- .get_cf(cf, "S_z")
  b1 <- .get_cf(cf, "S_z", "tr1")
  b2 <- .get_cf(cf, "S_z", "tr2")

  # Evaluate signed slopes at invader trait coordinates
  slope_C_i <- with(Q_inv, a0 + a1*tr1 + a2*tr2)
  slope_S_i <- with(Q_inv, b0 + b1*tr1 + b2*tr2)

  # Signed versions (before clamping to nonnegative penalties)
  # Note: alpha penalizes C_z with a negative learned slope (so alpha_signed = -slope_C)
  #       beta  penalizes S_z with a positive learned slope (so beta_signed  = +slope_S)
  alpha_signed_i <- stats::setNames(-as.numeric(slope_C_i), rownames(Q_inv))
  beta_signed_i  <- stats::setNames( as.numeric(slope_S_i), rownames(Q_inv))

  # Nonnegative penalties used in lambda (clamp at 0 if sign disagrees with prior)
  alpha_i <- stats::setNames(pmax(0, alpha_signed_i), rownames(Q_inv))
  beta_i  <- stats::setNames(pmax(0, -beta_signed_i),  rownames(Q_inv))

  # Flags indicating which invaders were clamped to 0
  clamped_alpha <- as.logical(slope_C_i >= 0)  # alpha_i set to 0 where slope_C_i >= 0
  clamped_beta  <- as.logical(slope_S_i >= 0)  # beta_i  set to 0 where slope_S_i >= 0

  # ---- Abiotic slope system (theta): main and interactions --------------------
  t0 <- .get_cf(cf, "r_z")
  t1 <- .get_cf(cf, "r_z", "tr1")
  t2 <- .get_cf(cf, "r_z", "tr2")

  theta0  <- if (is.finite(t0)) t0 else 1                     # conservative fallback
  theta_i <- stats::setNames(as.numeric(with(Q_inv, t0 + t1*tr1 + t2*tr2)),
                             rownames(Q_inv))

  # ---- Joint test for trait-varying r_z slope --------------------------------
  # We test H0: r_z:tr1 = r_z:tr2 = 0. If significant, use theta_i; else use theta0.
  test_tab   <- NULL
  use_theta  <- FALSE

  # Names for interaction coefficients (both orders) if present
  nm1 <- intersect(c("r_z:tr1","tr1:r_z"), names(cf))
  nm2 <- intersect(c("r_z:tr2","tr2:r_z"), names(cf))

  if (identical(test_type, "wald")) {
    if (length(nm1) && length(nm2) && !inherits(Vb, "try-error")) {
      b    <- c(cf[nm1[1]], cf[nm2[1]])
      Vsub <- Vb[c(nm1[1], nm2[1]), c(nm1[1], nm2[1]), drop = FALSE]
      ok   <- all(is.finite(Vsub)) && nrow(Vsub) == 2 && ncol(Vsub) == 2
      if (ok) {
        detV <- try(det(Vsub), silent = TRUE)
        if (!inherits(detV,"try-error") && is.finite(detV) && abs(detV) > .Machine$double.eps) {
          W <- as.numeric(t(b) %*% solve(Vsub) %*% b)
          p <- stats::pchisq(W, df = 2L, lower.tail = FALSE)
          test_tab <- data.frame(
            Test = "Wald (r_z:tr1, r_z:tr2)",
            Chisq = W, Df = 2L, `Pr(>Chisq)` = p,
            row.names = NULL, check.names = FALSE
          )
          use_theta <- is.finite(p) && p < 0.05
        }
      }
    }
  } else if (identical(test_type, "lrt")) {
    # Build a reduced model by dropping BOTH interaction terms.
    # Use update() to subtract both colon orders; if absent, update() silently ignores.
    form_full <- stats::formula(fit_coeffs)
    form_red  <- try(stats::update(form_full, . ~ . - r_z:tr1 - r_z:tr2 - tr1:r_z - tr2:r_z),
                     silent = TRUE)
    if (!inherits(form_red, "try-error")) {
      fit_red <- try(stats::update(fit_coeffs, formula = form_red), silent = TRUE)
      if (!inherits(fit_red, "try-error")) {
        a <- try(as.data.frame(stats::anova(fit_coeffs, fit_red, test = "Chisq")),
                 silent = TRUE)
        if (!inherits(a, "try-error") && nrow(a) >= 2L &&
            all(c("Chisq","Df","Pr(>Chisq)") %in% colnames(a))) {
          # The second row is the comparison full vs reduced
          test_tab <- data.frame(
            Test = "LRT (drop r_z:tr1, r_z:tr2)",
            Chisq = a$Chisq[2], Df = a$Df[2], `Pr(>Chisq)` = a$`Pr(>Chisq)`[2],
            row.names = NULL, check.names = FALSE
          )
          p <- as.numeric(a$`Pr(>Chisq)`[2])
          use_theta <- is.finite(p) && p < 0.05
        }
      }
    }
  } else {
    # test_type == "none": do nothing; fall back to constant theta0
  }

  # ---- Map test decision to gamma_i ------------------------------------------
  gamma_i <- if (isTRUE(use_theta)) stats::setNames(as.numeric(theta_i[inv_ids]), inv_ids)
  else                    stats::setNames(rep(theta0, length(inv_ids)), inv_ids)

  # ---- Per-invader tidy table -------------------------------------------------
  df <- data.frame(
    invader         = inv_ids,
    tr1             = Q_inv$tr1,
    tr2             = Q_inv$tr2,
    slope_C_i       = as.numeric(slope_C_i),
    alpha_signed_i  = as.numeric(alpha_signed_i[inv_ids]),
    alpha_i         = as.numeric(alpha_i[inv_ids]),
    slope_S_i       = as.numeric(slope_S_i),
    beta_signed_i   = as.numeric(beta_signed_i[inv_ids]),
    beta_i          = as.numeric(beta_i[inv_ids]),
    theta_i         = as.numeric(theta_i[inv_ids]),
    gamma_i         = as.numeric(gamma_i[inv_ids]),
    clamped_alpha   = as.logical(clamped_alpha),
    clamped_beta    = as.logical(clamped_beta),
    row.names       = NULL,
    check.names     = FALSE
  )

  # ---- Summary of clamping ----------------------------------------------------
  clamp_summary <- list(
    n_invaders          = length(inv_ids),
    n_alpha_clamped     = sum(clamped_alpha, na.rm = TRUE),
    n_beta_clamped      = sum(clamped_beta,  na.rm = TRUE),
    prop_alpha_clamped  = mean(clamped_alpha, na.rm = TRUE),
    prop_beta_clamped   = mean(clamped_beta,  na.rm = TRUE)
  )

  # ---- Return (with back-compat fields) --------------------------------------
  list(
    # Coefficient components (useful for site-varying post-processing)
    a0 = a0, a1 = a1, a2 = a2,
    b0 = b0, b1 = b1, b2 = b2,

    # Penalties (nonnegative) and their signed counterparts
    alpha_i        = alpha_i[inv_ids],
    alpha_signed_i = alpha_signed_i[inv_ids],
    beta_i         = beta_i[inv_ids],
    beta_signed_i  = beta_signed_i[inv_ids],

    # Abiotic slope systems
    theta0         = theta0,
    theta_i        = theta_i[inv_ids],
    gamma_i        = gamma_i,

    # Hypothesis test outputs
    lrt            = test_tab,      # back-compat alias
    test_tab       = test_tab,      # explicit name
    test_type      = test_type,     # "wald", "lrt", or "none"
    use_theta      = use_theta,

    # Diagnostics and notes
    df             = df,
    clamp_summary  = clamp_summary,
    prior_note     = "Biotic effects are modelled as nonnegative penalties; positive learned slopes are reported but not used to increase lambda."
  )
}
