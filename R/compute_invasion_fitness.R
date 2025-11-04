#' Compute invasion fitness (five options) with optional resident calibration
#'
#' @title Compute invasion fitness \eqn{\lambda_{is}} for multiple model options
#'
#' @description
#' `compute_invasion_fitness()` evaluates the invasion-fitness surface
#' \eqn{\lambda_{is}} over sites (\eqn{s}) and invaders (\eqn{i}) using the
#' three standardized predictors:
#' * **Abiotic suitability** \eqn{r^{(z)}_{is}} (alignment between invader traits and local environment),
#' * **Niche crowding** \eqn{C^{(z)}_{is}} (overlap with resident trait space weighted by composition),
#' * **Resident competition** \eqn{S^{(z)}_{is}} (site-only saturation).
#'
#' It supports five published/useful variants:
#' **A**: \eqn{\gamma=1}, **B**: global \eqn{\theta_0}, **C**: trait-varying \eqn{\theta_i},
#' **D**: site-varying \eqn{\Gamma_{is}} and \eqn{\alpha_{is}}, **E**: signed \eqn{S} effect.
#' Optionally, it calibrates an offset \eqn{\kappa} so the **mean resident** \eqn{\lambda}
#' is approximately zero on the same scale (when resident standardized matrices are supplied).
#'
#' @param r_is_z Matrix \eqn{S \times I} of standardized abiotic suitability for invaders.
#'               Row names are site IDs; column names are invader IDs.
#' @param C_is_z Matrix \eqn{S \times I} of standardized niche crowding for invaders.
#' @param S_is_z Matrix \eqn{S \times I} of standardized site saturation for invaders
#'               (site-only; same value down each site row).
#' @param option Character, one of `c("A","B","C","D","E")`. See Details.
#' @param alpha_i Optional named numeric vector (length \eqn{I}) of invader crowding
#'                sensitivities \eqn{\alpha_i}. Required for A–C and E; ignored if `alpha_is` is provided in D.
#' @param beta_i  Optional named numeric vector (length \eqn{I}) of saturation sensitivities
#'                \eqn{\beta_i}. Used in A–D. For E, use `beta_signed_i` for a signed effect.
#' @param theta0  Numeric scalar global slope for \eqn{r^{(z)}} (Option B; default `1`).
#' @param theta_i Optional named numeric vector (length \eqn{I}) of trait-varying slopes
#'                for \eqn{r^{(z)}} (Option C).
#' @param Gamma_is Optional matrix \eqn{S \times I} of site-varying slopes for \eqn{r^{(z)}}
#'                 (Option D). If NULL, falls back to `theta_i` or `theta0` as available.
#' @param alpha_is Optional matrix \eqn{S \times I} of site-varying penalties for crowding
#'                 (Option D). If NULL, falls back to broadcasting `alpha_i`.
#' @param beta_signed_i Optional named numeric vector (length \eqn{I}) providing a signed
#'                 \eqn{\beta} for Option E (facilitation allowed when \eqn{\beta<0}).
#' @param calibrate_kappa Logical; if `TRUE`, compute \eqn{\kappa} from **residents** so that
#'                        the mean resident \eqn{\lambda} is ~0 on the same scale (see below).
#' @param r_js_z,C_js_z,S_js_z Optional resident standardized matrices (site \eqn{\times} resident),
#'                        required only when `calibrate_kappa = TRUE`.
#' @param Q_res Optional data frame of resident trait-plane scores with columns `tr1`, `tr2`
#'              and rownames = resident IDs (only needed for \eqn{\kappa} if you want to
#'              derive resident analogs \eqn{\alpha_j}, \eqn{\beta_j} from trait slopes).
#' @param a0,a1,a2,b0,b1,b2 Optional numeric scalars: coefficients for
#'              resident analog slopes \eqn{\text{slope}_C = a_0 + a_1\,tr1 + a_2\,tr2}
#'              and \eqn{\text{slope}_S = b_0 + b_1\,tr1 + b_2\,tr2}. Used only when
#'              `calibrate_kappa=TRUE` and `Q_res` provided.
#' @param site_df Optional site table with columns `site`, `x`, `y` used to enrich the
#'                tidy output when `return_long=TRUE`.
#' @param return_long Logical; if `TRUE`, returns a tidy table with predictors and coefficients.
#' @param label Optional character to tag the option in the tidy output; defaults to a
#'              descriptive name per option (e.g., `"Option A (γ=1)"`).
#'
#' @details
#' **Options**
#' *A*: \eqn{\lambda_{is} = 1\cdot r^{(z)}_{is} - \alpha_i C^{(z)}_{is} - \beta_i S^{(z)}_{is} + \kappa}
#' *B*: \eqn{\lambda_{is} = \theta_0 r^{(z)}_{is} - \alpha_i C^{(z)}_{is} - \beta_i S^{(z)}_{is} + \kappa}
#' *C*: \eqn{\lambda_{is} = \theta_i r^{(z)}_{is} - \alpha_i C^{(z)}_{is} - \beta_i S^{(z)}_{is} + \kappa}
#' *D*: \eqn{\lambda_{is} = \Gamma_{is} r^{(z)}_{is} - \alpha_{is} C^{(z)}_{is} - \beta_i S^{(z)}_{is} + \kappa}
#' *E*: \eqn{\lambda_{is} = \theta_0 r^{(z)}_{is} - \alpha_i C^{(z)}_{is} + \beta^{(\mathrm{signed})}_i S^{(z)}_{is} + \kappa}
#'
#' **Resident calibration** (\eqn{\kappa}). When `calibrate_kappa=TRUE`, the function
#' computes resident analog penalties \eqn{\alpha_j,\beta_j} from `Q_res` and slopes
#' `{a0,a1,a2,b0,b1,b2}`, builds resident matrices `AJ`, `BJ`, evaluates
#' \eqn{\lambda^{(\mathrm{res})}_{js} = \theta_0 r^{(z)}_{js} - \alpha_j C^{(z)}_{js} - \beta_j S^{(z)}_{js}},
#' then sets \eqn{\kappa = -\mathrm{mean}(\lambda^{(\mathrm{res})}_{js})} (na-rm),
#' keeping residents and invaders on the same scale.
#'
#' @return A list with:
#' * `lambda_is` — matrix \eqn{S \times I} of invasion fitness,
#' * `GI` — the \eqn{\gamma} used (vector length \eqn{I} or matrix \eqn{S \times I}),
#' * `AI` — the \eqn{\alpha} used (vector length \eqn{I} or matrix \eqn{S \times I}),
#' * `BI` — the \eqn{\beta} used (vector length \eqn{I}),
#' * `kappa` — numeric calibration offset,
#' * `option` — option label,
#' * `lambda_long` — tidy tibble (if `return_long=TRUE`).
#'
#' @examples
#' # Minimal reproducible example (toy shapes only)
#' S = 5; I = 3
#' set.seed(1)
#' r_is_z = matrix(rnorm(S*I),  S, I, dimnames=list(paste0("s",1:S), paste0("i",1:I)))
#' C_is_z = matrix(rnorm(S*I),  S, I, dimnames=dimnames(r_is_z))
#' S_is_z = matrix(rep(scale(rnorm(S)), each=I), S, I, dimnames=dimnames(r_is_z)) # site-only broadcast
#' alpha_i = setNames(runif(ncol(r_is_z), 0.2, 1.0), colnames(r_is_z))
#' beta_i  = setNames(runif(ncol(r_is_z), 0.1, 0.5), colnames(r_is_z))
#'
#' # Option A (γ=1, κ=0)
#' outA = compute_invasion_fitness(r_is_z, C_is_z, S_is_z,
#'                                  option="A", alpha_i=alpha_i, beta_i=beta_i,
#'                                  theta0=1, return_long=FALSE)
#'
#' # Option B (γ = θ0)
#' outB = compute_invasion_fitness(r_is_z, C_is_z, S_is_z,
#'                                  option="B", alpha_i=alpha_i, beta_i=beta_i,
#'                                  theta0=0.8, return_long=FALSE)
#'
#' # Option C (γ_i = θ_i) — use I = ncol(r_is_z) to avoid length/name mismatch
#' I = ncol(r_is_z)
#' theta_i = setNames(runif(I, 0.5, 1.2), colnames(r_is_z))
#' outC = compute_invasion_fitness(r_is_z, C_is_z, S_is_z,
#'                                  option="C", alpha_i=alpha_i, beta_i=beta_i,
#'                                  theta_i=theta_i, return_long=FALSE)
#'
#' # Option D (site-varying Γ_is and α_is)
#' Gamma_is = matrix(rep(theta_i, each=nrow(r_is_z)), nrow=nrow(r_is_z), dimnames=dimnames(r_is_z))
#' alpha_is = pmax(0, matrix(rep(alpha_i, each=nrow(r_is_z)), nrow=nrow(r_is_z),
#'                            dimnames=dimnames(r_is_z)) + matrix(rnorm(prod(dim(r_is_z)), 0, 0.1),
#'                                                                nrow=nrow(r_is_z)))
#' outD = compute_invasion_fitness(r_is_z, C_is_z, S_is_z,
#'                                  option="D", alpha_is=alpha_is, beta_i=beta_i,
#'                                  Gamma_is=Gamma_is, return_long=FALSE)
#'
#' # Option E (signed S effect)
#' beta_signed_i = setNames(rnorm(ncol(r_is_z), 0, 0.3), colnames(r_is_z))
#' outE = compute_invasion_fitness(r_is_z, C_is_z, S_is_z,
#'                                  option="E", alpha_i=alpha_i, beta_signed_i=beta_signed_i,
#'                                  theta0=1, return_long=FALSE)
#'
#' @export
compute_invasion_fitness = function(
    r_is_z, C_is_z, S_is_z,
    option = c("A","B","C","D","E"),
    alpha_i = NULL, beta_i = NULL,
    theta0 = 1, theta_i = NULL,
    Gamma_is = NULL, alpha_is = NULL,
    beta_signed_i = NULL,
    calibrate_kappa = FALSE,
    r_js_z = NULL, C_js_z = NULL, S_js_z = NULL,
    Q_res = NULL, a0 = NULL, a1 = NULL, a2 = NULL, b0 = NULL, b1 = NULL, b2 = NULL,
    site_df = NULL, return_long = TRUE, label = NULL
){
  option = match.arg(option)

  # --- Coerce to matrix (robust to data.frame inputs) --------------------------
  r_is_z = if (is.matrix(r_is_z)) r_is_z else as.matrix(r_is_z)
  C_is_z = if (is.matrix(C_is_z)) C_is_z else as.matrix(C_is_z)
  S_is_z = if (is.matrix(S_is_z)) S_is_z else as.matrix(S_is_z)

  stopifnot(identical(dim(r_is_z), dim(C_is_z)), identical(dim(r_is_z), dim(S_is_z)))
  stopifnot(!is.null(rownames(r_is_z)), !is.null(colnames(r_is_z)))
  sites   = rownames(r_is_z)
  inv_ids = colnames(r_is_z)
  S = length(sites); I = length(inv_ids)

  # Flexible aligner for invader-vectors (γ_i, α_i, β_i, etc.)
  .align_inv_vec = function(x, nm_expected, allow_scalar = TRUE, arg = deparse(substitute(x))) {
    if (is.null(x)) return(NULL)
    if (!is.null(names(x))) {
      out = x[nm_expected]
      if (any(is.na(out))) {
        miss = nm_expected[is.na(out)]
        stop(sprintf("%s is missing names for: %s", arg, paste(miss, collapse=", ")))
      }
      return(out)
    } else {
      if (length(x) == length(nm_expected)) {
        names(x) = nm_expected
        return(x)
      }
      if (allow_scalar && length(x) == 1L) {
        return(stats::setNames(rep(x, length(nm_expected)), nm_expected))
      }
      stop(sprintf("%s must be named by invader IDs or have length 1 or %d.", arg, length(nm_expected)))
    }
  }

  # --- Choose coefficients per option -----------------------------------------
  if (option == "A") {
    GI = matrix(1, nrow = S, ncol = I, dimnames = list(sites, inv_ids))
    AI = matrix(rep(.align_inv_vec(alpha_i, inv_ids), each = S), nrow = S, dimnames = list(sites, inv_ids))
    BI = .align_inv_vec(beta_i, inv_ids)
    if (is.null(label)) label = "Option A (γ=1)"

  } else if (option == "B") {
    GI = matrix(theta0, nrow = S, ncol = I, dimnames = list(sites, inv_ids))
    AI = matrix(rep(.align_inv_vec(alpha_i, inv_ids), each = S), nrow = S, dimnames = list(sites, inv_ids))
    BI = .align_inv_vec(beta_i, inv_ids)
    if (is.null(label)) label = "Option B (γ=θ0)"

  } else if (option == "C") {
    GI = matrix(rep(.align_inv_vec(theta_i, inv_ids), each = S), nrow = S, dimnames = list(sites, inv_ids))
    AI = matrix(rep(.align_inv_vec(alpha_i, inv_ids),  each = S), nrow = S, dimnames = list(sites, inv_ids))
    BI = .align_inv_vec(beta_i, inv_ids)
    if (is.null(label)) label = "Option C (γ=θ_i)"

  } else if (option == "D") {
    if (!is.null(Gamma_is)) {
      Gamma_is = if (is.matrix(Gamma_is)) Gamma_is else as.matrix(Gamma_is)
      stopifnot(identical(dim(Gamma_is), dim(r_is_z)))
      GI = Gamma_is
    } else if (!is.null(theta_i)) {
      GI = matrix(rep(.align_inv_vec(theta_i, inv_ids), each = S), nrow = S, dimnames = list(sites, inv_ids))
    } else {
      GI = matrix(theta0, nrow = S, ncol = I, dimnames = list(sites, inv_ids))
    }

    if (!is.null(alpha_is)) {
      alpha_is = if (is.matrix(alpha_is)) alpha_is else as.matrix(alpha_is)
      stopifnot(identical(dim(alpha_is), dim(r_is_z)))
      AI = alpha_is
    } else {
      AI = matrix(rep(.align_inv_vec(alpha_i, inv_ids), each = S), nrow = S, dimnames = list(sites, inv_ids))
    }

    BI = .align_inv_vec(beta_i, inv_ids)
    if (is.null(label)) label = "Option D (Γ_is, α_is)"

  } else if (option == "E") {
    GI = matrix(theta0, nrow = S, ncol = I, dimnames = list(sites, inv_ids))
    AI = matrix(rep(.align_inv_vec(alpha_i, inv_ids), each = S), nrow = S, dimnames = list(sites, inv_ids))
    BI = .align_inv_vec(beta_signed_i, inv_ids)
    if (is.null(label)) label = "Option E (signed S)"
  }

  # --- Optional resident calibration (kappa) -----------------------------------
  kappa = 0
  if (isTRUE(calibrate_kappa)) {
    if (is.null(r_js_z) || is.null(C_js_z) || is.null(S_js_z))
      stop("calibrate_kappa=TRUE requires resident matrices: r_js_z, C_js_z, S_js_z.")
    r_js_z = if (is.matrix(r_js_z)) r_js_z else as.matrix(r_js_z)
    C_js_z = if (is.matrix(C_js_z)) C_js_z else as.matrix(C_js_z)
    S_js_z = if (is.matrix(S_js_z)) S_js_z else as.matrix(S_js_z)

    if (!is.null(Q_res) && all(vapply(list(a0,a1,a2,b0,b1,b2), function(x) !is.null(x), logical(1)))) {
      slope_C_j = with(Q_res, a0 + a1*tr1 + a2*tr2)
      slope_S_j = with(Q_res, b0 + b1*tr1 + b2*tr2)
      alpha_j = pmax(0, -as.numeric(slope_C_j)); names(alpha_j) = rownames(Q_res)
      beta_j  = pmax(0, -as.numeric(slope_S_j)); names(beta_j) = rownames(Q_res)
    } else {
      alpha_j = stats::setNames(rep(1, ncol(C_js_z)), colnames(C_js_z))
      beta_j  = stats::setNames(rep(1, ncol(S_js_z)), colnames(S_js_z))
    }

    AJ = matrix(alpha_j[colnames(C_js_z)], nrow=nrow(C_js_z), ncol=ncol(C_js_z),
                 byrow=TRUE, dimnames=dimnames(C_js_z))
    BJ = matrix(beta_j [colnames(S_js_z)], nrow=nrow(S_js_z), ncol=ncol(S_js_z),
                 byrow=TRUE, dimnames=dimnames(S_js_z))

    lambda_js_res = theta0 * r_js_z - AJ * C_js_z - BJ * S_js_z
    kappa = -mean(lambda_js_res, na.rm = TRUE)
  }

  # --- Assemble λ_is -----------------------------------------------------------
  lambda_is = GI * r_is_z - AI * C_is_z - sweep(S_is_z, 2, .align_inv_vec(BI, inv_ids), `*`) + kappa

  # --- Optional tidy table -----------------------------------------------------
  # --- Optional tidy table -----------------------------------------------------
  lambda_long = NULL
  if (isTRUE(return_long)) {
    # Try a user helper if present, but fall back robustly
    lambda_long_try = NULL
    if (exists("make_lambda_long")) {
      lambda_long_try = try(
        make_lambda_long(
          lambda_is, label, r_is_z, C_is_z, S_is_z,
          site_df = if (!is.null(site_df)) site_df else tibble::tibble(site = rownames(lambda_is)),
          gamma   = if (is.matrix(GI)) NULL else stats::setNames(as.numeric(GI[1, ]), colnames(GI)),
          alpha   = if (is.matrix(AI)) NULL else stats::setNames(as.numeric(AI[1, ]), colnames(AI)),
          beta    = .align_inv_vec(BI, inv_ids)
        ),
        silent = TRUE
      )
      if (inherits(lambda_long_try, "try-error")) lambda_long_try = NULL
    }

    if (is.null(lambda_long_try) || NROW(lambda_long_try) == 0L) {
      # Robust default
      lambda_long = tibble::tibble(
        site    = rep(rownames(lambda_is), times = ncol(lambda_is)),
        invader = rep(colnames(lambda_is), each  = nrow(lambda_is)),
        lambda  = as.vector(lambda_is),
        r_z     = as.vector(r_is_z),
        C_z     = as.vector(C_is_z),
        S_z     = as.vector(S_is_z),
        option  = label
      )
      if (!is.null(site_df) && all(c("site","x","y") %in% names(site_df))) {
        lambda_long = dplyr::left_join(lambda_long, site_df, by = "site")
      }
    } else {
      lambda_long = lambda_long_try
    }
  }


  list(
    lambda_is   = lambda_is,
    GI          = GI,
    AI          = AI,
    BI          = .align_inv_vec(BI, inv_ids),
    kappa       = kappa,
    option      = label,
    lambda_long = lambda_long
  )
}
