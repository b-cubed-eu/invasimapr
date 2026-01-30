#' Compute invasion fitness with multiple model options
#'
#' @title Compute invasion fitness \eqn{\lambda_{is}} for multiple model options
#'
#' @description
#' \code{compute_invasion_fitness()} evaluates the invasion-fitness surface
#' \eqn{\lambda_{is}} over sites (\eqn{s}) and invaders (\eqn{i}) using three
#' standardized predictors:
#'
#' \itemize{
#'   \item Abiotic suitability \eqn{r^{(z)}_{is}} (alignment between invader traits
#'   and the local environment),
#'   \item Niche crowding \eqn{C^{(z)}_{is}} (overlap with resident trait space,
#'   weighted by composition),
#'   \item Resident competition \eqn{S^{(z)}_{is}} (site-level saturation).
#' }
#'
#' Five model variants are supported:
#' \itemize{
#'   \item Option A: \eqn{\gamma = 1}
#'   \item Option B: global slope \eqn{\theta_0}
#'   \item Option C: invader-specific slopes \eqn{\theta_i}
#'   \item Option D: site-varying slopes \eqn{\Gamma_{is}} and penalties
#'     \eqn{\alpha_{is}}
#'   \item Option E: signed saturation effect
#' }
#'
#' Optionally, an offset \eqn{\kappa} can be calibrated so that the mean resident
#' invasion fitness is approximately zero when resident standardized matrices
#' are supplied.
#'
#' @param r_is_z Matrix of standardized abiotic suitability with dimensions
#'   \eqn{S} by \eqn{I}.
#' @param C_is_z Matrix of standardized niche crowding with dimensions
#'   \eqn{S} by \eqn{I}.
#' @param S_is_z Matrix of standardized site saturation with dimensions
#'   \eqn{S} by \eqn{I}.
#' @param option Character string specifying the model option. One of
#'   \code{"A"}, \code{"B"}, \code{"C"}, \code{"D"}, or \code{"E"}.
#' @param alpha_i Optional named vector of invader-level crowding sensitivities.
#' @param beta_i Optional named vector of invader-level saturation sensitivities.
#' @param theta0 Global abiotic slope used in options A, B, and E.
#' @param theta_i Optional invader-specific abiotic slopes used in option C.
#' @param Gamma_is Optional site by invader matrix of abiotic slopes used in
#'   option D.
#' @param alpha_is Optional site by invader matrix of crowding penalties used
#'   in option D.
#' @param beta_signed_i Optional signed saturation sensitivities used in
#'   option E.
#' @param calibrate_kappa Logical; if \code{TRUE}, compute a calibration offset
#'   using resident data.
#' @param r_js_z,C_js_z,S_js_z Optional resident standardized matrices required
#'   when \code{calibrate_kappa = TRUE}.
#' @param Q_res Optional data frame of resident trait-plane scores.
#' @param a0,a1,a2,b0,b1,b2 Optional numeric coefficients used to derive resident
#'   analog slopes when calibrating \eqn{\kappa}.
#' @param site_df Optional site metadata table with columns \code{site}, \code{x},
#'   and \code{y}.
#' @param return_long Logical; if \code{TRUE}, return a tidy long-format table.
#' @param label Optional character label attached to the output.
#'
#' @details
#' The invasion fitness is computed as:
#'
#' \itemize{
#'   \item Option A:
#'     \eqn{\lambda_{is} = r^{(z)}_{is} - \alpha_i C^{(z)}_{is}
#'     - \beta_i S^{(z)}_{is} + \kappa}
#'   \item Option B:
#'     \eqn{\lambda_{is} = \theta_0 r^{(z)}_{is} - \alpha_i C^{(z)}_{is}
#'     - \beta_i S^{(z)}_{is} + \kappa}
#'   \item Option C:
#'     \eqn{\lambda_{is} = \theta_i r^{(z)}_{is} - \alpha_i C^{(z)}_{is}
#'     - \beta_i S^{(z)}_{is} + \kappa}
#'   \item Option D:
#'     \eqn{\lambda_{is} = \Gamma_{is} r^{(z)}_{is} - \alpha_{is} C^{(z)}_{is}
#'     - \beta_i S^{(z)}_{is} + \kappa}
#'   \item Option E:
#'     \eqn{\lambda_{is} = \theta_0 r^{(z)}_{is} - \alpha_i C^{(z)}_{is}
#'     + \beta_i^{(signed)} S^{(z)}_{is} + \kappa}
#' }
#'
#' When \code{calibrate_kappa = TRUE}, the offset is chosen so that the mean
#' resident invasion fitness equals zero.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{lambda_is}: invasion fitness matrix
#'   \item \code{GI}: abiotic slope used
#'   \item \code{AI}: crowding penalties used
#'   \item \code{BI}: saturation penalties used
#'   \item \code{kappa}: calibration offset
#'   \item \code{option}: option label
#'   \item \code{lambda_long}: tidy table (optional)
#' }
#'
#' @examples
#' ## ---------------------------------------------------------------
#' ## Minimal reproducible example (toy dimensions, no real ecology)
#' ## ---------------------------------------------------------------
#'
#' ## S = number of sites, I = number of invaders
#' S = 5
#' I = 3
#' set.seed(1)
#'
#' ## r_is_z : intrinsic growth component
#' ## rows = sites (s), columns = invaders (i)
#' r_is_z = matrix(
#'   rnorm(S * I),
#'   nrow = S,
#'   ncol = I,
#'   dimnames = list(paste0("s", 1:S), paste0("i", 1:I))
#' )
#'
#' ## C_is_z : crowding / competition component (same shape as r_is_z)
#' C_is_z = matrix(
#'   rnorm(S * I),
#'   nrow = S,
#'   ncol = I,
#'   dimnames = dimnames(r_is_z)
#' )
#'
#' ## S_is_z : site-only environmental gradient
#' ## generated per-site, then broadcast across invaders
#' S_is_z = matrix(
#'   rep(scale(rnorm(S)), each = I),
#'   nrow = S,
#'   ncol = I,
#'   dimnames = dimnames(r_is_z)
#' )
#'
#' ## Invader-specific baseline parameters
#' alpha_i = setNames(runif(I, 0.2, 1.0), colnames(r_is_z))  # baseline crowding sensitivity
#' beta_i  = setNames(runif(I, 0.1, 0.5), colnames(r_is_z))  # strength of S effect
#'
#' ## ---------------------------------------------------------------
#' ## Option A: fixed gamma = 1, no S effect (kappa = 0)
#' ## ---------------------------------------------------------------
#' outA = compute_invasion_fitness(
#'   r_is_z, C_is_z, S_is_z,
#'   option      = "A",
#'   alpha_i     = alpha_i,
#'   beta_i      = beta_i,
#'   theta0      = 1,
#'   return_long = FALSE
#' )
#'
#' ## ---------------------------------------------------------------
#' ## Option B: gamma shared across invaders (gamma = theta_0)
#' ## ---------------------------------------------------------------
#' outB = compute_invasion_fitness(
#'   r_is_z, C_is_z, S_is_z,
#'   option      = "B",
#'   alpha_i     = alpha_i,
#'   beta_i      = beta_i,
#'   theta0      = 0.8,
#'   return_long = FALSE
#' )
#'
#' ## ---------------------------------------------------------------
#' ## Option C: invader-specific gamma_i
#' ## ---------------------------------------------------------------
#' theta_i = setNames(
#'   runif(I, 0.5, 1.2),
#'   colnames(r_is_z)
#' )
#'
#' outC = compute_invasion_fitness(
#'   r_is_z, C_is_z, S_is_z,
#'   option      = "C",
#'   alpha_i     = alpha_i,
#'   beta_i      = beta_i,
#'   theta_i     = theta_i,
#'   return_long = FALSE
#' )
#'
#' ## ---------------------------------------------------------------
#' ## Option D: fully site-varying Gamma_is and alpha_is
#' ## ---------------------------------------------------------------
#'
#' ## gamma_is : site x invader matrix of density-dependence scalars
#' ## constructed by repeating gamma_i across sites
#' Gamma_is = matrix(
#'   rep(theta_i, each = nrow(r_is_z)),
#'   nrow = nrow(r_is_z),
#'   ncol = ncol(r_is_z),
#'   dimnames = dimnames(r_is_z)
#' )
#'
#' ## alpha_is : site x invader crowding sensitivity
#' ## start from invader-specific alpha_i and allow small site-level variation
#' alpha_is = matrix(
#'   NA_real_,
#'   nrow = nrow(r_is_z),
#'   ncol = ncol(r_is_z),
#'   dimnames = dimnames(r_is_z)
#' )
#'
#' for (j in seq_len(ncol(r_is_z))) {
#'   alpha_is[, j] = alpha_i[j]
#' }
#'
#' ## Add site-level noise and enforce positivity.
#' ## pmax() drops matrix attributes, so we re-wrap explicitly.
#' alpha_is = {
#'   tmp = pmax(
#'     0,
#'     alpha_is + matrix(
#'       rnorm(length(alpha_is), 0, 0.1),
#'       nrow = nrow(r_is_z),
#'       ncol = ncol(r_is_z)
#'     )
#'   )
#'   matrix(
#'     tmp,
#'     nrow = nrow(r_is_z),
#'     ncol = ncol(r_is_z),
#'     dimnames = dimnames(r_is_z)
#'   )
#' }
#'
#' outD = compute_invasion_fitness(
#'   r_is_z, C_is_z, S_is_z,
#'   option      = "D",
#'   alpha_is    = alpha_is,
#'   beta_i      = beta_i,
#'   Gamma_is    = Gamma_is,
#'   return_long = FALSE
#' )
#'
#' ## ---------------------------------------------------------------
#' ## Option E: signed S effect (positive or negative beta_i)
#' ## ---------------------------------------------------------------
#' beta_signed_i = setNames(
#'   rnorm(I, 0, 0.3),
#'   colnames(r_is_z)
#' )
#'
#' outE = compute_invasion_fitness(
#'   r_is_z, C_is_z, S_is_z,
#'   option          = "E",
#'   alpha_i         = alpha_i,
#'   beta_signed_i   = beta_signed_i,
#'   theta0          = 1,
#'   return_long     = FALSE
#' )
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

  # Flexible aligner for invader-vectors (gamma_i, alpha_i, beta_i, etc.)
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
    if (is.null(label)) label = "Option A (gamma=1)"

  } else if (option == "B") {
    GI = matrix(theta0, nrow = S, ncol = I, dimnames = list(sites, inv_ids))
    AI = matrix(rep(.align_inv_vec(alpha_i, inv_ids), each = S), nrow = S, dimnames = list(sites, inv_ids))
    BI = .align_inv_vec(beta_i, inv_ids)
    if (is.null(label)) label = "Option B (gamma=theta_0)"

  } else if (option == "C") {
    GI = matrix(rep(.align_inv_vec(theta_i, inv_ids), each = S), nrow = S, dimnames = list(sites, inv_ids))
    AI = matrix(rep(.align_inv_vec(alpha_i, inv_ids),  each = S), nrow = S, dimnames = list(sites, inv_ids))
    BI = .align_inv_vec(beta_i, inv_ids)
    if (is.null(label)) label = "Option C (gamma=theta_i)"

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
    if (is.null(label)) label = "Option D (Gamma_is, alpha_is)"

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

  # --- Assemble lambda_is -----------------------------------------------------------
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
