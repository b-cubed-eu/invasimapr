#' Assemble site- and species-specific competition/impact matrices
#'
#' @description
#' Builds a 3D array of per-site interaction terms by combining:
#' (1) invader and resident abundances, (2) trait-based competition coefficients
#' \eqn{a_{ij}}, and (3) environmental matching between each site and residents.
#' The result is a tensor \code{I_raw[invader, resident, site]} suitable for
#' summarising total competitive pressure on each invader at each site.
#'
#' @details
#' For site \eqn{s}, invader \eqn{i}, resident \eqn{j}:
#' \deqn{I_{i j s} = r_{i s}\, r_{j s}\, a_{i j}\, K^{(env)}_{j s},}
#' where \eqn{r_{i s}} and \eqn{r_{j s}} are predicted (or expected) abundances,
#' \eqn{a_{i j} = \exp\!\big(-d_{i j}^2/(2\,\sigma_t^2)\big)} is the trait-based competition coefficient,
#' and \eqn{K^{(env)}_{j s}} is an environmental matching kernel for resident \eqn{j} at site \eqn{s}.
#'
#' @param a_ij numeric matrix (invaders × residents).
#'   Trait-based competition coefficients (e.g., from \code{compute_competition()}).
#'   Row names = invader IDs; column names = resident IDs.
#' @param Nstar numeric matrix (residents × sites).
#'   Resident abundances by site (e.g., from \code{compute_interaction_strength()}).
#'   Row names = resident IDs; column names = site IDs.
#' @param invader_pred_wide numeric matrix (sites × invaders) or \code{NULL}.
#'   Site × invader abundances on the response scale.
#'   If \code{NULL}, supply \code{predictions} instead.
#' @param predictions data.frame or \code{NULL}.
#'   Long table with columns \code{species}, \code{site_id}, \code{pred}.
#'   Used only if \code{invader_pred_wide} is \code{NULL}.
#' @param env_dist numeric matrix (sites × residents).
#'   Site–resident environmental distance (e.g., from \code{compute_environment_kernel()}).
#'   Row names = site IDs; column names = resident IDs.
#' @param sigma_e numeric or \code{NULL}.
#'   Bandwidth for the Gaussian environmental kernel;
#'   if \code{NULL}, compute the kernel upstream (e.g., \code{compute_environment_kernel(..., kernel = "gaussian")})
#'   and pass it via \code{K_env}.
#' @param K_env numeric matrix (sites × residents) or \code{NULL}.
#'   Optional precomputed environmental kernel (similarity) for residents by site.
#'   If supplied, \code{env_dist} and \code{sigma_e} are ignored.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{I_raw}: 3D array (invaders × residents × sites) of per-pair impact terms.
#'   \item \code{pressure_inv_site}: matrix (invaders × sites) with total pressure on each invader at each site
#'         (sum over residents).
#'   \item \code{meta}: list with matched IDs and dimensions.
#' }
#'
#' @examples
#' # comp$a_ij from compute_competition(), cis$Nstar from compute_interaction_strength(),
#' # ek$env_dist or ek$K_env from compute_environment_kernel(), and 'pred' from predict_invader_response()
#' am <- assemble_matrices(
#'   a_ij              = comp$a_ij,
#'   Nstar             = cis$Nstar,
#'   predictions       = pred$predictions,  # or invader_pred_wide = <sites × invaders matrix>
#'   env_dist          = ek$env_dist,
#'   sigma_e           = ek$sigma_e         # or pass K_env = ek$K_env and omit env_dist/sigma_e
#' )
#' dim(am$I_raw); dim(am$pressure_inv_site)
#'
#' @importFrom stats setNames
#' @export

assemble_matrices <- function(
    a_ij,
    Nstar,
    invader_pred_wide = NULL,
    predictions       = NULL,
    env_dist,
    sigma_e           = NULL,
    K_env             = NULL
){
  # ---- checks & alignment ----------------------------------------------------
  if (!is.matrix(a_ij)) stop("a_ij must be a matrix (invaders × residents).")
  if (!is.matrix(Nstar)) stop("Nstar must be a matrix (residents × sites).")
  inv_ids <- rownames(a_ij); res_ids <- colnames(a_ij)
  if (is.null(inv_ids) || is.null(res_ids)) stop("a_ij must have row/col names (invader/resident IDs).")
  if (!all(res_ids %in% rownames(Nstar))) stop("All residents in a_ij must appear as rows of Nstar.")
  # Reorder Nstar to match a_ij resident order
  Nstar <- Nstar[res_ids, , drop = FALSE]
  site_ids <- colnames(Nstar)
  if (is.null(site_ids)) stop("Nstar must have column names (site IDs).")

  # ---- invader abundances per site (sites × invaders) -----------------------
  if (is.null(invader_pred_wide)) {
    if (is.null(predictions)) stop("Provide 'invader_pred_wide' or 'predictions'.")
    need <- c("species","site_id","pred")
    if (!all(need %in% names(predictions))) stop("`predictions` must contain: species, site_id, pred.")
    inv_only <- predictions[predictions$species %in% inv_ids, c("species","site_id","pred"), drop = FALSE]
    # build sites × invaders matrix
    invader_pred_wide <- matrix(NA_real_, nrow = length(site_ids), ncol = length(inv_ids),
                                dimnames = list(site_ids, inv_ids))
    # fill by site
    split_by_site <- split(inv_only, inv_only$site_id)
    for (s in names(split_by_site)) {
      if (!s %in% site_ids) next
      sub <- split_by_site[[s]]
      col_idx <- match(sub$species, inv_ids, nomatch = 0L)
      ok <- col_idx > 0L
      invader_pred_wide[s, col_idx[ok]] <- as.numeric(sub$pred[ok])
    }
  } else {
    invader_pred_wide <- as.matrix(invader_pred_wide)
    # align rows/cols
    if (is.null(rownames(invader_pred_wide)) || is.null(colnames(invader_pred_wide)))
      stop("invader_pred_wide must have rownames (sites) and colnames (invaders).")
    # reorder to match site_ids and inv_ids
    invader_pred_wide <- invader_pred_wide[site_ids, inv_ids, drop = FALSE]
  }

  # ---- environmental kernel per site × resident -----------------------------
  if (is.null(K_env)) {
    if (is.null(env_dist) || is.null(sigma_e))
      stop("Provide either K_env, or env_dist + sigma_e for Gaussian env kernel.")
    # align env_dist
    if (is.null(rownames(env_dist)) || is.null(colnames(env_dist)))
      stop("env_dist must have rownames (sites) and colnames (residents).")
    env_dist <- as.matrix(env_dist)[site_ids, res_ids, drop = FALSE]
    K_env <- exp(-(env_dist^2)/(2 * sigma_e^2))
  } else {
    # align K_env
    if (is.null(rownames(K_env)) || is.null(colnames(K_env)))
      stop("K_env must have rownames (sites) and colnames (residents).")
    K_env <- as.matrix(K_env)[site_ids, res_ids, drop = FALSE]
  }

  # ---- build 3D tensor: [invader, resident, site] ---------------------------
  n_inv <- length(inv_ids); n_res <- length(res_ids); n_sites <- length(site_ids)
  I_raw <- array(0, dim = c(n_inv, n_res, n_sites),
                 dimnames = list(inv_ids, res_ids, site_ids))

  # Precompute resident abundance & env kernel slices by site
  # Loop over sites to keep memory predictable and clear
  for (s in seq_len(n_sites)) {
    sid <- site_ids[s]
    r_i <- invader_pred_wide[sid, , drop = TRUE]  # length n_inv
    r_j <- Nstar[ , sid, drop = TRUE]             # length n_res
    # outer product of abundances
    A <- outer(r_i, r_j)                          # n_inv × n_res
    # environmental kernel for residents at site s (replicated down rows)
    Evec <- K_env[sid, , drop = TRUE]             # length n_res
    Eker <- matrix(Evec, nrow = n_inv, ncol = n_res, byrow = TRUE)
    # combine with competition coefficients
    I_raw[ , , s] <- A * a_ij * Eker
  }

  # ---- summaries: total pressure on invaders per site ------------------------
  # Sum over residents (dimension 2)
  pressure_inv_site <- apply(I_raw, c(1,3), sum, na.rm = TRUE) # (invaders × sites)

  list(
    I_raw              = I_raw,
    pressure_inv_site  = pressure_inv_site,
    meta               = list(
      n_inv = n_inv, n_res = n_res, n_sites = n_sites,
      invaders = inv_ids, residents = res_ids, sites = site_ids
    )
  )
}
