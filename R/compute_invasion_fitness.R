#' Compute invasion fitness from competition/impact tensors and growth predictions
#'
#' @description
#' Aggregates the site- and species-specific impact tensor to quantify the **total
#' competitive penalty** experienced by each invader at each site, then computes
#' several invasion-fitness formulations by subtracting competition penalties from
#' predicted growth.
#'
#' @details
#' Let \eqn{I_{i j s}} be the per-pair impact at site \eqn{s} (e.g., from
#' \code{assemble_matrices()}), with invader \eqn{i}, resident \eqn{j}. The total
#' competition penalty at each invader–site is
#' \deqn{C_{i s} = \sum_{j} I_{i j s}.}
#' Given a matrix of invader growth predictions \eqn{r_{i s}}, we compute:
#' \itemize{
#'   \item \strong{Raw:} \eqn{\lambda^{raw}_{i s} = r_{i s} - C_{i s}}.
#'   \item \strong{Scaled:} \eqn{\lambda^{scaled}_{i s} = r_{i s} - C_{i s}/J} \,(where \eqn{J} is the number of residents).
#'   \item \strong{Relative-abundance:} \eqn{\lambda^{rel}_{i s} = r_{i s} - (A N^{rel})_{i s}}, where
#'         \eqn{A = [a_{i j}]} and \eqn{N^{rel}_{j s} = N_{j s} / \sum_{j'} N_{j' s}}.
#'   \item \strong{Logistic-capped:} \eqn{\lambda^{logis}_{i s} = r_{i s} - \frac{1}{1 + \exp\{-k\, (C^{*}_{i s} - x_{0})\}}},
#'         with \eqn{C^{*}_{i s}} taken as either \eqn{C_{i s}} (raw) or \eqn{(A N^{rel})_{i s}} (relative),
#'         selected via \code{logistic_on = "raw"} or \code{"rel"}.
#' }
#'
#' @param I_raw 3D numeric array \code{[invader, resident, site]} from \code{assemble_matrices()},
#'   or \code{NULL} if you pass \code{pressure_inv_site} instead.
#' @param pressure_inv_site numeric matrix \code{[invader, site]} (sum over residents),
#'   optional alternative to \code{I_raw}.
#' @param r_mat numeric matrix \code{[invader, site]} of invader growth on the response scale,
#'   or \code{NULL} if you pass \code{predictions}.
#' @param predictions data.frame or \code{NULL}. Long table with columns
#'   \code{species}, \code{site_id}, \code{pred}. Used to build \code{r_mat} if \code{r_mat} is \code{NULL}.
#' @param a_ij numeric matrix \code{[invader, resident]} of competition coefficients; required
#'   only for \code{lambda_rel} or for \code{lambda_logis} when \code{logistic_on = "rel"}.
#' @param Nstar numeric matrix \code{[resident, site]} of resident abundances; required for \code{lambda_rel}
#'   and for \code{lambda_logis} when \code{logistic_on = "rel"}.
#' @param logistic_on character. Use \code{"raw"} (default) to apply the logistic cap to \code{C_raw};
#'   use \code{"rel"} to apply it to \code{A \%*\% N_rel}.
#' @param k numeric. Logistic steepness parameter (default \code{1}).
#' @param x0 numeric or \code{NULL}. Logistic midpoint. If \code{NULL}, set to
#'   \code{median(C_target, na.rm = TRUE)}.
#' @param prefer \code{"logis"}|\code{"rel"}|\code{"raw"}|\code{"scaled"}. Which fitness to return as \code{$lambda}
#'   (default \code{"logis"}).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{C_raw} \code{[invader, site]}: total penalty summed over residents.
#'   \item \code{r_mat} \code{[invader, site]}: invader growth matrix.
#'   \item \code{lambda_raw}, \code{lambda_scaled}, \code{lambda_rel}, \code{lambda_logis}: fitness variants.
#'   \item \code{lambda}: the selected fitness (per \code{prefer}).
#'   \item \code{meta}: data about inputs, dimensions, and chosen options.
#' }
#'
#' @examples
#' fit <- compute_invasion_fitness(
#'   I_raw       = am$I_raw,                # or pressure_inv_site = am$pressure_inv_site
#'   predictions = pred$predictions,        # constructs r_mat internally
#'   a_ij        = comp$a_ij,               # for lambda_rel / logistic_on = "rel"
#'   Nstar       = cis$Nstar,
#'   logistic_on = "rel", k = 1, x0 = NULL, # logistic on A %*% N_rel
#'   prefer      = "logis"
#' )
#' str(fit$lambda)  # final fitness matrix [invader, site]
#'
#' @importFrom stats median
#' @export
compute_invasion_fitness <- function(
    I_raw               = NULL,
    pressure_inv_site   = NULL,
    r_mat               = NULL,
    predictions         = NULL,
    a_ij                = NULL,
    Nstar               = NULL,
    logistic_on         = c("raw","rel"),
    k                   = 1,
    x0                  = NULL,
    prefer              = c("logis","rel","raw","scaled")
){
  logistic_on <- match.arg(logistic_on)
  prefer      <- match.arg(prefer)

  # --- 1) Build C_raw (invader × site): sum over residents -------------------
  if (!is.null(I_raw)) {
    if (length(dim(I_raw)) != 3) stop("I_raw must be a 3D array [invader, resident, site].")
    # sum over residents (dim 2)
    C_raw <- apply(I_raw, c(1,3), sum, na.rm = TRUE)
    inv_ids  <- dimnames(I_raw)[[1]]; site_ids <- dimnames(I_raw)[[3]]
  } else {
    if (is.null(pressure_inv_site)) stop("Provide I_raw or pressure_inv_site.")
    C_raw <- as.matrix(pressure_inv_site)
    inv_ids  <- rownames(C_raw); site_ids <- colnames(C_raw)
  }
  if (is.null(inv_ids) || is.null(site_ids)) stop("C_raw needs row/colnames (invaders/sites).")

  # --- 2) Build r_mat (invader × site) ---------------------------------------
  if (is.null(r_mat)) {
    if (is.null(predictions)) stop("Provide r_mat or predictions to construct r_mat.")
    need <- c("species","site_id","pred")
    if (!all(need %in% names(predictions))) stop("predictions must contain: species, site_id, pred.")
    # build site × invader, then transpose to invader × site ordered by C_raw
    inv_set <- unique(inv_ids)
    site_set<- unique(site_ids)
    r_wide <- matrix(NA_real_, nrow = length(site_set), ncol = length(inv_set),
                     dimnames = list(site_set, inv_set))
    split_by_site <- split(predictions, predictions$site_id)
    for (s in names(split_by_site)) {
      if (!s %in% site_set) next
      sub <- split_by_site[[s]]
      col_idx <- match(sub$species, inv_set, nomatch = 0L)
      ok <- col_idx > 0L
      r_wide[s, col_idx[ok]] <- as.numeric(sub$pred[ok])
    }
    r_mat <- t(r_wide)[inv_ids, site_ids, drop = FALSE]
  } else {
    r_mat <- as.matrix(r_mat)[inv_ids, site_ids, drop = FALSE]
  }

  # --- 3) Auxiliary for relative-abundance penalty ---------------------------
  lambda_rel <- lambda_logis <- NULL
  if (!is.null(a_ij) && !is.null(Nstar)) {
    a_ij <- as.matrix(a_ij)
    # align dimnames
    if (is.null(rownames(a_ij)) || is.null(colnames(a_ij))) stop("a_ij must have row (invaders) and col (residents) names.")
    if (!all(rownames(a_ij) %in% inv_ids)) stop("Some invaders in a_ij not found in C_raw/r_mat.")
    # reorder to inv_ids
    a_ij <- a_ij[inv_ids, , drop = FALSE]
    # Nstar: residents × sites, align to a_ij cols and site_ids
    if (is.null(rownames(Nstar)) || is.null(colnames(Nstar))) stop("Nstar must have row (residents) and col (sites) names.")
    if (!all(colnames(a_ij) %in% rownames(Nstar))) stop("Residents in a_ij must be rows in Nstar.")
    Nstar_use <- Nstar[colnames(a_ij), site_ids, drop = FALSE]
    # normalise by site totals
    N_rel <- sweep(Nstar_use, 2, colSums(Nstar_use, na.rm = TRUE), "/")
    # a_ij %*% N_rel -> invader × site
    C_rel <- a_ij %*% N_rel
  } else {
    C_rel <- NULL
  }

  # --- 4) Fitness variants ----------------------------------------------------
  J <- if (!is.null(a_ij)) ncol(a_ij) else NA_integer_
  lambda_raw    <- r_mat - C_raw
  lambda_scaled <- r_mat - C_raw / ifelse(is.na(J), nrow(C_raw)*0 + nrow(C_raw), J)  # fallback J if missing

  if (!is.null(C_rel)) {
    lambda_rel <- r_mat - C_rel
  }

  # Logistic-capped (on chosen target)
  C_target <- if (logistic_on == "rel" && !is.null(C_rel)) C_rel else C_raw
  if (is.null(x0)) x0 <- stats::median(C_target, na.rm = TRUE)
  pen <- 1 / (1 + exp(-k * (C_target - x0)))
  lambda_logis <- r_mat - pen

  # --- 5) Pick preferred output ----------------------------------------------
  lambda <- switch(
    prefer,
    logis  = lambda_logis,
    rel    = { if (is.null(lambda_rel)) stop("lambda_rel requested but a_ij/Nstar not supplied."); lambda_rel },
    raw    = lambda_raw,
    scaled = lambda_scaled
  )

  list(
    C_raw         = C_raw,
    r_mat         = r_mat,
    lambda_raw    = lambda_raw,
    lambda_scaled = lambda_scaled,
    lambda_rel    = lambda_rel,
    lambda_logis  = lambda_logis,
    lambda        = lambda,
    meta          = list(
      logistic_on = logistic_on, k = k, x0 = x0,
      prefer = prefer,
      dims = list(invaders = length(inv_ids), sites = length(site_ids),
                  residents = if (!is.null(a_ij)) ncol(a_ij) else NA_integer_)
    )
  )
}
