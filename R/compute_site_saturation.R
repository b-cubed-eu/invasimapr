#' Compute site-only saturation (\eqn{S_s}) and global z-score (\eqn{S_s^{(z)}})
#'
#' Implements three site-only saturation definitions and returns (i) the raw
#' site statistic \eqn{S_s}, (ii) its global mean/SD, (iii) the global
#' z-scored site statistic \eqn{S_s^{(z)}}, and (iv) a **broadcast**
#' site–resident matrix \eqn{S_{js}^{(z)}} for downstream modelling.
#'
#' @section Definitions:
#' Let \(X_{js}\) be observed resident counts and \(\mu_{js}\) expected counts.
#' \describe{
#'   \item{\code{"opportunity_penalty"}}{
#'     Abundance penalty: \(\displaystyle S_s=\log(1+\sum_j X_{js})\).
#'   }
#'   \item{\code{"modelled_dominance"}}{
#'     Model-based crowding: \(\displaystyle S_s=\log(1+\sum_j \mu_{js})\).
#'     Requires `mu_js`.
#'   }
#'   \item{\code{"evenness_deficit"}}{
#'     Pielou’s evenness deficit:
#'     \(J_s = H_s / \log S_s^{(\#)}\) with Shannon \(H_s\) and richness
#'     \(S_s^{(\#)}\); then \(\displaystyle S_s = 1 - J_s\).
#'   }
#' }
#'
#' After computing \(S_s\), a **global** z-score is applied:
#' \(\displaystyle S_s^{(z)} = (S_s - \bar S)/\mathrm{sd}(S_s)\),
#' and this is broadcast across residents to form \(S_{js}^{(z)}\) (same value
#' for all \(j\) within a site \(s\)).
#'
#' @param mode One of \code{c("opportunity_penalty","modelled_dominance","evenness_deficit")}.
#' @param comm_res Site × resident matrix (or frame) of observed counts. Will be
#'   coerced to a numeric matrix; non-numeric entries become \code{NA}.
#' @param res_ids Character vector of resident IDs (column names to consider).
#'   Columns not found in \code{comm_res} are ignored with a warning.
#' @param mu_js Optional site × resident matrix (or frame) of expected abundances,
#'   required for \code{mode = "modelled_dominance"}. Will be coerced to numeric matrix.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{S_s}}{Numeric vector (sites) of saturation values.}
#'   \item{\code{S_mu}}{Global mean of \code{S_s}.}
#'   \item{\code{S_sd}}{Global SD of \code{S_s} (set to 1 if 0/NA).}
#'   \item{\code{S_s_z}}{Global z-scores of \code{S_s}.}
#'   \item{\code{S_js_z}}{Site × resident matrix broadcasting \code{S_s_z} across \code{res_ids}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Evenness deficit from observed community
#' sat <- compute_site_saturation("evenness_deficit", comm_res, colnames(comm_res))
#' head(sat$S_s_z)
#'
#' # Modelled dominance using expected abundances mu_js
#' sat2 <- compute_site_saturation("modelled_dominance", comm_res, colnames(comm_res), mu_js = mu_js)
#' dim(sat2$S_js_z)
#' }
#'
#' @seealso
#' \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_site_saturation.html}{Package reference: compute_site_saturation()}
#'
#' @importFrom vegan diversity specnumber
#' @export
compute_site_saturation = function(mode,
                                   comm_res,
                                   res_ids,
                                   mu_js = NULL) {
  # ---- Validate mode ----------------------------------------------------------
  mode <- match.arg(mode, c("opportunity_penalty", "modelled_dominance", "evenness_deficit"))

  # ---- Coerce and sanitise observed community matrix --------------------------
  if (!is.matrix(comm_res)) comm_res <- as.matrix(comm_res)
  if (is.null(rownames(comm_res))) rownames(comm_res) <- as.character(seq_len(nrow(comm_res)))

  # Ensure numeric storage (vegan expects numeric); non-numeric -> NA
  suppressWarnings(storage.mode(comm_res) <- "double")
  sites <- rownames(comm_res)

  # Restrict to requested resident columns (warn on misses)
  missing_cols <- setdiff(res_ids, colnames(comm_res))
  if (length(missing_cols)) {
    warning("compute_site_saturation(): dropping ", length(missing_cols),
            " res_ids not found in comm_res: ",
            paste(head(missing_cols, 5), collapse = ", "),
            if (length(missing_cols) > 5) " ..." else "")
  }
  keep_cols <- intersect(res_ids, colnames(comm_res))
  if (!length(keep_cols)) stop("No overlapping resident IDs between res_ids and comm_res.")
  X <- comm_res[, keep_cols, drop = FALSE]

  # ---- Mode-specific definitions ----------------------------------------------
  if (mode == "opportunity_penalty") {

    # Higher total resident abundance => larger S_s (penalty). Log(1 + sum) is
    # stable for zeros and reduces scale skew.
    S_s <- log1p(rowSums(X, na.rm = TRUE))

  } else if (mode == "modelled_dominance") {

    # Use expected abundances mu_js instead of observed X.
    if (is.null(mu_js)) stop("`mu_js` must be provided for mode = 'modelled_dominance'.")
    if (!is.matrix(mu_js)) mu_js <- as.matrix(mu_js)
    suppressWarnings(storage.mode(mu_js) <- "double")

    # Align mu_js rows to site order of comm_res when rownames are available.
    if (!is.null(rownames(mu_js))) {
      if (!all(sites %in% rownames(mu_js))) {
        warning("Some sites in comm_res not found in mu_js; aligning by rownames and allowing NAs.")
      }
      mu_js <- mu_js[sites, , drop = FALSE]
    }
    # S_s from expected totals
    S_s <- log1p(rowSums(mu_js, na.rm = TRUE))

  } else if (mode == "evenness_deficit") {

    # Pielou’s evenness deficit: 1 - J, with J = H / log(S#).
    # H: Shannon (base e); S#: richness. Guard S#<2 to avoid log(1)=0.
    H <- vegan::diversity(X, index = "shannon", MARGIN = 1)
    Snum <- vegan::specnumber(X, MARGIN = 1)
    J <- H / log(pmax(Snum, 2L))
    J[!is.finite(J)] <- 0
    S_s <- 1 - J

  } else {
    stop("Unhandled mode: ", mode)  # should be unreachable due to match.arg
  }

  # ---- Global z-score on site statistic and broadcast to site×resident --------
  S_mu <- mean(S_s, na.rm = TRUE)
  S_sd <- stats::sd(S_s, na.rm = TRUE)
  if (!is.finite(S_sd) || S_sd == 0) S_sd <- 1
  S_s_z <- (S_s - S_mu) / S_sd

  # Broadcast S_s_z across residents j for each site s
  S_js_z <- matrix(S_s_z,
                   nrow = length(sites),
                   ncol = length(res_ids),
                   dimnames = list(sites, res_ids))

  list(S_s = S_s, S_mu = S_mu, S_sd = S_sd, S_s_z = S_s_z, S_js_z = S_js_z)
}

# --- Small internal helper to fetch coefficients safely ------------------------
.get_cf = function(cf, main, by = NULL, default = 0) {
  if (is.null(by)) return(if (main %in% names(cf)) unname(cf[main]) else default)
  nm = c(paste0(main, ":", by), paste0(by, ":", main))
  nm = nm[nm %in% names(cf)]
  if (length(nm)) unname(cf[nm[1]]) else default
}
