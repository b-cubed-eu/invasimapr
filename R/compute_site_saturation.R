#' Compute site-only saturation S_s and global z-score S_s_z
#'
#' @description
#' Implements three site-only saturation options: "opportunity_penalty",
#' "modelled_dominance", and "evenness_deficit". Coerces non-matrix inputs with
#' `as.matrix()`. Returns S_s, global z, and a broadcast S_js_z for convenience.
#'
#' @param mode One of "opportunity_penalty", "modelled_dominance", "evenness_deficit".
#' @param comm_res Site × resident matrix (or data.frame/tibble) of counts for the
#'   observed-based modes (will be coerced to matrix if needed).
#' @param res_ids Character vector of resident IDs (columns of \code{comm_res}).
#' @param mu_js Optional site × resident matrix (or data.frame/tibble) of expected
#'   abundances for "modelled_dominance" (will be coerced to matrix if needed).
#' @return List with \code{S_s}, \code{S_mu}, \code{S_sd}, \code{S_s_z}, and \code{S_js_z}.
#' @examples
#' \dontrun{
#' sat = compute_site_saturation("evenness_deficit", comm_res, colnames(comm_res))
#' S_js_z = sat$S_js_z
#' }
#' @importFrom vegan diversity specnumber
#' @export
compute_site_saturation = function(mode,
                                    comm_res,
                                    res_ids,
                                    mu_js = NULL) {
  # -- Coerce comm_res to matrix if needed (e.g., tibble/data.frame) ------------
  if (!is.matrix(comm_res)) comm_res = as.matrix(comm_res)

  # Ensure we have row names (needed for broadcasting); fabricate if missing
  if (is.null(rownames(comm_res))) {
    rownames(comm_res) = as.character(seq_len(nrow(comm_res)))
  }

  # Ensure numeric storage (vegan expects numeric); this may turn non-numeric into NA
  suppressWarnings(storage.mode(comm_res) == "double")

  # Check that requested resident columns exist; subset, but warn on any missings
  missing_cols = setdiff(res_ids, colnames(comm_res))
  if (length(missing_cols)) {
    warning("compute_site_saturation(): dropping ", length(missing_cols),
            " res_ids not found in comm_res: ",
            paste(head(missing_cols, 5), collapse = ", "),
            if (length(missing_cols) > 5) " ..." else "")
  }
  keep_cols = intersect(res_ids, colnames(comm_res))
  if (!length(keep_cols)) stop("No overlapping resident IDs between res_ids and comm_res.")
  X = comm_res[, keep_cols, drop = FALSE]

  sites = rownames(comm_res)

  # -- Handle mode-specific definitions ----------------------------------------
  if (mode == "opportunity_penalty") {
    # Higher total resident abundance → higher saturation (penalty)
    # (S_s is increasing in summed counts; log-scaled)
    S_s = log1p(rowSums(X, na.rm = TRUE))

  } else if (mode == "modelled_dominance") {
    # Use expected abundances mu_js instead of observed X
    if (is.null(mu_js)) stop("`mu_js` must be provided for mode = 'modelled_dominance'.")
    if (!is.matrix(mu_js)) mu_js = as.matrix(mu_js)
    # Align rows to sites if possible
    if (!is.null(rownames(mu_js)) && !is.null(sites) && !all(sites %in% rownames(mu_js))) {
      warning("Some sites in comm_res not found in mu_js; aligning by intersection.")
    }
    # Reorder/align rows to 'sites' where possible; otherwise use current order
    if (!is.null(rownames(mu_js))) {
      mu_js = mu_js[intersect(sites, rownames(mu_js)), , drop = FALSE]
      # if alignment changed, also subset sites accordingly
      sites = rownames(mu_js)
      # Re-broadcast base X to aligned sites for consistent dimnames later
      X = X[sites, , drop = FALSE]
    }
    suppressWarnings(storage.mode(mu_js) == "double")
    S_s = log1p(rowSums(mu_js, na.rm = TRUE))

  } else if (mode == "evenness_deficit") {
    # Pielou’s evenness deficit: 1 - J, where J = H / log(S)
    H = vegan::diversity(X, index = "shannon", MARGIN = 1)
    S = vegan::specnumber(X, MARGIN = 1)
    J = H / log(pmax(S, 2))   # guard small S; uses log(2) when S=1
    J[!is.finite(J)] = 0
    S_s = 1 - J

  } else {
    stop("Unknown mode for site saturation: '", mode, "'.")
  }

  # -- Global z-scoring (site-only) and broadcast to site×resident -------------
  S_mu = mean(S_s, na.rm = TRUE)
  S_sd = stats::sd(S_s, na.rm = TRUE)
  if (!is.finite(S_sd) || S_sd == 0) S_sd = 1
  S_s_z = (S_s - S_mu) / S_sd

  # Broadcast S_s_z to each resident column (site × resident)
  S_js_z = matrix(S_s_z,
                   nrow = length(sites),
                   ncol = length(res_ids),
                   dimnames = list(sites, res_ids))

  list(S_s = S_s, S_mu = S_mu, S_sd = S_sd, S_s_z = S_s_z, S_js_z = S_js_z)
}

# -- small internal helper to fetch coefficients safely -------------------------
.get_cf = function(cf, main, by = NULL, default = 0) {
  if (is.null(by)) return(if (main %in% names(cf)) unname(cf[main]) else default)
  nm = c(paste0(main, ":", by), paste0(by, ":", main))
  nm = nm[nm %in% names(cf)]
  if (length(nm)) unname(cf[nm[1]]) else default
}
