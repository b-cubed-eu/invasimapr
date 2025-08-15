#' Environmental optima and site-species environmental distance (with optional kernel)
#'
#' @description
#' Estimates each species’ **environmental optimum** as the abundance-weighted mean of
#' site-level environmental covariates, then computes the **distance** (default: Gower)
#' between every **site** and every **species optimum**. Optionally returns an
#' **environmental kernel** (similarity) by transforming distances (e.g., Gaussian).
#'
#' @details
#' Let \eqn{E_s} be the vector of environmental covariates at site \eqn{s}, and
#' \eqn{w_{s,j}} the abundance (weight) of species \eqn{j} at site \eqn{s}. The
#' environmental **optimum** for species \eqn{j} is
#' \deqn{ \mu_j \;=\; \frac{\sum_s w_{s,j} \, E_s}{\sum_s w_{s,j}}. }
#' The **environmental distance** between site \eqn{s} and species \eqn{j} is then
#' \eqn{d_{s j} = \mathrm{Gower}(E_s, \mu_j)} (or another supported metric).
#'
#' The optional **kernel** converts distances to similarity, e.g., Gaussian
#' \eqn{K_{sj} = \exp\{-d_{sj}^2/(2 \sigma_e^2)\}}. The bandwidth \eqn{\sigma_e}
#' controls how quickly suitability decays with mismatch; by default we estimate
#' \eqn{\sigma_e} from the distribution of \eqn{d_{sj}}.
#'
#' @param site_env data.frame. One row per site; must include `site_col` and the
#'   environmental variables (numeric/factor/ordered). Coordinate columns can be present.
#' @param abundance_wide matrix/data.frame or NULL. **Sites × species** numeric
#'   abundance/weight table (rows = sites, columns = species). If `NULL`, provide
#'   `predictions` instead.
#' @param predictions data.frame or NULL. Long table with columns `species`,
#'   `site_id`, and `pred` to construct `abundance_wide` internally.
#' @param site_col character. Site identifier column in `site_env` (and `predictions`). Default `"site_id"`.
#' @param env_cols NULL or character. Environmental columns in `site_env`. If `NULL`,
#'   auto-detect as all non-`site_col` non-coordinate columns.
#' @param coord_cols character. Columns to exclude from env detection (e.g., coords). Default `c("x","y")`.
#' @param method character. Distance metric for site-optimum comparison: `"gower"` (default),
#'   `"euclidean"`, or `"manhattan"`.
#' @param gower_stand logical. Standardise numeric vars inside Gower (default `TRUE`).
#' @param kernel character. One of `"distance"` (no transform), `"similarity"` (1 - scaled distance),
#'   or `"gaussian"` (Gaussian kernel). Default `"distance"`.
#' @param sigma_e numeric or NULL. Bandwidth for Gaussian kernel. If `NULL`, estimated from
#'   non-zero `env_dist` using `sigma_method`.
#' @param sigma_method character. How to estimate `sigma_e` when missing: `"sd"` (default),
#'   `"median"`, or `"iqr"`.
#' @param scale_01 logical. Min-max scale distances to [0,1] **before** `"similarity"`;
#'   for non-Gower metrics only, this is applied automatically if needed.
#'
#' @return A list with:
#' \itemize{
#'   \item `env_opt`: **species × env** matrix of abundance-weighted optima.
#'   \item `env_dist`: **sites × species** distance matrix (site ↔ species-optimum).
#'   \item `K_env`: **sites × species** kernel (if `kernel != "distance"`).
#'   \item `sigma_e`: bandwidth used (if Gaussian).
#'   \item `meta`: list of settings and detected columns.
#' }
#'
#' @examples
#' # --- Simulate sites, environments, and species weights (no extra packages) ---
#' set.seed(123)
#'
#' # Sites and species
#' n_sites <- 6
#' n_spp   <- 4
#' sites   <- paste0("s", seq_len(n_sites))
#' spp     <- paste0("sp", seq_len(n_spp))
#'
#' # Site-level environment table (include optional coord columns)
#' site_env <- data.frame(
#'   site_id = sites,
#'   x       = rnorm(n_sites),
#'   y       = rnorm(n_sites),
#'   temp    = runif(n_sites, 5, 25),     # numeric env var
#'   precip  = runif(n_sites, 400, 900),  # numeric env var
#'   check.names = FALSE
#' )
#'
#' # Sites × species abundance/weights (any non-negative numbers)
#' abundance_wide <- matrix(
#'   rexp(n_sites * n_spp, rate = 1),
#'   nrow = n_sites, ncol = n_spp,
#'   dimnames = list(sites, spp)
#' )
#'
#' # Compute environmental optima and site–species distances
#' ek <- compute_environment_kernel(
#'   site_env       = site_env,
#'   abundance_wide = abundance_wide,   # avoids needing tidyr
#'   site_col       = "site_id",
#'   method         = "euclidean",      # uses stats::dist (no extra deps)
#'   kernel         = "gaussian",       # also returns K_env
#'   sigma_method   = "sd"
#' )
#'
#' # Inspect results
#' str(ek$env_opt)       # species × env optima
#' dim(ek$env_dist)      # sites × species distance matrix
#' if (!is.null(ek$K_env)) range(ek$K_env, na.rm = TRUE)
#'
#' # (Optional) If you prefer Gower distances, set method = "gower".
#' # This uses cluster::daisy internally and may require the 'cluster' package:
#' # ek_gower <- compute_environment_kernel(
#' #   site_env       = site_env,
#' #   abundance_wide = abundance_wide,
#' #   site_col       = "site_id",
#' #   method         = "gower",
#' #   kernel         = "similarity"
#' # )
#'
#' @importFrom stats dist sd median IQR
#' @export
compute_environment_kernel <- function(
    site_env,
    abundance_wide = NULL,
    predictions = NULL,
    site_col = "site_id",
    env_cols = NULL,
    coord_cols = c("x", "y"),
    method = c("gower", "euclidean", "manhattan"),
    gower_stand = TRUE,
    kernel = c("distance", "similarity", "gaussian"),
    sigma_e = NULL,
    sigma_method = c("sd", "median", "iqr"),
    scale_01 = TRUE) {
  # ---- checks ----------------------------------------------------------------
  stopifnot(is.data.frame(site_env))
  if (!site_col %in% names(site_env)) stop("`site_col` not found in `site_env`.")
  method <- match.arg(method)
  kernel <- match.arg(kernel)
  sigma_method <- match.arg(sigma_method)

  # ---- detect environmental columns ------------------------------------------
  if (is.null(env_cols)) {
    env_cols <- setdiff(names(site_env), c(site_col, coord_cols))
  }
  env_cols <- intersect(env_cols, names(site_env))
  if (!length(env_cols)) stop("No environmental columns detected in `site_env`.")

  # ---- build abundance_wide if needed (sites × species) ----------------------
  if (is.null(abundance_wide)) {
    if (is.null(predictions)) stop("Provide either `abundance_wide` or `predictions`.")
    need <- c("species", site_col, "pred")
    if (!all(need %in% names(predictions))) {
      stop("`predictions` must contain: ", paste(need, collapse = ", "))
    }

    # pivot to sites × species
    if (requireNamespace("tidyr", quietly = TRUE)) {
      aw <- tidyr::pivot_wider(
        predictions,
        id_cols = tidyselect::all_of(site_col),
        names_from = "species",
        values_from = "pred"
      )
      rn <- aw[[site_col]]
      aw[[site_col]] <- NULL
      abundance_wide <- as.matrix(aw)
      rownames(abundance_wide) <- as.character(rn)
    } else {
      # base fallback
      sp <- sort(unique(as.character(predictions$species)))
      si <- sort(unique(as.character(predictions[[site_col]])))
      abundance_wide <- matrix(NA_real_,
        nrow = length(si), ncol = length(sp),
        dimnames = list(si, sp)
      )
      for (i in seq_len(nrow(predictions))) {
        s <- as.character(predictions[[site_col]][i])
        sp_i <- as.character(predictions$species[i])
        abundance_wide[s, sp_i] <- predictions$pred[i]
      }
    }
  } else {
    abundance_wide <- as.matrix(abundance_wide)
  }

  # ---- align sites between site_env and abundance ----------------------------
  # Keep only sites present in both
  common_sites <- intersect(as.character(site_env[[site_col]]), rownames(abundance_wide))
  if (!length(common_sites)) stop("No overlapping sites between `site_env` and `abundance_wide`.")
  site_env_use <- site_env[match(common_sites, site_env[[site_col]]), c(site_col, env_cols), drop = FALSE]
  abundance_wide <- abundance_wide[common_sites, , drop = FALSE]

  # ---- compute species environmental optima (weighted means over sites) ------
  E <- as.matrix(site_env_use[, env_cols, drop = FALSE]) # S × P (env)
  W <- as.matrix(abundance_wide) # S × J (weights)
  # Handle species with zero total weight: return NA row
  col_sums_W <- colSums(W, na.rm = TRUE)
  env_opt <- t(vapply(
    seq_len(ncol(W)),
    function(j) {
      w <- W[, j]
      tot <- sum(w, na.rm = TRUE)
      if (is.na(tot) || tot == 0) rep(NA_real_, ncol(E)) else colSums(E * w, na.rm = TRUE) / tot
    },
    numeric(ncol(E))
  ))
  rownames(env_opt) <- colnames(W) # species
  colnames(env_opt) <- colnames(E) # env vars

  # ---- build combined table for site-optimum distances -----------------------
  # rows: sites then species optima; columns: env variables only
  site_block <- E
  rownames(site_block) <- site_env_use[[site_col]]
  env_all <- rbind(site_block, env_opt)

  # ---- compute distances (site ↔ species-optimum) ----------------------------
  if (method == "gower") {
    d_all <- cluster::daisy(env_all, metric = "gower", stand = isTRUE(gower_stand))
    D <- as.matrix(d_all)
  } else {
    # numeric-only fallback for Euclidean/Manhattan
    num_only <- vapply(as.data.frame(env_all), is.numeric, logical(1))
    if (!all(num_only)) {
      warning(
        "Non-numeric env columns dropped for method = '", method, "': ",
        paste(colnames(env_all)[!num_only], collapse = ", ")
      )
      env_all <- as.matrix(env_all[, num_only, drop = FALSE])
    } else {
      env_all <- as.matrix(env_all)
    }
    D <- as.matrix(stats::dist(env_all, method = method))
  }

  # extract site (rows) × species (cols) block
  site_ids <- rownames(site_block)
  spp_ids <- rownames(env_opt)
  env_dist <- D[site_ids, spp_ids, drop = FALSE]

  # ---- optional scaling for similarity kernel --------------------------------
  # Ensure [0,1] for non-Gower if needed (safe for Gower too when scale_01 = TRUE)
  if ((kernel == "similarity" && isTRUE(scale_01)) ||
    (kernel == "similarity" && method != "gower")) {
    v <- env_dist[is.finite(env_dist)]
    rng <- range(v, na.rm = TRUE)
    if (diff(rng) > 0) env_dist <- (env_dist - rng[1]) / diff(rng) else env_dist[, ] <- 0
  }

  # ---- kernelisation ---------------------------------------------------------
  K_env <- NULL
  if (kernel == "similarity") {
    # similarity in [0,1]
    K_env <- 1 - env_dist
  } else if (kernel == "gaussian") {
    if (is.null(sigma_e)) {
      vals <- env_dist[is.finite(env_dist) & env_dist > 0]
      sigma_e <- if (length(vals)) {
        switch(sigma_method,
          sd     = stats::sd(vals),
          median = stats::median(vals),
          iqr    = stats::IQR(vals)
        )
      } else {
        1
      }
      if (!is.finite(sigma_e) || sigma_e <= 0) sigma_e <- 1
    }
    K_env <- exp(-(env_dist^2) / (2 * sigma_e^2))
  }

  # ---- return ----------------------------------------------------------------
  list(
    env_opt = env_opt, # species × env
    env_dist = env_dist, # sites × species (distance)
    K_env = K_env, # sites × species (similarity), if requested
    sigma_e = if (!is.null(K_env) && kernel == "gaussian") sigma_e else NULL,
    meta = list(
      site_col = site_col,
      env_cols = env_cols,
      method = method,
      gower_stand = gower_stand,
      kernel = kernel,
      sigma_method = if (is.null(sigma_e)) NA_character_ else sigma_method
    )
  )
}
