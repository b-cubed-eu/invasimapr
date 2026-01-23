#' Compute resident crowding from community composition and trait similarity
#'
#' @title Resident crowding \eqn{C_{js}} via Hellinger composition and Gower–Gaussian trait kernel
#'
#' @description
#' Given a site × resident community matrix and a resident traits table, this function:
#' \enumerate{
#'   \item normalizes site composition with Hellinger transformation (\eqn{W_{s\cdot}}),
#'   \item computes Gower distances among residents (mixed types supported),
#'   \item converts distances to a Gaussian similarity kernel \eqn{K} with bandwidth \eqn{\sigma_\alpha},
#'   \item returns the resident crowding matrix \eqn{C = W \times K^\top} (sites × residents),
#'   \item optionally row-z-scores \eqn{C} for within-site contrasts,
#'   \item optionally draws a heatmap of \eqn{C} (row-z) and a geographic map of site-level summaries.
#' }
#'
#' @param comm_res A numeric matrix or data.frame of \strong{sites × resident species}
#'   (non-negative abundances or transformed counts). Row names must be site IDs.
#' @param traits_res A data.frame of \strong{resident species × traits}. Row names must be
#'   resident IDs that match the column names of \code{comm_res}. Mixed types allowed.
#' @param site_df Optional data.frame with columns \code{site}, \code{x}, \code{y} for mapping.
#'   If provided and \code{do_map=TRUE}, a site-level map is produced.
#' @param overlay_sf Optional \code{sf} object to outline on the map (e.g., study boundary).
#' @param method_comp Character method for \code{vegan::decostand}. Default \code{"hellinger"}.
#' @param distance_metric Character distance name for \code{cluster::daisy}. Default \code{"gower"}.
#'   You usually want \code{"gower"} for mixed traits.
#' @param sigma_alpha Optional numeric bandwidth for the Gaussian kernel. If \code{NULL}
#'   (default), uses the median of positive pairwise Gower distances. If that is missing or zero,
#'   falls back to \code{0.5}.
#' @param row_z Logical; if \code{TRUE} (default) return a row-z-scored version \code{C_js_z} for
#'   within-site contrasts.
#' @param quantile_probs Numeric vector of probabilities to summarize site-level crowding quantiles
#'   (e.g., \code{c(0.95)}). Use \code{NULL} to skip quantiles. Default \code{0.95}.
#' @param do_heatmap Logical; if \code{TRUE} (default) draw a clustered heatmap of \code{C_js_z}
#'   using \pkg{pheatmap} if available.
#' @param show_names Logical; show row/column names on the heatmap. Default \code{TRUE}.
#' @param cluster_rows Logical; cluster sites on the heatmap. Default \code{TRUE}.
#' @param cluster_cols Logical; cluster residents on the heatmap. Default \code{TRUE}.
#' @param heatmap_palette Function returning colors (e.g., \code{viridisLite::viridis}).
#'   Default \code{viridisLite::viridis}.
#' @param do_map Logical; if \code{TRUE}, draw a ggplot map of mean crowding per site when
#'   \code{site_df} is provided. Default \code{TRUE}.
#'
#' @details
#' The Hellinger transformation places compositions in a Euclidean-friendly space.
#' Gower distance handles mixed trait types. Distances are converted to similarities with
#' \eqn{K_{ij}=\exp(-D_{ij}^2/(2\sigma_\alpha^2))}, and the diagonal is set to 0 to remove
#' self-crowding. Crowding for resident \eqn{j} at site \eqn{s} is the weighted sum of trait-similar
#' resident composition at that site: \eqn{C_{sj} = \sum_{r} W_{sr} K_{jr}}.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{W_site} — Hellinger composition matrix (sites × residents).
#'   \item \code{D_res} — resident × resident Gower distance matrix.
#'   \item \code{sigma_alpha} — Gaussian bandwidth used.
#'   \item \code{K_res_res} — Gaussian similarity kernel with zero diagonal.
#'   \item \code{C_js} — resident crowding matrix (sites × residents).
#'   \item \code{C_js_z} — row-z-scored version of \code{C_js} (if \code{row_z=TRUE}).
#'   \item \code{site_summary} — data.frame with per-site summaries (mean and requested quantiles).
#'   \item \code{heatmap} — \code{pheatmap} object or \code{NULL} if not drawn.
#'   \item \code{map} — \code{ggplot} object or \code{NULL} if not drawn.
#' }
#'
#' @examples
#' \dontrun{
#' # Fake data: 6 sites × 5 residents
#' set.seed(42)
#' comm_res = matrix(rpois(30, lambda = 5), nrow = 6,
#'                    dimnames = list(paste0("s", 1:6), paste0("sp", 1:5)))
#' traits_res = data.frame(
#'   size = rnorm(5, 10, 2),
#'   diet = factor(sample(c("A","B"), 5, TRUE)),
#'   row.names = paste0("sp", 1:5)
#' )
#' site_df = data.frame(site = paste0("s", 1:6),
#'                       x = rep(1:3, each=2),
#'                       y = rep(1:2, 3))
#'
#' out = compute_resident_crowding(
#'   comm_res   = comm_res,
#'   traits_res = traits_res,
#'   site_df    = site_df,
#'   quantile_probs = 0.95,
#'   do_heatmap = TRUE,
#'   do_map     = TRUE
#' )
#'
#' str(out, max.level = 1)
#' out$heatmap  # pheatmap object
#' out$map      # ggplot map of mean C per site
#' }
#'
#' @importFrom stats median
#' @importFrom matrixStats rowQuantiles
#' @export
compute_resident_crowding = function(
    comm_res,
    traits_res,
    site_df = NULL,
    overlay_sf = NULL,
    method_comp = "hellinger",
    distance_metric = "gower",
    sigma_alpha = NULL,
    row_z = TRUE,
    quantile_probs = 0.95,
    do_heatmap = TRUE,
    show_names = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    heatmap_palette = viridisLite::viridis,
    do_map = TRUE
) {
  # --- 0) Defensive checks and alignment ---------------------------------------
  stopifnot(is.matrix(comm_res) || is.data.frame(comm_res))
  comm_res = as.matrix(comm_res)
  if (is.null(rownames(comm_res))) stop("comm_res must have rownames = site IDs.")
  if (is.null(colnames(comm_res))) stop("comm_res must have colnames = resident IDs.")

  stopifnot(is.data.frame(traits_res))
  if (is.null(rownames(traits_res))) stop("traits_res must have rownames = resident IDs.")

  # Match residents present in both comm_res (columns) and traits_res (rows)
  common_res = intersect(colnames(comm_res), rownames(traits_res))
  if (length(common_res) < 2) stop("Need at least 2 residents present in both comm_res and traits_res.")
  # Subset and reorder to the intersected resident set
  comm_res = comm_res[, common_res, drop = FALSE]
  traits_res = traits_res[common_res, , drop = FALSE]

  sites = rownames(comm_res)

  # --- 1) Hellinger composition W_site (sites × residents) ---------------------
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for decostand (Hellinger). Please install it.")
  }
  # decostand(method="hellinger") expects non-negative. We assume counts/abundances.
  W_site = as.matrix(vegan::decostand(comm_res, method = method_comp))

  # Sanity: preserve site order
  stopifnot(identical(rownames(W_site), sites))

  # --- 2) Gower distances among residents (handles mixed traits) ---------------
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Package 'cluster' is required for daisy (Gower). Please install it.")
  }
  # Mark ordered factors for daisy if present
  ord_cols = which(vapply(traits_res, is.ordered, logical(1)))
  type_arg = if (length(ord_cols)) list(ordered = ord_cols) else NULL

  D_res = cluster::daisy(traits_res, metric = distance_metric, type = type_arg)
  D_res = as.matrix(D_res)
  diag(D_res) = 0

  # --- 3) Gaussian kernel from distances, bandwidth sigma_alpha ----------------
  if (is.null(sigma_alpha)) {
    d_vals = D_res[upper.tri(D_res)]
    # Use the median of strictly positive distances to avoid zero collapse
    sigma_alpha = stats::median(d_vals[d_vals > 0], na.rm = TRUE)
    if (!is.finite(sigma_alpha) || sigma_alpha <= 0) sigma_alpha = 0.5
  } else {
    stopifnot(is.numeric(sigma_alpha), length(sigma_alpha) == 1L, sigma_alpha > 0)
  }

  # Convert Gower distances to Gaussian similarities; suppress self-crowding
  K_res_res = exp(-(D_res^2) / (2 * sigma_alpha^2))
  diag(K_res_res) = 0

  # --- 4) Resident crowding: C = W %*% t(K)  (sites × residents) ---------------
  C_js = W_site %*% t(K_res_res)
  rownames(C_js) = sites
  colnames(C_js) = rownames(traits_res)

  # Optional row-wise z-score for within-site contrasts
  C_js_z = NULL
  if (isTRUE(row_z)) {
    # row-wise z: subtract row mean and divide by row sd; guard non-finite
    C_js_z = t(scale(t(C_js)))
    C_js_z[!is.finite(C_js_z)] = 0
  }

  # --- 5) Site-level summaries --------------------------------------------------
  site_summary = data.frame(
    site   = rownames(C_js),
    C_mean = rowMeans(C_js, na.rm = TRUE),
    row.names = NULL
  )

  if (!is.null(quantile_probs)) {
    # Compute requested quantiles across residents per site
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
      stop("Package 'matrixStats' is required for rowQuantiles. Please install it or set quantile_probs = NULL.")
    }
    qmat = matrixStats::rowQuantiles(C_js, probs = quantile_probs, na.rm = TRUE)
    qmat = as.data.frame(qmat)
    # Name quantile columns clearly, e.g., C_q50, C_q95
    qnames = paste0("C_q", sprintf("%02d", round(100 * quantile_probs)))
    names(qmat) = qnames
    site_summary = cbind(site_summary, qmat)
  }

  # Join coordinates for mapping if site_df provided
  map_plot = NULL
  if (isTRUE(do_map) && !is.null(site_df)) {
    stopifnot(all(c("site", "x", "y") %in% names(site_df)))
    dfC = merge(site_summary, site_df[, c("site", "x", "y")], by = "site", all.x = TRUE)

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' not available; skipping map plot.")
    } else {
      p = ggplot2::ggplot(dfC, ggplot2::aes(x, y, fill = C_mean)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(name = "mean C_js", direction = -1) +
        ggplot2::labs(
          title = expression("Mean resident crowding " * C[js] * " per site"),
          x = "Longitude", y = "Latitude"
        ) +
        ggplot2::theme_minimal()

      # Optional overlay boundary
      if (!is.null(overlay_sf)) {
        if (requireNamespace("sf", quietly = TRUE)) {
          p = p + ggplot2::geom_sf(
            data = overlay_sf, inherit.aes = FALSE,
            fill = NA, color = "black", linewidth = 0.4
          )
        } else {
          warning("Package 'sf' not available; overlay_sf will be ignored.")
        }
      }

      map_plot = p
    }
  }

  # --- 6) Heatmap of row-z-scored C_js -----------------------------------------
  heatmap_obj = NULL
  if (isTRUE(do_heatmap)) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      warning("Package 'pheatmap' not available; skipping heatmap.")
    } else {
      M = if (!is.null(C_js_z)) C_js_z else C_js
      heatmap_obj = pheatmap::pheatmap(
        M,
        color        = heatmap_palette(100),
        cluster_rows = cluster_rows,
        cluster_cols = cluster_cols,
        show_rownames = show_names,
        show_colnames = show_names,
        main = expression("Resident crowding " * C[js] * ifelse(1==1, " (row-z)", "")) # label hint
      )
    }
  }

  # --- 7) Return bundle ---------------------------------------------------------
  list(
    W_site       = W_site,
    D_res        = D_res,
    sigma_alpha  = sigma_alpha,
    K_res_res    = K_res_res,
    C_js         = C_js,
    C_js_z       = C_js_z,
    site_summary = site_summary,
    heatmap      = heatmap_obj,
    map          = map_plot
  )
}
