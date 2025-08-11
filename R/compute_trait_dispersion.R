#' Trait Dispersion Pipeline: Gower → clustering → PCoA → density → metrics
#'
#' @description
#' Computes a Gower dissimilarity matrix on mixed-type traits, performs hierarchical
#' clustering, runs PCoA, visualizes density/centrality, and summarizes functional
#' dispersion metrics (FDis, FRic, Rao's Q).
#'
#' @param trait_df data.frame. One row per species (or unit); mixed types allowed.
#' @param species_col integer or character (default = 1). Column to exclude from distance.
#' @param k integer (default = 4). Number of clusters in dendrogram.
#' @param pcoa_dims integer (default = 2). PCoA axes retained (≥ 2).
#' @param abundance numeric vector or NULL. Optional species weights; normalized internally.
#' @param kde_n integer (default = 100). KDE grid resolution for density.
#' @param viridis_option character (default = "D"). Palette option for \code{viridisLite::viridis}.
#' @param show_density_plot logical (default = TRUE). Also compute a base \code{filled.contour}.
#' @param show_plots logical (default = FALSE). If TRUE, prints a portrait patchwork of all ggplots.
#' @param seed integer or NULL. Optional RNG seed.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{distance_matrix} (matrix), \code{hc} (\code{hclust}),
#'   \item \code{pcoa} (list from \code{cmdscale}), \code{scores} (data.frame with \code{centrality}),
#'   \item \code{centroid} (numeric), \code{metrics_df} (data.frame with FDis/FRic/RaoQ),
#'   \item \code{plots} (list of ggplots: \code{$dend}, \code{$density_gg}, \code{$centrality_hist}, \code{$metrics_bar}).
#' }
#'
#' @details
#' \strong{FRic} is computed as convex-hull volume/area in PCoA space via \code{geometry::convhulln}.
#' If degenerate, returns \code{NA} with a warning. \strong{Rao's Q} uses Euclidean distances in
#' reduced PCoA space: \eqn{0.5\sum_{i,j} p_i p_j d_{ij}}.
#'
#' @section Simulated example:
#' Use \code{\link{tdp_simulate_traits}} to generate a small, mixed-type trait table.
#'
#' @examples
#' # --- Built-in simulated example for reproducibility (mixed continuous + categorical traits) ---
#' tdp_simulate_traits = function(n = 30, seed = NULL) {
#'   if (!is.null(seed)) set.seed(seed)
#'   species = paste0("sp_", seq_len(n))
#'   # continuous
#'   t1 = rnorm(n, 0, 1)
#'   t2 = runif(n, -1, 1)
#'   # ordinal
#'   t3 = factor(sample(1:3, n, TRUE), ordered = TRUE)
#'   # binary
#'   t4 = factor(sample(c(0, 1), n, TRUE))
#'   # categorical
#'   t5 = factor(sample(LETTERS[1:4], n, TRUE))
#'   traits = data.frame(
#'     species = species,
#'     t_len   = t1,
#'     t_mass  = t2,
#'     t_rank  = t3,
#'     t_bin   = t4,
#'     t_cat   = t5,
#'     check.names = FALSE
#'   )
#'   abundance = rexp(n, rate = 1)
#'   list(traits = traits, abundance = abundance)
#' }
#'
#' # Generate simulated traits + abundance
#' sim_data = tdp_simulate_traits(n = 25, seed = 123)
#'
#' # Run pipeline and show combined patchwork
#' res = compute_trait_dispersion(
#'   trait_df    = sim_data$traits,
#'   species_col = "species",
#'   abundance   = sim_data$abundance,
#'   k           = 4,
#'   pcoa_dims   = 2,
#'   show_density_plot = FALSE,
#'   show_plots   = TRUE
#' )
#'
#' # Access metrics
#' res$metrics_df
#'
#' # Access a specific plot
#' res$plots$density_gg
#'
#' @importFrom cluster daisy
#' @importFrom stats as.dist hclust cmdscale dist weighted.mean
#' @importFrom ggplot2 ggplot aes geom_histogram geom_col labs theme_bw theme_classic guides
#' @importFrom MASS kde2d
#' @importFrom geometry convhulln
#' @importFrom viridisLite viridis
#' @importFrom patchwork plot_layout
#' @export
compute_trait_dispersion = function(trait_df,
                                      species_col = 1,
                                      k = 4,
                                      pcoa_dims = 2,
                                      abundance = NULL,
                                      kde_n = 100,
                                      viridis_option = "D",
                                      show_density_plot = TRUE,
                                      show_plots = FALSE,
                                      seed = NULL) {
  # ---- input checks ---------------------------------------------------------
  if (!is.data.frame(trait_df)) stop("trait_df must be a data.frame.")
  sp_col_idx = if (is.character(species_col)) {
    match(species_col, names(trait_df))
  } else {
    as.integer(species_col)
  }
  if (is.na(sp_col_idx) || sp_col_idx < 1L || sp_col_idx > ncol(trait_df))
    stop("species_col not found / out of range.")
  if (pcoa_dims < 2L) stop("pcoa_dims must be >= 2.")
  if (!is.null(seed)) set.seed(seed)

  n = nrow(trait_df)

  # ---- weights --------------------------------------------------------------
  if (is.null(abundance)) {
    p = rep(1 / n, n)
  } else {
    if (!is.numeric(abundance) || length(abundance) != n || any(abundance < 0))
      stop("abundance must be non-negative numeric, length == nrow(trait_df).")
    s = sum(abundance); if (s == 0) stop("sum(abundance) is zero.")
    p = as.numeric(abundance) / s
  }

  # ---- Gower dissimilarity --------------------------------------------------
  gower_obj  = cluster::daisy(trait_df[, -sp_col_idx, drop = FALSE], metric = "gower")
  trait_dist = as.matrix(gower_obj)

  # ---- clustering & ordination ---------------------------------------------
  hc   = stats::hclust(stats::as.dist(gower_obj))
  pcoa = stats::cmdscale(gower_obj, eig = TRUE, k = max(2L, pcoa_dims))
  scores = as.data.frame(pcoa$points)[, seq_len(pcoa_dims), drop = FALSE]
  colnames(scores) = paste0("PCoA", seq_len(pcoa_dims))

  # ---- centroid & centrality -----------------------------------------------
  centroid = vapply(seq_len(pcoa_dims), function(j) stats::weighted.mean(scores[[j]], w = p), numeric(1))
  scores$centrality = sqrt(rowSums((as.matrix(scores[, seq_len(pcoa_dims)]) -
                                       matrix(centroid, nrow(scores), pcoa_dims, byrow = TRUE))^2))

  # ---- metrics --------------------------------------------------------------
  FDis = sum(p * scores$centrality)

  FRic = NA_real_
  if (n >= (pcoa_dims + 1L)) {
    ch_try = try(geometry::convhulln(as.matrix(scores[, seq_len(pcoa_dims)]), options = "FA"), silent = TRUE)
    if (!inherits(ch_try, "try-error") && !is.null(ch_try$vol)) {
      FRic = as.numeric(ch_try$vol)
    } else {
      warning("FRic not computed: convex hull failed; returning NA.")
    }
  } else {
    warning(sprintf("FRic requires at least %d points; have %d. Returning NA.", pcoa_dims + 1L, n))
  }

  dmat = as.matrix(stats::dist(scores[, seq_len(pcoa_dims)]))
  RaoQ = 0.5 * sum(outer(p, p) * dmat)

  metrics_df = data.frame(Metric = c("FDis", "FRic", "RaoQ"),
                           Value  = c(FDis, FRic, RaoQ),
                           row.names = NULL, check.names = FALSE)

  # ---- plots (always constructed) ------------------------------------------
  k_cols = viridisLite::viridis(k, option = viridis_option)
  dend_gg = factoextra::fviz_dend(
    hc, k = k, cex = 0.5,
    k_colors = k_cols, color_labels_by_k = TRUE,
    rect = TRUE, rect_border = "grey40",
    main = "Gower Cluster Dendrogram"
  ) + ggplot2::guides(scale = "none")

  scores2 = scores[, c("PCoA1", "PCoA2")]
  xrange = range(scores2$PCoA1); xpad = 0.1 * diff(xrange); xlims = xrange + c(-xpad, xpad)
  yrange = range(scores2$PCoA2); ypad = 0.1 * diff(yrange); ylims = yrange + c(-ypad, ypad)
  kd = MASS::kde2d(scores2$PCoA1, scores2$PCoA2, n = kde_n, lims = c(xlims, ylims))

  if (isTRUE(show_density_plot)) {
    filled.contour(
      kd,
      color.palette = function(n) viridisLite::viridis(n, option = viridis_option),
      xlim = xlims, ylim = ylims,
      plot.title = title(main = "Trait Space Density Contours", xlab = "PCoA1", ylab = "PCoA2"),
      plot.axes = {
        axis(1); axis(2)
        points(scores2, pch = 19, cex = 0.5)
        contour(kd$x, kd$y, kd$z, add = TRUE, drawlabels = FALSE, lwd = 0.7, col = "grey60")
        contour(kd$x, kd$y, kd$z, add = TRUE, drawlabels = FALSE,
                levels = max(kd$z) * 0.5, lwd = 2, col = "black")
      },
      key.title = title(main = "Density")
    )
  }

  density_gg = ggplot2::ggplot(scores2, ggplot2::aes(PCoA1, PCoA2)) +
    ggplot2::stat_density_2d_filled(contour = TRUE) +
    ggplot2::geom_point(size = 0.7) +
    ggplot2::labs(title = "Trait Space Density (PCoA1–PCoA2)", x = "PCoA1", y = "PCoA2") +
    ggplot2::theme_bw()

  centrality_hist = ggplot2::ggplot(scores, ggplot2::aes(x = centrality)) +
    ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "white") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Distance to community centroid", y = "Number of species",
                  title = "Trait Centrality (Community Edge vs Core)")

  metrics_bar = ggplot2::ggplot(metrics_df, ggplot2::aes(x = Metric, y = Value)) +
    ggplot2::geom_col(width = 0.6, fill = "firebrick") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Community-Level Trait Dispersion", y = "Metric value")

  plots = list(
    dend = dend_gg,
    density_gg = density_gg,
    centrality_hist = centrality_hist,
    metrics_bar = metrics_bar
  )

  # ---- optional console display (patchwork) --------------------------------
  if (isTRUE(show_plots)) {
    combined = plots$dend /
      (plots$centrality_hist | plots$metrics_bar) /
      plots$density_gg +
      patchwork::plot_layout(heights = c(1, 2, 1))
    print(combined)
  }

  # ---- return ---------------------------------------------------------------
  list(
    distance_matrix = trait_dist,
    hc = hc,
    pcoa = pcoa,
    scores = scores,
    centroid = centroid,
    metrics_df = metrics_df,
    plots = plots
  )
}
