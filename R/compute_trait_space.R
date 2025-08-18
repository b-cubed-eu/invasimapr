#' Compute Trait Space: per-trait similarity + community dispersion
#'
#' @description
#' A unified workflow that (1) computes within-trait similarity (numeric + categorical)
#' and (2) computes community-level dispersion in trait space
#' (Gower → clustering → PCoA → density → metrics).
#'
#' @param trait_df data.frame. One row per species (or unit); mixed types allowed.
#' @param species_col integer or character (default = 1).
#'   Column indicating the species ID; excluded from distances.
#' @param do_similarity logical (default = TRUE). If TRUE, compute per-trait similarity.
#' @param similarity_cols NULL, character, or integer. Which columns to use for similarity;
#'   default = all columns except \code{species_col}.
#' @param do_dispersion logical (default = TRUE). If TRUE, run the dispersion pipeline.
#' @param k integer (default = 4). Number of clusters for dendrogram.
#' @param pcoa_dims integer (default = 2). Number of PCoA axes retained (≥ 2).
#' @param abundance numeric vector or NULL. Optional species weights; normalized internally.
#' @param kde_n integer (default = 100). KDE grid resolution for density.
#' @param viridis_option character (default = "D"). Palette for \code{viridisLite::viridis()}.
#' @param show_density_plot logical (default = TRUE). Also emit a base \code{filled.contour}.
#' @param show_plots logical (default = FALSE). If TRUE, prints a patchwork of ggplots.
#' @param seed integer or NULL. Optional RNG seed.
#'
#' @return A list with (present elements depend on flags):
#' \itemize{
#'   \item \code{similarity}: data.frame with columns \code{Trait}, \code{Similarity} (0-100).
#'   \item \code{dispersion}: list containing
#'         \code{distance_matrix}, \code{hc}, \code{pcoa}, \code{scores}, \code{centroid},
#'         \code{metrics_df}, and \code{plots} (ggplots: \code{$dend}, \code{$density_gg},
#'         \code{$centrality_hist}, \code{$metrics_bar}).
#' }
#'
#' @examples
#' # Simulate a small mixed-type table
#' set.seed(123)
#' n <- 20
#' trait_df <- data.frame(
#'   species = paste0("sp_", seq_len(n)),
#'   height = rnorm(n),
#'   mass = runif(n, -1, 1),
#'   rank = factor(sample(1:3, n, TRUE), ordered = TRUE),
#'   bin = factor(sample(c(0, 1), n, TRUE)),
#'   cat = factor(sample(LETTERS[1:4], n, TRUE)),
#'   check.names = FALSE
#' )
#' abundance <- rexp(n)
#'
#' out <- compute_trait_space(trait_df,
#'   species_col = "species", abundance = abundance,
#'   k = 3, pcoa_dims = 2, show_density_plot = FALSE
#' )
#' out$similarity
#' out$dispersion$metrics_df
#'
#' @importFrom cluster daisy
#' @importFrom stats as.dist hclust cmdscale dist weighted.mean
#' @importFrom ggplot2 ggplot aes geom_histogram geom_col labs theme_bw theme_classic guides
#' @importFrom MASS kde2d
#' @importFrom geometry convhulln
#' @importFrom viridisLite viridis
#' @importFrom patchwork plot_layout
#' @export
compute_trait_space <- function(trait_df,
                                species_col = 1,
                                do_similarity = TRUE,
                                similarity_cols = NULL,
                                do_dispersion = TRUE,
                                k = 4,
                                pcoa_dims = 2,
                                abundance = NULL,
                                kde_n = 100,
                                viridis_option = "D",
                                show_density_plot = TRUE,
                                show_plots = FALSE,
                                seed = NULL) {
  # ---- checks ----
  if (!is.data.frame(trait_df)) stop("trait_df must be a data.frame.")
  sp_col_idx <- if (is.character(species_col)) match(species_col, names(trait_df)) else as.integer(species_col)
  if (is.na(sp_col_idx) || sp_col_idx < 1L || sp_col_idx > ncol(trait_df)) {
    stop("species_col not found / out of range.")
  }
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(trait_df)
  res <- list()

  # ======================= PART 1: similarity ===============================
  if (isTRUE(do_similarity)) {
    # pick columns for similarity
    if (is.null(similarity_cols)) {
      sim_idx <- setdiff(seq_along(trait_df), sp_col_idx)
    } else {
      sim_idx <- if (is.character(similarity_cols)) {
        match(similarity_cols, names(trait_df))
      } else {
        as.integer(similarity_cols)
      }
      if (anyNA(sim_idx)) stop("similarity_cols: some names/indices not found.")
    }

    # helpers
    sim_numeric <- function(x) {
      x2 <- x[!is.na(x)]
      if (length(x2) <= 1) {
        return(1)
      }
      r <- diff(range(x2))
      if (r == 0) {
        return(1)
      }
      x_scaled <- (x2 - min(x2)) / r
      1 - mean(stats::dist(x_scaled))
    }
    sim_categ <- function(x) {
      x <- as.character(x)
      freqs <- table(x)
      if (sum(freqs) < 2) {
        return(1)
      }
      as.numeric(sum(choose(freqs, 2)) / choose(sum(freqs), 2))
    }

    trait_names <- names(trait_df)[sim_idx]
    sim_vals <- vapply(sim_idx, function(j) {
      col <- trait_df[[j]]
      if (is.numeric(col)) sim_numeric(col) else sim_categ(col)
    }, numeric(1))

    res$similarity <- data.frame(
      Trait = trait_names,
      Similarity = 100 * sim_vals,
      row.names = NULL, check.names = FALSE
    )
  }

  # ======================= PART 2: dispersion ===============================
  if (isTRUE(do_dispersion)) {
    if (pcoa_dims < 2L) stop("pcoa_dims must be >= 2.")
    # weights
    p <- if (is.null(abundance)) {
      rep(1 / n, n)
    } else {
      if (!is.numeric(abundance) || length(abundance) != n || any(abundance < 0)) {
        stop("abundance must be non-negative numeric, length == nrow(trait_df).")
      }
      s <- sum(abundance)
      if (s == 0) stop("sum(abundance) is zero.")
      as.numeric(abundance) / s
    }

    # Gower
    gower_obj <- cluster::daisy(trait_df[, -sp_col_idx, drop = FALSE], metric = "gower")
    trait_dist <- as.matrix(gower_obj)

    # clustering + PCoA
    hc <- stats::hclust(stats::as.dist(gower_obj))
    pcoa <- stats::cmdscale(gower_obj, eig = TRUE, k = max(2L, pcoa_dims))
    scores <- as.data.frame(pcoa$points)[, seq_len(pcoa_dims), drop = FALSE]
    colnames(scores) <- paste0("PCoA", seq_len(pcoa_dims))

    # centroid + centrality
    centroid <- vapply(seq_len(pcoa_dims), function(j) stats::weighted.mean(scores[[j]], w = p), numeric(1))
    scores$centrality <- sqrt(rowSums((as.matrix(scores[, seq_len(pcoa_dims)]) -
      matrix(centroid, nrow(scores), pcoa_dims, byrow = TRUE))^2))

    # metrics
    FDis <- sum(p * scores$centrality)

    FRic <- NA_real_
    if (n >= (pcoa_dims + 1L)) {
      ch_try <- try(geometry::convhulln(as.matrix(scores[, seq_len(pcoa_dims)]), options = "FA"), silent = TRUE)
      if (!inherits(ch_try, "try-error") && !is.null(ch_try$vol)) {
        FRic <- as.numeric(ch_try$vol)
      } else {
        warning("FRic not computed: convex hull failed; returning NA.")
      }
    } else {
      warning(sprintf("FRic requires at least %d points; have %d. Returning NA.", pcoa_dims + 1L, n))
    }

    dmat <- as.matrix(stats::dist(scores[, seq_len(pcoa_dims)]))
    RaoQ <- 0.5 * sum(outer(p, p) * dmat)

    metrics_df <- data.frame(
      Metric = c("FDis", "FRic", "RaoQ"),
      Value = c(FDis, FRic, RaoQ),
      row.names = NULL, check.names = FALSE
    )

    # plots
    k_cols <- viridisLite::viridis(k, option = viridis_option)
    dend_gg <- factoextra::fviz_dend(
      hc,
      k = k, cex = 0.5,
      k_colors = k_cols, color_labels_by_k = TRUE,
      rect = TRUE, rect_border = "grey40",
      main = "Gower Cluster Dendrogram"
    ) + ggplot2::guides(scale = "none")

    scores2 <- scores[, c("PCoA1", "PCoA2")]
    xrange <- range(scores2$PCoA1)
    xpad <- 0.1 * diff(xrange)
    xlims <- xrange + c(-xpad, xpad)
    yrange <- range(scores2$PCoA2)
    ypad <- 0.1 * diff(yrange)
    ylims <- yrange + c(-ypad, ypad)
    kd <- MASS::kde2d(scores2$PCoA1, scores2$PCoA2, n = kde_n, lims = c(xlims, ylims))

    if (isTRUE(show_density_plot)) {
      filled.contour(
        kd,
        color.palette = function(n) viridisLite::viridis(n, option = viridis_option),
        xlim = xlims, ylim = ylims,
        plot.title = title(main = "Trait Space Density Contours", xlab = "PCoA1", ylab = "PCoA2"),
        plot.axes = {
          axis(1)
          axis(2)
          points(scores2, pch = 19, cex = 0.5)
          contour(kd$x, kd$y, kd$z, add = TRUE, drawlabels = FALSE, lwd = 0.7, col = "grey60")
          contour(kd$x, kd$y, kd$z,
            add = TRUE, drawlabels = FALSE,
            levels = max(kd$z) * 0.5, lwd = 2, col = "black"
          )
        },
        key.title = title(main = "Density")
      )
    }

    density_gg <- ggplot2::ggplot(scores2, ggplot2::aes(PCoA1, PCoA2)) +
      ggplot2::stat_density_2d_filled(contour = TRUE) +
      ggplot2::geom_point(size = 0.7) +
      ggplot2::labs(title = "Trait Space Density (PCoA1-PCoA2)", x = "PCoA1", y = "PCoA2") +
      ggplot2::theme_bw()

    centrality_hist <- ggplot2::ggplot(scores, ggplot2::aes(x = centrality)) +
      ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "white") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        x = "Distance to community centroid", y = "Number of species",
        title = "Trait Centrality (Community Edge vs Core)") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

    metrics_bar <- ggplot2::ggplot(metrics_df, ggplot2::aes(x = Metric, y = Value)) +
      ggplot2::geom_col(width = 0.6, fill = "firebrick") +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Community-Level Trait Dispersion", y = "Metric value") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

    plots <- list(
      dend = dend_gg,
      density_gg = density_gg,
      centrality_hist = centrality_hist,
      metrics_bar = metrics_bar
    )

    if (isTRUE(show_plots)) {
      combined <- plots$dend /
        (plots$centrality_hist | plots$metrics_bar) /
        plots$density_gg +
        patchwork::plot_layout(heights = c(1, 2, 1))
      print(combined)
    }

    res$dispersion <- list(
      distance_matrix = trait_dist,
      hc = hc,
      pcoa = pcoa,
      scores = scores,
      centroid = centroid,
      metrics_df = metrics_df,
      plots = plots
    )
  }

  res
}
