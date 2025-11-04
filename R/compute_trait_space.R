#' Compute and plot shared trait space for residents and invaders
#'
#' @title Shared trait-space construction (Gower + PCoA), resident hull, centroid, and density plot
#'
#' @description
#' Builds a unified trait space for resident and invader species by:
#' (1) computing a Gower distance on the stacked trait table,
#' (2) projecting to 2D via classical MDS/PCoA (`stats::cmdscale`),
#' (3) deriving the resident convex hull (realised niche region),
#' (4) estimating a kernel density over the 2D space, and
#' (5) optionally rendering a base R `filled.contour` map that overlays
#' the hull, all points (residents black, invaders red), and the cloud centroid (white square with black outline).
#' Optionally it also computes a functional dendrogram (`hclust`) and a
#' pretty dendrogram plot (`factoextra::fviz_dend`) when available.
#'
#' @param traits_res A data.frame of resident species (rows) by traits (columns).
#'   Mixed types are supported (numeric, integer, factor, ordered) via Gower distance.
#'   Row names are used as resident IDs; if absent, sequential IDs are created.
#' @param traits_inv A data.frame of invader species (rows) by traits (columns).
#'   Must share the same trait columns (names and types) as `traits_res`.
#'   Row names are used as invader IDs; if absent, sequential IDs are created.
#'   If you have no invaders, pass a 0-row data.frame with matching columns.
#' @param k Integer embedding dimension for PCoA; default 2.
#' @param kde_n Integer grid size for 2D kernel density estimation (per axis). Default `100`.
#' @param pad_prop Numeric padding proportion added to the plotting range on each axis. Default `0.10` (10%).
#' @param main_title Character main title for the contour plot.
#' @param legend_line Character subtitle used as an inline legend in the title area.
#' @param cex_main,cex_sub,cex_lab,cex_axis Numeric text scaling for main title, subtitle,
#'   axis labels, and axis tick labels, respectively. Defaults `1, 0.72, 0.85, 0.75`.
#' @param highlight_level Numeric in (0,1]; draws a bold contour at this proportion
#'   of the maximum density (e.g., `0.5`=50% of max). Set `NA` to skip.
#' @param do_plot Logical; if `TRUE` renders the `filled.contour` figure. Default `TRUE`.
#' @param do_dend Logical; if `TRUE`, computes `hclust` on Gower and attempts to produce a
#'   dendrogram plot with `factoextra::fviz_dend` when available. Default `FALSE`.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{gower}  — Gower distance matrix (all species).
#'   \item \code{scores} — Data frame of 2D PCoA scores (\code{Q_all}) with columns \code{tr1,tr2}.
#'   \item \code{Q_res}, \code{Q_inv} — Subsets of \code{scores} for residents and invaders.
#'   \item \code{hull_res} — Data frame (tr1,tr2) of the resident convex hull ring (closed), or \code{NULL} if < 3 residents.
#'   \item \code{centroid} — Numeric vector of the overall (res+inv) centroid \code{c(tr1, tr2)}.
#'   \item \code{density} — List \code{list(x, y, z)} from \code{MASS::kde2d}.
#'   \item \code{colors}  — Named vector with color per species ID (residents black, invaders red).
#'   \item \code{hc}      — \code{hclust} object (if \code{do_dend=TRUE}), else \code{NULL}.
#'   \item \code{dend_plot} — A ggplot object from \code{factoextra::fviz_dend} when available, else \code{NULL}.
#'   \item \code{xlim}, \code{ylim} — Numeric length-2 ranges used for plotting with padding.
#' }
#'
#' @details
#' This function assumes that the trait columns in \code{traits_inv} are compatible with those in
#' \code{traits_res}. Gower distance (\code{cluster::daisy}) supports mixed data types and handles
#' factors/ordered factors appropriately. The PCoA (`cmdscale`) is computed on the Gower distances
#' and the first two axes are returned by default for visualisation. The resident convex hull is
#' computed on the resident scores only and is drawn as the "realised niche region".
#'
#' @examples
#' \dontrun{
#' # Minimal reproducible example with numeric traits
#' set.seed(1)
#' traits_res = data.frame(
#'   trait_1 = rnorm(20),
#'   trait_2 = rnorm(20, 2),
#'   row.names = paste0("sp", 1:20)
#' )
#' traits_inv = data.frame(
#'   trait_1 = rnorm(5, 0.5),
#'   trait_2 = rnorm(5, 2.5),
#'   row.names = paste0("inv", 1:5)
#' )
#'
#' out = compute_trait_space(
#'   traits_res = traits_res,
#'   traits_inv = traits_inv,
#'   do_plot    = TRUE,
#'   do_dend    = TRUE,
#'   main_title = "Trait space density with hull and centroid",
#'   legend_line = "Hull = realised niche region; white square = centroid; black = residents; red = invaders"
#' )
#'
#' # Access returned objects
#' head(out$scores)
#' out$centroid
#' if (!is.null(out$dend_plot)) print(out$dend_plot)
#' }
#'
#' @importFrom cluster daisy
#' @importFrom MASS kde2d
#' @importFrom stats cmdscale hclust as.dist
#' @importFrom grDevices chull
#' @importFrom graphics filled.contour title axis points lines contour par
#' @export
compute_trait_space = function(traits_res,
                            traits_inv,
                            k = 2,
                            kde_n = 100,
                            pad_prop = 0.10,
                            main_title = "Trait space density with convex hull and centroid",
                            legend_line = "Hull = realised niche region; white square = centroid; black dots = residents; red dots = invaders",
                            cex_main = 1, cex_sub = 0.72, cex_lab = 0.85, cex_axis = 0.75,
                            highlight_level = 0.5,
                            do_plot = TRUE,
                            do_dend = FALSE) {

  # ---- 0) Input checks & ID hygiene -------------------------------------------
  stopifnot(is.data.frame(traits_res), is.data.frame(traits_inv))
  if (is.null(rownames(traits_res))) rownames(traits_res) = paste0("res_", seq_len(nrow(traits_res)))
  if (is.null(rownames(traits_inv))) rownames(traits_inv) = paste0("inv_", seq_len(nrow(traits_inv)))

  # Coerce factors/ordered to ensure Gower handles them properly
  # (cluster::daisy is tolerant; we still keep original types.)
  # Ensure same column set / order
  common_cols = intersect(colnames(traits_res), colnames(traits_inv))
  if (length(common_cols) == 0L && nrow(traits_inv) > 0) {
    stop("traits_res and traits_inv share no common trait columns.")
  }
  # Align both tables to the common trait set (if invaders present)
  if (nrow(traits_inv) > 0) {
    traits_res2 = traits_res[, common_cols, drop = FALSE]
    traits_inv2 = traits_inv[, common_cols, drop = FALSE]
  } else {
    traits_res2 = traits_res
    traits_inv2 = traits_inv
  }

  # ---- 1) Stack species and compute Gower distance ----------------------------
  # all_traits = rbind(traits_res2, traits_inv2)
  # Ensure categorical columns are factors (not character) and align levels
  to_factor = function(df) {
    df = as.data.frame(df, stringsAsFactors = FALSE)
    for (nm in names(df)) if (is.character(df[[nm]])) df[[nm]] = factor(df[[nm]])
    df
  }

  # FIX:
  traits_res = to_factor(traits_res)
  traits_inv = to_factor(traits_inv)

  # Harmonize factor levels across res + inv so Gower treats them consistently
  common = intersect(names(traits_res), names(traits_inv))
  for (nm in common) {
    if (is.factor(traits_res[[nm]]) || is.factor(traits_inv[[nm]])) {
      lev = union(levels(factor(traits_res[[nm]])), levels(factor(traits_inv[[nm]])))
      if (is.factor(traits_res[[nm]])) traits_res[[nm]] = factor(traits_res[[nm]], levels = lev)
      if (is.factor(traits_inv [[nm]])) traits_inv [[nm]] = factor(traits_inv [[nm]], levels = lev)
    }
  }

  # continue with your existing code, e.g.:
  all_traits = rbind(
    cbind(.role = "res", traits_res[ , common, drop = FALSE]),
    cbind(.role = "inv", traits_inv [ , common, drop = FALSE])
  )

  # Gower distance supports mixed types; converts internally as needed
  # gower_all = cluster::daisy(all_traits, metric = "gower")
  # drop the .role column before daisy if needed
  gower_all = cluster::daisy(all_traits[ , setdiff(names(all_traits), ".role"), drop = FALSE], metric = "gower")
  gower_mat = as.matrix(gower_all)

  # ---- 2) 2D PCoA (cmdscale) --------------------------------------------------
  # Classical MDS on distances; returns coordinates (scores)
  pcoa_xy = stats::cmdscale(gower_mat, k = k)
  if (k < 2) stop("k must be >= 2 to plot the trait plane.")
  Q_all = data.frame(tr1 = pcoa_xy[, 1], tr2 = pcoa_xy[, 2],
                      row.names = rownames(all_traits))

  # Split back to residents / invaders
  res_ids = rownames(traits_res2)
  inv_ids = rownames(traits_inv2)
  Q_res = Q_all[res_ids, , drop = FALSE]
  Q_inv = Q_all[inv_ids, , drop = FALSE]

  # ---- 3) Resident convex hull and centroid -----------------------------------
  hull_res = NULL
  if (nrow(Q_res) >= 3) {
    # chull returns vertex indices; close the polygon by repeating the first
    h_idx = grDevices::chull(Q_res$tr1, Q_res$tr2)
    hull_res = Q_res[c(h_idx, h_idx[1]), c("tr1", "tr2")]
  }
  # Overall centroid is computed on all species (res+inv)
  cloud_centroid = colMeans(Q_all[, c("tr1", "tr2"), drop = FALSE])

  # ---- 4) 2D Kernel density estimate in the trait plane -----------------------
  # Plot limits padded to avoid clipping aesthetics at the edges
  xr = range(Q_all$tr1, na.rm = TRUE); xw = diff(xr); xpad = ifelse(xw > 0, pad_prop * xw, 0.5)
  yr = range(Q_all$tr2, na.rm = TRUE); yw = diff(yr); ypad = ifelse(yw > 0, pad_prop * yw, 0.5)
  xlims = xr + c(-xpad, xpad)
  ylims = yr + c(-ypad, ypad)

  dens = MASS::kde2d(Q_all$tr1, Q_all$tr2, n = kde_n, lims = c(xlims, ylims))

  # ---- 5) Colors: residents (black), invaders (red) ---------------------------
  species_ids = rownames(Q_all)
  is_inv = species_ids %in% inv_ids
  col_vec = ifelse(is_inv, "red", "black")
  names(col_vec) = species_ids

  # ---- 6) Density plot (base) -------------------------------------------------
  dens_plot = NULL
  if (isTRUE(do_plot)) {
    # draw the plot
    filled.contour(
      dens,
      color.palette = viridisLite::viridis,
      xlim = xlims, ylim = ylims,

      plot.title = {
        graphics::title(
          main = main_title,
          sub  = legend_line,
          xlab = "trPC1", ylab = "trPC2",
          cex.main = cex_main,
          cex.sub  = cex_sub,
          cex.lab  = cex_lab
        )
      },

      plot.axes = {
        oldpar = graphics::par(cex.axis = cex_axis)
        graphics::axis(1); graphics::axis(2)
        graphics::par(oldpar)

        graphics::points(Q_all$tr1, Q_all$tr2, pch = 19, cex = 0.55, col = col_vec)

        graphics::contour(dens$x, dens$y, dens$z,
                          add = TRUE, drawlabels = FALSE, lwd = 0.7, col = "grey60")

        if (is.finite(highlight_level) && length(highlight_level) == 1L) {
          lvl = max(dens$z, na.rm = TRUE) * highlight_level
          graphics::contour(dens$x, dens$y, dens$z,
                            levels = lvl, add = TRUE,
                            drawlabels = FALSE, lwd = 2, col = "black")
        }

        if (!is.null(hull_res)) {
          graphics::lines(hull_res$tr1, hull_res$tr2, lwd = 2, col = "lightgrey")
        }

        graphics::points(cloud_centroid[1], cloud_centroid[2],
                         pch = 22, bg = "white", col = "black", lwd = 1.2, cex = 1.8)
      },

      key.title = graphics::title(main = "Density", cex.main = 0.8)
    )

    # =- NEW: capture the fully rendered plot as a recordedplot
    dens_plot = grDevices::recordPlot()

    # (optional) show it now as well
    # grDevices::replayPlot(dens_plot)
  }

  # ---- 7) Optional dendrogram -------------------------------------------------
  hc = dend_plot = NULL
  if (isTRUE(do_dend)) {
    hc = stats::hclust(stats::as.dist(gower_mat))
    if (requireNamespace("factoextra", quietly = TRUE) &&
        requireNamespace("viridis", quietly = TRUE)) {
      dend_plot = factoextra::fviz_dend(
        hc, k = 4, cex = 0.5,
        k_colors = viridis::viridis(4, option = "D"),
        color_labels_by_k = TRUE, rect = TRUE, rect_border = "grey40",
        main = "Gower Cluster Dendrogram"
      )
    }
    if (!is.null(dend_plot)) print(dend_plot)
  }

  # ---- 8) Return --------------------------------------------------------------
  list(
    gower     = gower_all,
    scores    = Q_all,
    Q_res     = Q_res,
    Q_inv     = Q_inv,
    hull_res  = hull_res,
    centroid  = cloud_centroid,
    density   = dens,
    colors    = col_vec,
    hc        = hc,
    dens_plot = dens_plot,   # recordedplot (replay with grDevices::replayPlot)
    dend_plot = dend_plot,   # ggplot object (print directly)
    xlim      = xlims,
    ylim      = ylims
  )
}
