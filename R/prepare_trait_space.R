# ======================================================================
# 1) TRAIT SPACE & RESIDENT CROWDING (wrapper) — NO site-z here
# ======================================================================

#' Prepare trait space and resident crowding (no site-z)
#'
#' @description
#' Optionally standardises trait inputs, builds the joint trait space, computes
#' centrality/hull, and computes **raw resident crowding \eqn{C_{js}}** (no per-site
#' standardisation here). Per-site standardisation (row-z) is deferred to
#' [model_residents()]. All plots are **always created** (and returned), but are
#' only **displayed** when `show_plots = TRUE`.
#'
#' @param fit An `invasimapr_fit` produced by `prepare_inputs()`. Must contain
#'   `fit$inputs$traits_res`, `fit$inputs$comm_res`, and (optionally) `fit$inputs$site_df`.
#' @param traits_inv Data frame (or matrix) of invader traits
#'   (rows = invaders; columns aligned to resident trait columns).
#' @param crowding_sigma Optional numeric bandwidth (Gaussian kernel) for crowding.
#'   If `NULL`, `compute_resident_crowding()` chooses a default/optimized value.
#' @param show_plots Logical master switch. If `TRUE`, plots produced by helper
#'   functions are displayed in the console. If `FALSE`, plots are still computed
#'   and returned in the output lists but are **not shown**.
#' @param do_standardise Logical; if `TRUE` and `standardise_model_inputs()` is
#'   available, standardise resident and invader trait inputs (environmental
#'   predictors handled later). Default `TRUE`.
#' @param row_z Logical; defer site-wise z-scoring (row standardisation) to
#'   `model_residents()`. Leave `FALSE` here to keep **raw** `C_js`. Default `FALSE`.
#'
#' @return The input `invasimapr_fit` with:
#' \itemize{
#'   \item `$traits`: list with Gower distances, ordinations (`Q_res`, `Q_inv`),
#'         hull, centroid, density, stored plots, and centrality summaries.
#'   \item `$crowding`: list with `W_site`, `D_res`, chosen `sigma_alpha`,
#'         `K_res_res`, and **raw** `C_js` (no z-scoring, no mean/sd here).
#'   \item `$meta`: residents, sites, invaders, and `n_invaders`.
#' }
#'
#' @details
#' **Display control:** To ensure plots are created but optionally hidden, when
#' `show_plots = FALSE` this wrapper temporarily opens a null graphics device
#' (`grDevices::pdf(file = NULL)`) around helper calls that would otherwise print
#' plots, then closes it. This suppresses on-screen output without altering any
#' returned plot objects.
#'
#' **Standardisation:** When `do_standardise = TRUE` *and* the function
#' `standardise_model_inputs()` exists, trait inputs are replaced by the
#' standardised versions if returned; otherwise, raw traits are used with a
#' message. Environmental variables (if any) are handled elsewhere.
#'
#' @seealso
#' `compute_trait_space()`, `compute_centrality_hull()`,
#' `compute_resident_crowding()`, `standardise_model_inputs()`, `model_residents()`
#'
#' @examples
#' \dontrun{
#' # Create plots but don't display them:
#' fit2 <- prepare_trait_space(fit, traits_inv, show_plots = FALSE)
#'
#' # Create and display plots:
#' fit3 <- prepare_trait_space(fit, traits_inv, show_plots = TRUE)
#' }
#'
#' @importFrom grDevices pdf dev.off
#' @export
prepare_trait_space = function(fit,
                               traits_inv,
                               crowding_sigma = NULL,
                               show_plots     = TRUE,
                               do_standardise = TRUE,
                               row_z          = FALSE) {

  # ---- fast input checks -------------------------------------------------------
  if (!inherits(fit, "invasimapr_fit")) stop("`fit` must be class 'invasimapr_fit'.")
  inputs <- fit$inputs
  if (is.null(inputs$traits_res)) stop("traits_res not found in fit$inputs.")
  if (is.null(inputs$comm_res))   stop("comm_res not found in fit$inputs.")

  # lightweight null-coalescing helper
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- 0) optional standardisation (traits only; env handled later) -----------
  tr_res_use <- inputs$traits_res
  tr_inv_use <- traits_inv

  if (isTRUE(do_standardise) && exists("standardise_model_inputs", mode = "function")) {
    std <- try(standardise_model_inputs(traits_res = tr_res_use, traits_inv = tr_inv_use),
               silent = TRUE)
    if (!inherits(std, "try-error")) {
      tr_res_use <- std$traits_res_glmm %||% tr_res_use
      tr_inv_use <- std$traits_inv_glmm %||% tr_inv_use
      fit$inputs_std <- (fit$inputs_std %||% list())
      fit$inputs_std$traits_res_glmm <- tr_res_use
      fit$inputs_std$traits_inv_glmm <- tr_inv_use
    } else {
      message("standardise_model_inputs() failed; proceeding with raw traits.")
    }
  }

  # ---- helper to optionally suppress display while still generating plots ------
  .with_optional_suppression <- function(expr, show) {
    if (isTRUE(show)) {
      eval.parent(substitute(expr))
    } else {
      grDevices::pdf(file = NULL)
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      eval.parent(substitute(expr))
    }
  }

  # ---- 1) shared trait space (always request plots from helper) ----------------
  ts <- .with_optional_suppression({
    compute_trait_space(
      traits_res = tr_res_use,
      traits_inv = tr_inv_use,
      do_plot    = TRUE,   # always create plot objects
      do_dend    = TRUE    # always create dendrogram plot
    )
  }, show = show_plots)

  # ---- 1b) ensure a replayable density "plot" exists even when hidden ----------
  make_ts_replotter <- function(ts) {
    dens      <- ts$density
    Q_res     <- ts$Q_res
    Q_inv     <- ts$Q_inv
    hull_res  <- ts$hull_res
    centroid  <- ts$centroid
    xlims     <- ts$xlim
    ylims     <- ts$ylim
    stopifnot(!is.null(dens))
    local({
      function(main_title = "Trait space density with convex hull and centroid",
               legend_line = "Hull = realised niche region; white square = centroid; black = residents; red = invaders",
               cex_main = 1, cex_sub = 0.72, cex_lab = 0.85, cex_axis = 0.75,
               highlight_level = 0.5) {

        graphics::filled.contour(
          dens,
          color.palette = viridisLite::viridis,
          xlim = xlims, ylim = ylims,
          plot.title = {
            graphics::title(
              main = main_title,
              sub  = legend_line,
              xlab = "trPC1", ylab = "trPC2",
              cex.main = cex_main, cex.sub = cex_sub, cex.lab = cex_lab
            )
          },
          plot.axes = {
            oldpar <- graphics::par(cex.axis = cex_axis)
            graphics::axis(1); graphics::axis(2)
            graphics::par(oldpar)

            if (nrow(Q_res)) graphics::points(Q_res$tr1, Q_res$tr2, pch = 19, cex = 0.55, col = "black")
            if (nrow(Q_inv)) graphics::points(Q_inv$tr1, Q_inv$tr2, pch = 19, cex = 0.55, col = "red")

            graphics::contour(dens$x, dens$y, dens$z,
                              add = TRUE, drawlabels = FALSE, lwd = 0.7, col = "grey60")

            if (is.finite(highlight_level) && length(highlight_level) == 1L) {
              lvl <- max(dens$z, na.rm = TRUE) * highlight_level
              graphics::contour(dens$x, dens$y, dens$z,
                                levels = lvl, add = TRUE,
                                drawlabels = FALSE, lwd = 2, col = "black")
            }

            if (!is.null(hull_res)) {
              graphics::lines(hull_res$tr1, hull_res$tr2, lwd = 2, col = "lightgrey")
            }

            graphics::points(centroid[1], centroid[2],
                             pch = 22, bg = "white", col = "black", lwd = 1.2, cex = 1.8)
          },
          key.title = graphics::title(main = "Density", cex.main = 0.8)
        )
        invisible(NULL)
      }
    })
  }

  if (is.null(ts$dens_plot) || inherits(ts$dens_plot, "recordedplot")) {
    ts$dens_plot <- make_ts_replotter(ts)
  }

  # ---- 2) centrality & hull ----------------------------------------------------
  ch <- compute_centrality_hull(Q_res = ts$Q_res, Q_inv = ts$Q_inv)

  # ---- 3) resident crowding (RAW C_js only; NO site-z here) -------------------
  cr <- .with_optional_suppression({
    compute_resident_crowding(
      comm_res    = inputs$comm_res,
      traits_res  = tr_res_use,
      site_df     = inputs$site_df,
      sigma_alpha = crowding_sigma,
      row_z       = row_z,
      do_heatmap  = TRUE,
      do_map      = TRUE
    )
  }, show = show_plots)

  # ---- 4) stash ---------------------------------------------------------------
  fit$traits <- list(
    gower      = ts$gower,
    Q_res      = ts$Q_res,
    Q_inv      = ts$Q_inv,
    hull       = ts$hull_res,
    centroid   = ts$centroid,
    density    = ts$density,
    plots_ts   = list(dens_plot = ts$dens_plot, dend_plot = ts$dend_plot),
    centrality = ch$df,
    hull_df    = ch$hull_df,
    plots_ch   = ch$plots
  )

  fit$crowding <- list(
    W_site      = cr$W_site,
    D_res       = cr$D_res,
    sigma_alpha = cr$sigma_alpha,
    K_res_res   = cr$K_res_res,
    C_js        = cr$C_js
  )

  # cache simple meta with minimal repeated lookups
  cm <- inputs$comm_res
  fit$meta$residents  <- colnames(cm)
  fit$meta$sites      <- rownames(cm)
  fit$meta$invaders   <- rownames(traits_inv)
  fit$meta$n_invaders <- nrow(traits_inv)

  fit
}
