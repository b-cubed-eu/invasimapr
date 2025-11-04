#' Prepare trait space and resident crowding (no site-z)
#'
#' Orchestrates the *traits → space → centrality/hull → crowding* pipeline for
#' residents and invaders. Optionally standardises trait inputs, constructs a
#' joint trait space, computes convex hull / centroid / centrality, and derives
#' **raw** resident crowding \eqn{C_{js}} using a Gaussian kernel (no per-site
#' standardisation here). Row-wise z-scoring is deferred to
#' \href{https://b-cubed-eu.github.io/invasimapr/reference/model_residents.html}{`model_residents()`}.
#' Plot objects are *always* created and returned; they are only displayed when
#' `show_plots = TRUE`.
#'
#' @param fit An `invasimapr_fit` produced by
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.html}{`prepare_inputs()`/`assemble_matrices()`}.
#'   Must contain `fit$inputs$traits_res`, `fit$inputs$comm_res`, and
#'   (optionally) `fit$inputs$site_df`.
#' @param traits_inv Data frame (or matrix) of **raw** invader traits
#'   (rows = invaders; columns aligned to resident trait columns).
#' @param crowding_sigma Optional numeric bandwidth for the resident crowding
#'   kernel. If `NULL`, a default/optimised value is chosen inside
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.html}{`compute_resident_crowding()`}.
#' @param show_plots Logical. If `TRUE`, display plots as they are created.
#'   If `FALSE`, plots are still created/returned but not shown. Default `TRUE`.
#' @param do_standardise Logical. If `TRUE` and
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.html}{`standardise_model_inputs()`}
#'   exists, standardise resident *and* invader trait inputs (environment handled later).
#'   Default `TRUE`.
#' @param row_z Logical. Whether to perform **row-wise** (site) z-scoring inside
#'   `compute_resident_crowding()`. Leave `FALSE` here to keep **raw** `C_js`;
#'   row-z is typically applied in
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/model_residents.html}{`model_residents()`}.
#'
#' @details
#' **Pipeline**
#' 1. *(Optional)* Standardise `traits_res/traits_inv` via
#'    \href{https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.html}{`standardise_model_inputs()`}
#'    when available; otherwise pass through raw traits.
#' 2. Build a **joint trait space** with
#'    \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_trait_space.html}{`compute_trait_space()`}
#'    (returns Gower distances, ordination scores for residents `Q_res` and invaders
#'    `Q_inv`, a density surface, dendrograms, etc.). Plots are requested but may
#'    be hidden depending on `show_plots`.
#' 3. Compute **centrality & hull** (resident convex hull, centroid, centrality
#'    scores) using
#'    \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_centrality_hull.html}{`compute_centrality_hull()`}.
#' 4. Compute **resident crowding** with
#'    \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.html}{`compute_resident_crowding()`},
#'    returning kernel weights `K_res_res`, distances `D_res`, site kernels `W_site`,
#'    chosen `sigma_alpha`, and **raw** `C_js` (unless `row_z = TRUE` is requested).
#'
#' **Display control**
#' To ensure reproducible plot objects without cluttering the console,
#' the function temporarily opens a null graphics device when `show_plots = FALSE`.
#'
#' **Assumptions & safeguards**
#' - Requires resident community and trait tables in `fit$inputs`.
#' - If standardisation fails or is unavailable, the function proceeds with raw traits.
#' - Site-wise z-scoring is intentionally **not** applied here by default.
#'
#' @return The input `invasimapr_fit` with updated components:
#' \describe{
#'   \item{`$traits`}{List with `gower`, `Q_res`, `Q_inv`, `hull`, `centroid`,
#'         `density`, stored plots (`plots_ts`, `plots_ch`), centrality table, and
#'         hull vertices.}
#'   \item{`$crowding`}{List with `W_site`, `D_res`, `sigma_alpha`, `K_res_res`,
#'         and **raw** `C_js` (unless `row_z = TRUE`).}
#'   \item{`$meta`}{Lightweight cache of `residents`, `sites`, `invaders`,
#'         and `n_invaders`.}
#' }
#'
#' @seealso
#' - Trait space: \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_trait_space.html}{`compute_trait_space()`}
#' - Centrality & hull: \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_centrality_hull.html}{`compute_centrality_hull()`}
#' - Resident crowding: \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_resident_crowding.html}{`compute_resident_crowding()`}
#' - Standardisation: \href{https://b-cubed-eu.github.io/invasimapr/reference/standardise_model_inputs.html}{`standardise_model_inputs()`}
#' - Downstream modelling: \href{https://b-cubed-eu.github.io/invasimapr/reference/model_residents.html}{`model_residents()`}
#'
#' @examples
#' \dontrun{
#' # Compute all plots but do not display them
#' fit2 <- prepare_trait_space(fit, traits_inv, show_plots = FALSE)
#'
#' # Display plots during construction
#' fit3 <- prepare_trait_space(fit, traits_inv, show_plots = TRUE)
#'
#' # Use a fixed kernel bandwidth and request row-z inside the crowding step
#' fit4 <- prepare_trait_space(fit, traits_inv, crowding_sigma = 0.35, row_z = TRUE)
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

  # ---- Fast input checks ------------------------------------------------------
  if (!inherits(fit, "invasimapr_fit"))
    stop("`fit` must be class 'invasimapr_fit'.")

  inputs <- fit$inputs
  if (is.null(inputs$traits_res)) stop("`traits_res` not found in fit$inputs.")
  if (is.null(inputs$comm_res))   stop("`comm_res` not found in fit$inputs.")

  # Lightweight null-coalescing helper
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- 0) Optional standardisation (traits only; env handled later) ----------
  tr_res_use <- inputs$traits_res
  tr_inv_use <- traits_inv

  if (isTRUE(do_standardise) && exists("standardise_model_inputs", mode = "function")) {
    std <- try(
      standardise_model_inputs(traits_res = tr_res_use, traits_inv = tr_inv_use),
      silent = TRUE
    )
    if (!inherits(std, "try-error")) {
      tr_res_use <- std$traits_res_glmm %||% tr_res_use
      tr_inv_use <- std$traits_inv_glmm %||% tr_inv_use
      # cache for downstream steps (predictors / residents model)
      fit$inputs_std <- (fit$inputs_std %||% list())
      fit$inputs_std$traits_res_glmm <- tr_res_use
      fit$inputs_std$traits_inv_glmm <- tr_inv_use
    } else {
      message("standardise_model_inputs() failed; proceeding with raw traits.")
    }
  }

  # ---- Helper: evaluate with optional plot suppression -----------------------
  .with_optional_suppression <- function(expr, show) {
    if (isTRUE(show)) {
      eval.parent(substitute(expr))
    } else {
      grDevices::pdf(file = NULL)
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      eval.parent(substitute(expr))
    }
  }

  # ---- 1) Joint trait space (plots always created) ---------------------------
  ts <- .with_optional_suppression({
    compute_trait_space(
      traits_res = tr_res_use,
      traits_inv = tr_inv_use,
      do_plot    = TRUE,   # always create plot objects
      do_dend    = TRUE    # always create dendrogram plot
    )
  }, show = show_plots)

  # ---- 1b) Ensure a replayable density plot function exists ------------------
  # Some devices return recorded plots; provide a callable re-plotter for safety.
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

  # ---- 2) Centrality & hull ---------------------------------------------------
  ch <- compute_centrality_hull(Q_res = ts$Q_res, Q_inv = ts$Q_inv)

  # ---- 3) Resident crowding (RAW C_js here; row-z typically later) -----------
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

  # ---- 4) Stash structured outputs into `fit` --------------------------------
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

  # Cache minimal metadata for downstream steps
  cm <- inputs$comm_res
  fit$meta$residents  <- colnames(cm)
  fit$meta$sites      <- rownames(cm)
  fit$meta$invaders   <- rownames(traits_inv)
  fit$meta$n_invaders <- nrow(traits_inv)

  fit
}
