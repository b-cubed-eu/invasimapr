#' Centrality and convex-hull membership in trait space
#'
#' @title Compute trait-space centrality (robust Mahalanobis) and hull status, with visuals
#'
#' @description
#' Given 2D trait coordinates for **residents** and **invaders** (typically PCoA scores),
#' this function:
#' \itemize{
#'   \item estimates a robust centre and covariance for residents (MCD via \code{MASS::cov.rob}, fallback to classical),
#'   \item computes Mahalanobis distances for residents and invaders,
#'   \item converts distances to **centrality** using the resident distance CDF
#'         (\code{centrality = 1 - F(d)}; 1 = core, 0 = peripheral),
#'   \item determines whether invaders sit **inside** the resident convex hull (familiar) or **outside** (novel),
#'   \item returns a tidy table and three ggplot figures.
#' }
#'
#' @param Q_res Data frame with resident coordinates; must contain columns \code{tr1}, \code{tr2}.
#'              Row names are treated as resident IDs.
#' @param Q_inv Data frame with invader coordinates; must contain \code{tr1}, \code{tr2}.
#'              Row names are treated as invader IDs. Can be 0-row.
#' @param ellipse_level Numeric in (0,1); confidence level for the resident normal ellipse
#'                      in the trait-plane plot. Default \code{0.50}.
#' @param point_size Numeric; point size in the trait-plane plot. Default \code{2.8}.
#' @param alpha Numeric in (0,1]; point alpha. Default \code{0.95}.
#' @param stroke Numeric; outline stroke width for points. Default \code{0.7}.
#' @param rank_by Character; one of \code{"centrality"} or \code{"distance"}.
#'               Controls the invader ranking plot. Default \code{"centrality"}.
#' @param peripheral_first Logical; if \code{TRUE} and \code{rank_by = "centrality"},
#'                         sort by increasing centrality (peripheral first). If ranking by
#'                         distance, \code{TRUE} sorts by decreasing distance (far first).
#'                         Default \code{TRUE}.
#' @param palette Centrality fill palette; passed to \code{scale_fill_viridis_c}. Default \code{"viridis"}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{df} — tidy table with columns:
#'         \code{id}, \code{grp} (\code{"resident"|"invader"}), \code{tr1}, \code{tr2},
#'         \code{d_md} (Mahalanobis), \code{d_eu} (Euclidean), \code{centrality} (0–1),
#'         \code{in_hull} (logical for invaders; residents = TRUE).
#'   \item \code{center}, \code{cov} — robust centre and covariance used.
#'   \item \code{hull_df} — closed ring (tr1,tr2) of resident convex hull, or \code{NULL} if <3 residents.
#'   \item \code{plots} — list with:
#'         \code{p_trait} (trait-plane scatter),
#'         \code{p_dist} (distance distributions),
#'         \code{p_rank} (invader ranking with hull flags; \code{NULL} if no invaders).
#' }
#'
#' @details
#' Mahalanobis distances may fail if the covariance is nearly singular; a tiny ridge is added
#' internally if needed. Hull membership is computed with \pkg{sp} when available,
#' otherwise with \pkg{sf} if available; if neither is installed and residents < 3,
#' hull membership is returned as \code{NA}.
#'
#' @examples
#' \dontrun{
#' # Assume you already computed PCoA coordinates:
#' # Q_res, Q_inv each with columns tr1, tr2 and rownames as IDs.
#' out = compute_centrality_hull(Q_res, Q_inv)
#' head(out$df)
#' out$plots$p_trait
#' out$plots$p_dist
#' out$plots$p_rank
#' }
#'
#' @importFrom MASS cov.rob
#' @importFrom stats cov mahalanobis ecdf
#' @importFrom grDevices chull
#' @importFrom ggplot2 ggplot aes geom_point geom_polygon stat_ellipse scale_fill_viridis_c
#' @importFrom ggplot2 scale_color_manual labs theme_minimal theme element_text geom_histogram
#' @importFrom ggplot2 geom_density coord_flip geom_col scale_fill_manual guides guide_legend
#' @importFrom dplyr bind_rows mutate arrange desc tibble left_join
#' @export
compute_centrality_hull = function(Q_res,
                                    Q_inv,
                                    ellipse_level = 0.50,
                                    point_size    = 2.8,
                                    alpha         = 0.95,
                                    stroke        = 0.7,
                                    rank_by       = c("centrality","distance"),
                                    peripheral_first = TRUE,
                                    palette       = "viridis") {

  # ---------------------------
  # 0) Input checks & shaping
  # ---------------------------
  stopifnot(is.data.frame(Q_res), all(c("tr1","tr2") %in% colnames(Q_res)))
  stopifnot(is.data.frame(Q_inv), all(c("tr1","tr2") %in% colnames(Q_inv)))

  R = as.data.frame(Q_res[, c("tr1","tr2"), drop = FALSE])
  I = as.data.frame(Q_inv[, c("tr1","tr2"), drop = FALSE])

  rid = rownames(R); if (is.null(rid)) rid = paste0("res_", seq_len(nrow(R)))
  iid = rownames(I); if (is.null(iid)) iid = paste0("inv_", seq_len(nrow(I)))

  # ---------------------------
  # 1) Robust centre & covariance (MCD -> classical fallback)
  # ---------------------------
  rob = try(MASS::cov.rob(R, method = "mcd"), silent = TRUE)
  if (inherits(rob, "try-error")) {
    rob = list(center = colMeans(R), cov = stats::cov(R))
  }
  center = as.numeric(rob$center); names(center) = c("tr1","tr2")
  Sigma  = as.matrix(rob$cov)

  # Safe Mahalanobis with tiny ridge if needed
  .mahal_safe = function(X, center, cov) {
    out = try(stats::mahalanobis(X, center = center, cov = cov), silent = TRUE)
    if (inherits(out, "try-error")) {
      # add a tiny ridge to stabilise
      p = ncol(cov)
      out = stats::mahalanobis(X, center = center, cov = cov + diag(1e-8, p))
    }
    as.numeric(out)
  }

  # Distances (Mahalanobis + Euclidean)
  dR_md = .mahal_safe(R, center, Sigma)
  dI_md = if (nrow(I)) .mahal_safe(I, center, Sigma) else numeric(0)

  dR_eu = sqrt(rowSums((R - matrix(center, nrow(R), 2, byrow = TRUE))^2))
  dI_eu = if (nrow(I)) sqrt(rowSums((I - matrix(center, nrow(I), 2, byrow = TRUE))^2)) else numeric(0)

  # ---------------------------
  # 2) Centrality = 1 - F_resident(d_md)
  # ---------------------------
  F_R   = stats::ecdf(dR_md)
  centR = 1 - F_R(dR_md)
  centI = if (nrow(I)) 1 - F_R(dI_md) else numeric(0)

  # ---------------------------
  # 3) Resident convex hull & invader in-hull test
  # ---------------------------
  hull_df = NULL
  in_hull_I = if (nrow(I)) rep(NA, nrow(I)) else logical(0)

  if (nrow(R) >= 3) {
    H_idx = grDevices::chull(R$tr1, R$tr2)
    H_xy  = R[c(H_idx, H_idx[1]), , drop = FALSE]  # closed ring
    hull_df = data.frame(tr1 = H_xy$tr1, tr2 = H_xy$tr2)

    if (nrow(I)) {
      if (requireNamespace("sp", quietly = TRUE)) {
        in_hull_I = sp::point.in.polygon(I$tr1, I$tr2, H_xy$tr1, H_xy$tr2) > 0
      } else if (requireNamespace("sf", quietly = TRUE)) {
        poly = sf::st_polygon(list(as.matrix(H_xy)))
        in_hull_I = as.logical(
          sf::st_within(
            sf::st_as_sf(I, coords = c("tr1","tr2"), crs = NA),
            sf::st_sfc(poly),
            sparse = FALSE
          )
        )
      } else {
        in_hull_I = rep(NA, nrow(I))  # no geometry backend available
      }
    }
  } else if (nrow(I)) {
    in_hull_I = rep(NA, nrow(I)) # cannot define a hull with < 3 points
  }

  # ---------------------------
  # 4) Tidy table (residents + invaders)
  # ---------------------------
  df_res = dplyr::tibble(
    id = rid, grp = "resident",
    tr1 = R$tr1, tr2 = R$tr2,
    d_md = dR_md, d_eu = dR_eu,
    centrality = centR,
    in_hull = TRUE
  )
  df_inv = if (nrow(I)) dplyr::tibble(
    id = iid, grp = "invader",
    tr1 = I$tr1, tr2 = I$tr2,
    d_md = dI_md, d_eu = dI_eu,
    centrality = centI,
    in_hull = in_hull_I
  ) else dplyr::tibble(
    id = character(), grp = character(),
    tr1 = numeric(), tr2 = numeric(),
    d_md = numeric(), d_eu = numeric(),
    centrality = numeric(), in_hull = logical()
  )

  df = dplyr::bind_rows(df_res, df_inv)

  # ---------------------------
  # 5) Plots
  # ---------------------------
  # (a) Trait-plane scatter (hull + 50% ellipse + centrality)
  hull_layer = if (!is.null(hull_df)) {
    ggplot2::geom_polygon(
      data = hull_df, ggplot2::aes(tr1, tr2),
      inherit.aes = FALSE, fill = NA, color = "grey30", linewidth = 0.7
    )
  } else NULL

  p_trait = ggplot2::ggplot(df, ggplot2::aes(tr1, tr2)) +
    (hull_layer %||% ggplot2::geom_blank()) +
    ggplot2::stat_ellipse(
      data = subset(df, grp == "resident"),
      type = "norm", level = ellipse_level,
      linetype = 2, linewidth = 0.6, color = "grey40"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = grp,
        fill = centrality,
        color = ifelse(grp == "invader" & !is.na(in_hull) & !in_hull, "outside", "inside")
      ),
      size = point_size, stroke = stroke, alpha = alpha
    ) +
    ggplot2::scale_shape_manual(values = c(resident = 21, invader = 24)) +
    ggplot2::scale_fill_viridis_c(direction = -1, option = palette,
                                  name = "Centrality\n(1 = core)") +
    ggplot2::scale_color_manual(values = c(inside = "black", outside = "#cc2b2b"),
                                guide = "none") +
    ggplot2::labs(
      title = "Resident vs invader centrality and hull status",
      x = "Trait axis 1", y = "Trait axis 2"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # (b) Distance distributions (Mahalanobis)
  df_dist = rbind(
    data.frame(grp = "resident", d_md = dR_md),
    if (nrow(I)) data.frame(grp = "invader",  d_md = dI_md) else NULL
  )

  p_dist = ggplot2::ggplot(df_dist, ggplot2::aes(d_md, fill = grp)) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::labs(
      title = "Mahalanobis distance distributions",
      x = "Mahalanobis distance (to resident centre)", y = "Density"
    ) +
    ggplot2::scale_fill_manual(values = c(resident = "grey40", invader = "#cc2b2b")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # (c) Invader ranking (centrality or distance), with hull flags
  p_rank = NULL
  if (nrow(I)) {
    inv = subset(df, grp == "invader")

    rb = match.arg(rank_by)
    if (rb == "centrality") {
      ord = if (isTRUE(peripheral_first)) order(inv$centrality, inv$d_md, decreasing = FALSE)
      else                          order(inv$centrality, inv$d_md, decreasing = TRUE)
      inv = inv[ord, ]
      inv$rank_val = inv$centrality
      xlab = "Centrality (1 = core)"
    } else {
      ord = if (isTRUE(peripheral_first)) order(inv$d_md, inv$centrality, decreasing = TRUE)
      else                          order(inv$d_md, inv$centrality, decreasing = FALSE)
      inv = inv[ord, ]
      inv$rank_val = inv$d_md
      xlab = "Mahalanobis distance (higher = more novel)"
    }

    p_rank = ggplot2::ggplot(inv,
                              ggplot2::aes(y = stats::reorder(id, rank_val), x = rank_val, fill = in_hull)
    ) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(
        values = c(`TRUE` = "grey60", `FALSE` = "#cc2b2b"),
        name = "In resident hull"
      ) +
      ggplot2::labs(
        title = sprintf("Invader ranking by %s %s",
                        rb, if (peripheral_first) "(peripheral first)" else "(core first)"),
        x = xlab, y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }

  # ---------------------------
  # 6) Return
  # ---------------------------
  list(
    df      = df,
    center  = center,
    cov     = Sigma,
    hull_df = hull_df,
    plots   = list(p_trait = p_trait, p_dist = p_dist, p_rank = p_rank)
  )
}

# tiny helper for null-coalescing in plot construction
`%||%` = function(a, b) if (!is.null(a)) a else b
