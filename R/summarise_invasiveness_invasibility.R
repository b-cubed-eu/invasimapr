#' Summarize invasion fitness into species, trait, and site metrics (with plots)
#'
#' @title Summaries of invasion fitness: species invasiveness, trait effects, and site invasibility
#'
#' @description
#' Invasion fitness \eqn{\lambda_{is}} integrates structure in **trait space**
#' (distances, overlaps, **convex hull**, **cloud centroid**) with **abiotic suitability**
#' (alignment to environment), **niche crowding** (overlap with residents weighted by composition),
#' and **resident competition** (site saturation). `summarise_invasiveness_invasibility()` collapses the
#' site × invader surface into **species**, **trait**, and **site** summaries that are actionable
#' and map-readyeither using a **probabilistic** measure \(P(F>0)\) or a **hard rule**
#' \eqn{I\{\lambda>0\}}.
#'
#' \strong{Species invasiveness} (per invader) summarises breadth of establishment across sites:
#' \deqn{V_i = |S|^{-1}\sum_s \mathbb{I}\{\lambda_{is}>0\}\quad \text{or}\quad
#' \tilde V_i = |S|^{-1}\sum_s P(F>0\mid i,s).}
#'
#' \strong{Site invasibility} (per site) quantifies openness to newcomers:
#' \deqn{V_s = |I|^{-1}\sum_i \mathbb{I}\{\lambda_{is}>0\}\quad \text{or}\quad
#' \tilde V_s = |I|^{-1}\sum_i P(F>0\mid i,s).}
#'
#' \strong{Trait invasiveness} scores which invader traits explain variation in \(V_i\),
#' via standardized slopes for continuous traits and ANOVA \(R^2\) for categorical traits.
#'
#' @param lambda_is Matrix \eqn{S\times I} of invasion fitness (rows = sites, cols = invaders).
#'   If `NULL`, you must supply `p_is`. Provide row/column names.
#' @param p_is Optional matrix \eqn{S\times I} of establishment probabilities
#'   (e.g., from `compute_establishment_probability()`).
#' @param use_probabilistic Logical. If `TRUE`, summaries use `p_is` (expected values).
#'   If `FALSE`, summaries use the hard rule \eqn{I\{\lambda>0\}}.
#' @param prob_threshold Numeric in (0,1). If `use_probabilistic=TRUE` and you still want a
#'   binary view for selected outputs/maps, cells with `p_is >= prob_threshold` count as 1.
#' @param site_df Optional data frame with columns `site, x, y` for mapping.
#'   If missing, maps are omitted but tabular summaries are still returned.
#' @param traits_inv Optional data frame of invader traits for trait summaries;
#'   rownames must be invader IDs (matching `colnames(lambda_is)` or `colnames(p_is)`).
#'   May contain numeric and factor columns; non‐trait columns can be present.
#' @param comm_res Optional site × resident matrix for relative metrics; used to compute
#'   per-site resident richness for extra normalization (optional).
#' @param return_long Logical. If `TRUE`, include a tidy long table for the site×invader surface.
#' @param make_plots Logical. If `TRUE`, return ggplot objects for maps/bars/heatmaps.
#' @param label Optional character label used in plot titles and output metadata.
#'
#' @return A list with:
#' \itemize{
#'   \item `species`: data frame with species-level invasiveness metrics (probabilistic and/or hard).
#'   \item `site`: data frame with site-level invasibility metrics (probabilistic and/or hard + map coords if provided).
#'   \item `trait_effects`: data frame of trait effect sizes (|β| for continuous; ANOVA \(R^2\) for categorical).
#'   \item `establish_long`: tidy long table of the working surface (`val` = probability or 0/1; includes `lambda` if available).
#'   \item `plots`: list of ggplot objects (may be `NULL` if `make_plots=FALSE` or `ggplot2` not installed), including:
#'         `site_map`, `invader_rank`, `heatmap`, `trait_effects`.
#'   \item `meta`: list with `mode` ("probabilistic" or "hard"), `threshold`, and `label`.
#' }
#'
#' @details
#' \strong{Working surface.} If `use_probabilistic=TRUE`, the core surface is `p_is`
#' (expected establishments). If `FALSE`, it is the binary matrix \eqn{I\{\lambda>0\}}.
#' If both `p_is` and `lambda_is` are supplied, both are used: the probabilistic summaries
#' use `p_is`, the hard-rule summaries come from `lambda_is`.
#'
#' \strong{Trait effects.} For each numeric trait \eqn{T_k}, we fit a simple regression
#' \deqn{\tilde V_i \sim T_{ik}}{Vtilde_i ~ T_{ik}} and report the standardized slope \eqn{|\beta_k|}{|beta_k|}.
#' For each factor trait, we report ANOVA \(R^2\) from a one-way model.
#' These are quick effect sizes to rank traits; more elaborate models can be layered later.
#'
#' @examples
#' set.seed(42)
#' S = 8; I = 5
#' sites = paste0("s", 1:S)
#' inv   = paste0("i", 1:I)
#'
#' # Fake fitness and probabilities
#' lambda_is = matrix(rnorm(S*I, sd=1), S, I, dimnames=list(sites, inv))
#' p_is      = pnorm(lambda_is)  # crude probit just for the example
#'
#' # Minimal site coordinates for plotting (optional)
#' site_df = data.frame(site = sites,
#'                       x = rep(1:4, each=2)[1:S],
#'                       y = rep(1:2, times=4)[1:S])
#'
#' # Minimal trait table for invaders (one numeric, one factor)
#' traits_inv = data.frame(trait_size = runif(I, 0, 1),
#'                          trait_type = factor(sample(c("A","B"), I, TRUE)),
#'                          row.names = inv)
#'
#' # 1) Probabilistic summaries (use p_is)
#' outP = summarise_invasiveness_invasibility(
#'   lambda_is = lambda_is,
#'   p_is      = p_is,
#'   use_probabilistic = TRUE,
#'   site_df   = site_df,
#'   traits_inv = traits_inv,
#'   make_plots = FALSE
#' )
#' names(outP)
#' head(outP$species)
#'
#' # 2) Hard rule summaries (use I{lambda>0})
#' outH = summarise_invasiveness_invasibility(
#'   lambda_is = lambda_is,
#'   use_probabilistic = FALSE,
#'   site_df = site_df,
#'   traits_inv = traits_inv,
#'   make_plots = FALSE
#' )
#' head(outH$site)
#'
#' @export
summarise_invasiveness_invasibility = function(
    lambda_is = NULL,
    p_is = NULL,
    use_probabilistic = FALSE,
    prob_threshold = 0.5,
    site_df = NULL,
    traits_inv = NULL,
    comm_res = NULL,
    return_long = TRUE,
    make_plots = TRUE,
    label = NULL
){
  if (is.null(lambda_is) && is.null(p_is)) {
    stop("Provide at least one of `lambda_is` or `p_is`.")
  }
  if (!is.null(lambda_is)) {
    stopifnot(is.matrix(lambda_is), !is.null(rownames(lambda_is)), !is.null(colnames(lambda_is)))
  }
  if (!is.null(p_is)) {
    stopifnot(is.matrix(p_is), !is.null(rownames(p_is)), !is.null(colnames(p_is)))
  }
  if (!is.null(lambda_is) && !is.null(p_is)) {
    if (!identical(dimnames(lambda_is), dimnames(p_is))) {
      stop("When both `lambda_is` and `p_is` are provided, they must have identical dimnames.")
    }
  }

  mode = if (isTRUE(use_probabilistic)) "probabilistic" else "hard"

  rn = if (!is.null(lambda_is)) rownames(lambda_is) else rownames(p_is)
  cn = if (!is.null(lambda_is)) colnames(lambda_is) else colnames(p_is)
  S  = length(rn); I = length(cn)

  W_hard = NULL
  if (!is.null(lambda_is)) W_hard = (lambda_is > 0) * 1L

  W =
    if (mode == "probabilistic") {
      if (is.null(p_is)) stop("`use_probabilistic=TRUE` requires `p_is`.")
      p_is
    } else {
      if (is.null(W_hard)) stop("`use_probabilistic=FALSE` requires `lambda_is`.")
      W_hard
    }

  # Long table
  establish_long = NULL
  if (isTRUE(return_long)) {
    establish_long = tibble::tibble(
      site    = rep(rn, times = I),
      invader = rep(cn, each  = S),
      val     = as.vector(W),
      lambda  = if (!is.null(lambda_is)) as.vector(lambda_is) else NA_real_,
      mode    = mode,
      label   = label %||% mode
    )
    if (!is.null(site_df) && all(c("site","x","y") %in% names(site_df))) {
      establish_long = dplyr::left_join(establish_long, site_df, by = "site")
    }
  }

  # Species / site summaries
  species = tibble::tibble(
    invader   = cn,
    V_i       = colMeans(W, na.rm = TRUE),
    n_sites   = S
  )
  if (!is.null(W_hard)) {
    species$V_i_hard = colMeans(W_hard, na.rm = TRUE)
    species$n_hard   = colSums(W_hard, na.rm = TRUE)
  }

  site = tibble::tibble(
    site     = rn,
    V_s      = rowMeans(W, na.rm = TRUE),
    n_inv    = I,
    total_expected = rowSums(W, na.rm = TRUE)
  )
  if (!is.null(site_df) && all(c("site","x","y") %in% names(site_df))) {
    site = dplyr::left_join(site, site_df, by = "site")
  }
  if (!is.null(W_hard)) {
    site$V_s_hard = rowMeans(W_hard, na.rm = TRUE)
    site$n_est    = rowSums(W_hard, na.rm = TRUE)
    site$n_fail   = rowSums(1 - W_hard, na.rm = TRUE)
  }

  # Optional: resident richness per site
  if (!is.null(comm_res) && is.matrix(comm_res)) {
    n_res_site = rowSums(comm_res > 0, na.rm = TRUE)
    site$n_res_site = n_res_site[match(site$site, names(n_res_site))]
  }

  # ---- Trait effects (optional; requires traits_inv) --------------------------
  trait_effects = NULL
  if (!is.null(traits_inv)) {
    dat = as.data.frame(traits_inv)
    if (is.null(rownames(dat))) stop("`traits_inv` must have rownames = invader IDs.")
    dat = tibble::as_tibble(dat, rownames = "invader")

    # Remove PCA columns so we only use interpretable traits
    pca_cols = grep("^(ENV_|TR_)PC\\d+$", names(dat), value = TRUE)
    if (length(pca_cols)) dat = dat[, setdiff(names(dat), pca_cols), drop = FALSE]

    # Merge species metric (prob or hard) with traits
    inv_metric = tibble::tibble(invader = colnames(W), V_i = colMeans(W, na.rm = TRUE))
    dat = dplyr::left_join(dat, inv_metric, by = "invader")

    # If no usable trait cols left, skip gracefully
    cols = setdiff(colnames(dat), c("invader","V_i"))
    if (length(cols)) {
      is_num = vapply(dat[, cols, drop = FALSE], is.numeric, logical(1))
      cont_cols = cols[is_num]
      cat_cols  = cols[!is_num]

      eff_cont = if (length(cont_cols)) {
        purrr::map_dfr(cont_cols, function(tr){
          df = dat[, c("V_i", tr)]
          df = tidyr::drop_na(df)
          if (nrow(df) < 3 || stats::sd(df[[tr]], na.rm = TRUE) == 0) return(NULL)
          x = scale(df[[tr]])[,1]; y = scale(df$V_i)[,1]
          fit = stats::lm(y ~ x)
          tibble::tibble(
            trait = tr, type = "continuous",
            effect = unname(stats::coef(fit)[2]),
            p = summary(fit)$coefficients[2,4]
          )
        })
      } else tibble::tibble(trait=character(), type=character(), effect=numeric(), p=numeric())

      eff_cat = if (length(cat_cols)) {
        purrr::map_dfr(cat_cols, function(tr){
          df = dat[, c("V_i", tr)]
          df = tidyr::drop_na(df)
          if (nrow(df) < 3 || dplyr::n_distinct(df[[tr]]) < 2) return(NULL)
          df[[tr]] = as.factor(df[[tr]])
          fit = stats::lm(V_i ~ ., data = df)
          aov_tab = stats::anova(fit)
          ss_between = aov_tab$`Sum Sq`[1]
          ss_total   = sum(aov_tab$`Sum Sq`, na.rm = TRUE)
          r2 = if (is.finite(ss_between) && ss_total > 0) ss_between/ss_total else NA_real_
          tibble::tibble(trait = tr, type = "categorical", effect = r2, p = aov_tab$`Pr(>F)`[1])
        })
      } else tibble::tibble(trait=character(), type=character(), effect=numeric(), p=numeric())

      trait_effects = dplyr::bind_rows(eff_cont, eff_cat) |>
        dplyr::mutate(
          trait_clean = gsub("^trait_", "", .data$trait),
          plot_effect = dplyr::if_else(.data$type=="continuous", abs(.data$effect), .data$effect),
          signif = .data$p < 0.05
        ) |>
        dplyr::arrange(dplyr::desc(.data$plot_effect))

      if (nrow(trait_effects) == 0) trait_effects = NULL
    }
  }


  # Plots (no boundary here; added by summarise_results())
  plots = list(site_map=NULL, invader_rank=NULL, heatmap=NULL, trait_effects=NULL)
  if (isTRUE(make_plots) && requireNamespace("ggplot2", quietly = TRUE)) {
    if (!is.null(site_df) && all(c("x","y") %in% names(site))) {
      plots$site_map =
        ggplot2::ggplot(site, ggplot2::aes(x, y, fill = .data$total_expected)) +
        ggplot2::geom_tile(color="grey70") +
        ggplot2::scale_fill_viridis_c(name = if (mode=="probabilistic") "Expected # establishing"
                                      else "# establishing") +
        ggplot2::labs(
          title = paste0("Invasibility map — ", mode, if (!is.null(label)) paste0(" (", label, ")") else ""),
          x = "Longitude", y = "Latitude"
        ) +
        ggplot2::theme_minimal()
    }

    plots$invader_rank =
      ggplot2::ggplot(species,
                      ggplot2::aes(x = stats::reorder(.data$invader, .data$V_i), y = .data$V_i, fill = .data$V_i)) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::percent_format()) +
      ggplot2::scale_fill_viridis_c(name="Mean across sites") +
      ggplot2::labs(
        x = "Invader", y = if (mode=="probabilistic") "Mean P(establish)" else "Share of sites (λ>0)",
        title = paste0("Invasiveness — ", mode, if (!is.null(label)) paste0(" (", label, ")") else "")
      ) +
      ggplot2::theme_minimal()

    # Heatmap (site × invader)
    d_heat =
      if (!is.null(establish_long)) establish_long[, c("site","invader","val")]
    else tibble::tibble(
      site = rep(rn, times = I),
      invader = rep(cn, each = S),
      val = as.vector(W)
    )

    # If the surface is binary, use discrete darkgrey (0) and darkred (1)
    is_binary = all(is.finite(d_heat$val)) && all(unique(d_heat$val) %in% c(0,1))

    if (is_binary) {
      d_heat$val_f = factor(d_heat$val, levels = c(0,1), labels = c("0","1"))
      plots$heatmap =
        ggplot2::ggplot(d_heat, ggplot2::aes(.data$invader, .data$site, fill = .data$val_f)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_manual(values = c("0" = "darkgrey", "1" = "darkred"),
                                   name = "Establish (0/1)") +
        ggplot2::labs(title = paste0("Establishment matrix — ", mode),
                      x = "Invader", y = "Site") +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
    } else {
      plots$heatmap =
        ggplot2::ggplot(d_heat, ggplot2::aes(.data$invader, .data$site, fill = .data$val)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(name = if (mode=="probabilistic") "P(establish)" else "Establish (0/1)") +
        ggplot2::labs(title = paste0("Establishment matrix — ", mode),
                      x = "Invader", y = "Site") +
        ggplot2::theme_minimal() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
    }

    # Trait-effects plot (unchanged)
    if (!is.null(trait_effects) && nrow(trait_effects)) {
      plots$trait_effects =
        ggplot2::ggplot(trait_effects,
                        ggplot2::aes(x = .data$plot_effect,
                                     y = forcats::fct_reorder(.data$trait_clean, .data$plot_effect))) +
        ggplot2::geom_segment(
          ggplot2::aes(x = 0, xend = .data$plot_effect,
                       yend = forcats::fct_reorder(.data$trait_clean, .data$plot_effect)),
          color = "grey70") +
        ggplot2::geom_point(ggplot2::aes(color = .data$type, shape = .data$signif), size = 3) +
        ggplot2::scale_x_continuous(name = "Effect size (|β| for continuous; R² for categorical)") +
        ggplot2::scale_color_manual(values = c(continuous="#1b9e77", categorical="#d95f02"), name="Trait type") +
        ggplot2::scale_shape_manual(values = c(`TRUE`=16, `FALSE`=1), name="p < 0.05") +
        ggplot2::labs(y = "Trait", title = "Trait invasiveness effects") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position="right")
    }

  }

  list(
    species         = species,
    site            = site,
    trait_effects   = trait_effects,
    establish_long  = establish_long,
    plots           = plots,
    meta = list(mode = mode, threshold = prob_threshold, label = label)
  )
}

