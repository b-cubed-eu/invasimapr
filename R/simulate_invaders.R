#' Simulate hypothetical invader trait profiles from a resident trait pool
#'
#' @description
#' Generates trait rows for hypothetical invaders by resampling the empirical
#' distribution of resident traits. Row names are set to the invader IDs.
#'
#' @param resident_traits data.frame with a species ID column and trait columns.
#' @param n_inv integer; number of invaders to simulate.
#' @param species_col character; species ID column name in resident_traits.
#' @param trait_cols NULL or character; which trait columns to use (default: all except species_col).
#' @param mode "columnwise" (new combinations; default) or "rowwise" (preserve covariance).
#' @param numeric_method for columnwise numeric traits: "bootstrap" (default), "normal", or "uniform".
#' @param keep_bounds logical; constrain normal/uniform draws to observed [min,max] (default TRUE).
#' @param inv_prefix character; prefix for invader IDs (default "inv").
#' @param keep_species_column logical; keep species ID column after setting row names (default TRUE).
#' @param seed NULL or integer; RNG seed.
#'
#' @return data.frame of simulated invaders with row names == species IDs.
#' @export
simulate_invaders = function(
    resident_traits,
    n_inv              = 10,
    species_col        = "species",
    trait_cols         = NULL,
    mode               = c("columnwise", "rowwise"),
    numeric_method     = c("bootstrap", "normal", "uniform"),
    keep_bounds        = TRUE,
    inv_prefix         = "inv",
    keep_species_column= TRUE,
    seed               = NULL
){
  stopifnot(is.data.frame(resident_traits))
  if (!species_col %in% names(resident_traits)) stop("'", species_col, "' not found.")
  mode = match.arg(mode)
  numeric_method = match.arg(numeric_method)
  if (is.null(trait_cols)) trait_cols = setdiff(names(resident_traits), species_col)
  trait_cols = intersect(trait_cols, names(resident_traits))
  if (!length(trait_cols)) stop("No trait columns detected.")
  if (!is.null(seed)) set.seed(seed)

  # --- build invader ID vector, ensuring uniqueness vs residents -------------
  existing_ids = as.character(resident_traits[[species_col]])
  proposed_ids = paste0(inv_prefix, seq_len(n_inv))
  if (length(intersect(existing_ids, proposed_ids))) {
    # make IDs unique relative to residents and within themselves
    proposed_ids = make.unique(c(existing_ids, proposed_ids), sep = "_dup_")
    proposed_ids = tail(proposed_ids, n_inv)
  }

  # --- helper for truncated normal ------------------------------------------
  draw_trunc_normal = function(mu, sd, a, b, n) {
    if (is.na(sd) || sd == 0) return(rep(mu, n))
    out = numeric(0L)
    while (length(out) < n) {
      z = stats::rnorm(n, mean = mu, sd = sd)
      z = z[z >= a & z <= b]
      out = c(out, z)
    }
    out[seq_len(n)]
  }

  if (mode == "rowwise") {
    idx = sample(seq_len(nrow(resident_traits)), size = n_inv, replace = TRUE)
    inv = resident_traits[idx, c(species_col, trait_cols), drop = FALSE]
    inv[[species_col]] = proposed_ids
  } else {
    sim_list = lapply(trait_cols, function(col) {
      x = resident_traits[[col]]
      if (is.factor(x)) {
        factor(sample(levels(x), n_inv, replace = TRUE), levels = levels(x), ordered = is.ordered(x))
      } else if (is.character(x)) {
        vals = unique(x)
        factor(sample(vals, n_inv, replace = TRUE), levels = vals)
      } else if (is.numeric(x)) {
        xmin = min(x, na.rm = TRUE); xmax = max(x, na.rm = TRUE)
        if (numeric_method == "bootstrap") {
          sample(x, n_inv, replace = TRUE)
        } else if (numeric_method == "uniform") {
          stats::runif(n_inv, xmin, xmax)
        } else {
          mu = mean(x, na.rm = TRUE); sd = stats::sd(x, na.rm = TRUE)
          if (keep_bounds) draw_trunc_normal(mu, sd, xmin, xmax, n_inv) else stats::rnorm(n_inv, mu, sd)
        }
      } else {
        sample(x, n_inv, replace = TRUE)
      }
    })
    names(sim_list) = trait_cols
    inv = as.data.frame(sim_list, stringsAsFactors = FALSE, check.names = FALSE)
    # restore factor structures to match residents
    for (col in trait_cols) {
      res_col = resident_traits[[col]]
      if (is.factor(res_col) && !is.factor(inv[[col]])) {
        inv[[col]] = factor(inv[[col]], levels = levels(res_col), ordered = is.ordered(res_col))
      }
    }
    inv[[species_col]] = proposed_ids
    inv = inv[, c(species_col, trait_cols), drop = FALSE]
  }

  # --- set row names to species IDs ------------------------------------------
  rownames(inv) = as.character(inv[[species_col]])
  if (!keep_species_column) inv[[species_col]] = NULL

  inv
}
