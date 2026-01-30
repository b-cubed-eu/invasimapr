#' Simulate hypothetical invader trait profiles from a resident trait pool
#'
#' @description
#' Generate \eqn{n_{\mathrm{inv}}} rows of trait values for hypothetical invaders
#' by resampling the empirical distribution of resident traits. By default,
#' traits are sampled independently by column, creating novel trait
#' combinations. Alternatively, entire rows can be resampled to preserve the
#' covariance structure of the resident trait space.
#'
#' Row names of the returned object are set to the invader identifiers.
#'
#' @details
#' \strong{Species identifiers}
#' Species identifiers can be supplied via a dedicated column specified by
#' `species_col`, or taken from the row names of `resident_traits` when
#' `species_col = NULL`. Newly simulated invaders receive fresh, unique
#' identifiers constructed from `inv_prefix`, which become the row names of
#' the returned data frame. When `species_col` is not `NULL`, the same identifiers
#' are also stored in that column unless `keep_species_column = FALSE`.
#'
#' \strong{Trait selection}
#' If `trait_cols` is `NULL`, all columns except `species_col` (when present)
#' are treated as traits. Otherwise, only the intersection of `trait_cols` and
#' existing column names is used.
#'
#' \strong{Sampling modes}
#' \itemize{
#'   \item \code{mode = "columnwise"}: Each trait is generated independently.
#'         Numeric traits can be drawn by bootstrap sampling, from a normal
#'         distribution parameterised by the empirical mean and standard
#'         deviation, or from a uniform distribution on the observed range.
#'         Factor and character traits are sampled from observed values.
#'   \item \code{mode = "rowwise"}: Entire rows are resampled with replacement
#'         from the resident trait table, preserving joint structure across
#'         traits.
#' }
#'
#' \strong{Identifier collisions}
#' If proposed invader identifiers collide with existing resident identifiers,
#' they are made unique using \link[base]{make.unique}.
#'
#' @param resident_traits A data frame containing either a species identifier
#'   column (specified by `species_col`) or species identifiers as row names,
#'   plus one or more trait columns.
#' @param n_inv Integer; number of invaders to simulate (default `10`).
#' @param species_col `NULL` or character; name of the species identifier column
#'   in `resident_traits`. If `NULL`, species identifiers are taken from row names.
#' @param trait_cols `NULL` or character vector specifying which trait columns
#'   to use. Defaults to all columns except `species_col` when present.
#' @param mode Character; either `"columnwise"` or `"rowwise"`.
#' @param numeric_method Character; for numeric traits in columnwise mode,
#'   one of `"bootstrap"`, `"normal"`, or `"uniform"`.
#' @param keep_bounds Logical; if `TRUE`, normal or uniform draws are constrained
#'   to the observed minimum and maximum values.
#' @param inv_prefix Character; prefix used to construct invader identifiers.
#' @param keep_species_column Logical; if `TRUE` and `species_col` is not `NULL`,
#'   retain the species identifier column after setting row names.
#' @param seed `NULL` or integer; optional random seed for reproducibility.
#'
#' @return
#' A data frame of simulated invader traits. Row names correspond to invader
#' identifiers. If `species_col` is not `NULL` and `keep_species_column = TRUE`,
#' the species identifier column contains the same identifiers.
#'
#' @seealso
#' \link[base]{sample},
#' \link[base]{make.unique},
#' \link[stats]{rnorm},
#' \link[stats]{runif},
#' \link[stats]{sd}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' residents <- data.frame(
#'   species  = paste0("sp", 1:5),
#'   height   = c(10.2, 9.8, 11.1, 10.5, 9.9),
#'   SLA      = c(15.0, 15.2, 14.7, 15.5, 15.1),
#'   lifeform = factor(c("tree", "shrub", "shrub", "tree", "herb"))
#' )
#'
#' inv <- simulate_invaders(
#'   residents,
#'   n_inv = 4,
#'   species_col = "species",
#'   mode = "columnwise",
#'   numeric_method = "bootstrap"
#' )
#' head(inv)
#' }
#'
#' @keywords simulation resampling bootstrap traits
#' @importFrom stats rnorm runif sd
#' @export
simulate_invaders = function(
    resident_traits,
    n_inv = 10,
    species_col = "species",
    trait_cols = NULL,
    mode = c("columnwise", "rowwise"),
    numeric_method = c("bootstrap", "normal", "uniform"),
    keep_bounds = TRUE,
    inv_prefix = "inv",
    keep_species_column = TRUE,
    seed = NULL) {

  stopifnot(is.data.frame(resident_traits))
  mode           = match.arg(mode)
  numeric_method = match.arg(numeric_method)

  # --- Guard against empty resident table -------------------------------------
  if (nrow(resident_traits) < 1L) {
    stop("`resident_traits` has 0 rows; nothing to resample. ",
         "Check your filtering/inputs.")
  }

  # --- Species IDs -------------------------------------------------------------
  if (is.null(species_col)) {
    existing_ids = rownames(resident_traits)
    if (is.null(existing_ids)) stop("No `species_col` and row names are NULL.")
  } else {
    if (!species_col %in% names(resident_traits)) {
      stop("'", species_col, "' not found in resident_traits.")
    }
    existing_ids = as.character(resident_traits[[species_col]])
  }

  # --- Trait columns -----------------------------------------------------------
  if (is.null(trait_cols)) {
    trait_cols = if (is.null(species_col)) names(resident_traits) else
      setdiff(names(resident_traits), species_col)
  }
  trait_cols = intersect(trait_cols, names(resident_traits))
  if (!length(trait_cols)) stop("No trait columns detected.")

  if (!is.null(seed)) set.seed(seed)

  # --- New invader IDs (unique vs. residents) ---------------------------------
  proposed_ids = paste0(inv_prefix, seq_len(n_inv))
  if (length(intersect(existing_ids, proposed_ids))) {
    proposed_ids = make.unique(c(existing_ids, proposed_ids), sep = "_dup_")
    proposed_ids = tail(proposed_ids, n_inv)
  }

  # --- Helper: truncated normal draw ------------------------------------------
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

  # --- Helper: safe sampling domain -------------------------------------------
  safe_sample_vals = function(vals, n, as_factor_levels = NULL, ordered = FALSE) {
    # Remove NAs for sampling domain
    vals_nonna = vals[!is.na(vals)]
    if (length(vals_nonna) < 1L) {
      # Nothing to sample from → return NAs (as factor if requested)
      if (is.null(as_factor_levels)) return(rep(NA, n))
      return(factor(rep(NA_character_, n),
                    levels = as_factor_levels, ordered = ordered))
    }
    smp = sample(vals_nonna, n, replace = TRUE)
    if (is.null(as_factor_levels)) return(smp)
    factor(smp, levels = as_factor_levels, ordered = ordered)
  }

  # --- Simulation core ---------------------------------------------------------
  if (mode == "rowwise") {
    # Resample entire rows to preserve covariance.
    idx = sample(seq_len(nrow(resident_traits)), size = n_inv, replace = TRUE)
    if (is.null(species_col)) {
      inv = resident_traits[idx, trait_cols, drop = FALSE]
    } else {
      inv = resident_traits[idx, c(species_col, trait_cols), drop = FALSE]
      inv[[species_col]] = proposed_ids
    }

  } else {
    # Columnwise: independently sample each trait column with robust guards.
    sim_list = lapply(trait_cols, function(col) {
      x = resident_traits[[col]]

      if (is.factor(x)) {
        # Sample from factor levels; if zero levels, fall back to non-NA values.
        levs = levels(x)
        if (length(levs) < 1L) {
          # No levels defined; use unique non-NA character values as pseudo-levels
          vals = unique(as.character(x))
          vals = vals[!is.na(vals)]
          if (length(vals) < 1L) {
            return(factor(rep(NA_character_, n_inv), levels = character(0),
                          ordered = is.ordered(x)))
          }
          return(factor(sample(vals, n_inv, replace = TRUE),
                        levels = vals, ordered = is.ordered(x)))
        } else {
          return(safe_sample_vals(levs, n_inv,
                                  as_factor_levels = levs,
                                  ordered = is.ordered(x)))
        }

      } else if (is.character(x)) {
        vals = unique(x)
        return(safe_sample_vals(vals, n_inv, as_factor_levels = vals))

      } else if (is.numeric(x)) {
        xmin = suppressWarnings(min(x, na.rm = TRUE))
        xmax = suppressWarnings(max(x, na.rm = TRUE))
        if (!is.finite(xmin) || !is.finite(xmax)) {
          # All NA numeric column → return NA numerics
          return(rep(NA_real_, n_inv))
        }
        if (numeric_method == "bootstrap") {
          # Bootstrap from non-NA values only
          return(safe_sample_vals(x, n_inv))
        } else if (numeric_method == "uniform") {
          return(stats::runif(n_inv, xmin, xmax))
        } else {
          mu = mean(x, na.rm = TRUE)
          sd = stats::sd(x, na.rm = TRUE)
          if (keep_bounds) {
            return(draw_trunc_normal(mu, sd, xmin, xmax, n_inv))
          } else {
            return(stats::rnorm(n_inv, mu, sd))
          }
        }

      } else {
        # Fallback for other types (e.g., Date) → bootstrap non-NA
        return(safe_sample_vals(x, n_inv))
      }
    })
    names(sim_list) = trait_cols
    inv = as.data.frame(sim_list, stringsAsFactors = FALSE, check.names = FALSE)

    # Restore factor structure to match resident columns when needed.
    for (col in trait_cols) {
      res_col = resident_traits[[col]]
      if (is.factor(res_col) && !is.factor(inv[[col]])) {
        inv[[col]] = factor(inv[[col]],
                             levels = levels(res_col),
                             ordered = is.ordered(res_col))
      }
    }

    if (!is.null(species_col)) {
      inv[[species_col]] = proposed_ids
      inv = inv[, c(species_col, trait_cols), drop = FALSE]
    }
  }

  # --- Set row names; optionally drop species column --------------------------
  if (is.null(species_col)) {
    rownames(inv) = proposed_ids
  } else {
    rownames(inv) = as.character(inv[[species_col]])
    if (!keep_species_column) inv[[species_col]] = NULL
  }

  inv
}
