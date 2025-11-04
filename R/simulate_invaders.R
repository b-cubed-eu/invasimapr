#' Simulate hypothetical invader trait profiles from a resident trait pool
#'
#' @description
#' Generate \eqn{n_{\mathrm{inv}}} trait rows for hypothetical invaders by resampling
#' the empirical distribution of resident traits. By default, traits are sampled
#' independently per column (“columnwise”), creating novel combinations across
#' traits. Alternatively, entire rows can be bootstrapped (“rowwise”) to preserve
#' the resident covariance structure. Row names of the returned object are set to
#' the invader IDs.
#'
#' @details
#' **Species identifiers.** Supply species IDs in a dedicated column via
#' `species_col`, or set `species_col = NULL` to use `row.names(resident_traits)`
#' as the species IDs. In both cases, newly simulated invaders receive fresh,
#' unique IDs (`inv_prefix1`, `inv_prefix2`, …), which become the **row names**
#' of the returned `data.frame`. When `species_col` is not `NULL`, the same IDs
#' are also stored in that column (unless `keep_species_column = FALSE`).
#'
#' **Trait selection.** If `trait_cols` is `NULL`, all columns except
#' `species_col` (when present) are considered traits. Otherwise, only the
#' intersection of `trait_cols` and existing column names is used.
#'
#' **Sampling modes.**
#' - `mode = "columnwise"`: Each trait is generated independently. Numeric traits
#'   can be drawn by bootstrap (empirical sampling), uniform within observed
#'   bounds, or from a normal distribution parameterized by the empirical mean
#'   and SD (optionally truncated to observed min–max if `keep_bounds = TRUE`).
#'   Factor and character traits are sampled from their observed values/levels.
#' - `mode = "rowwise"`: Entire rows are resampled with replacement from
#'   `resident_traits`, preserving the joint structure. In this mode, if
#'   `species_col` is provided, its values are overwritten with the new invader
#'   IDs. If `species_col = NULL`, only the trait columns are returned and the
#'   row names are replaced by the invader IDs.
#'
#' **ID collisions.** If any proposed invader IDs would collide with existing
#' resident IDs, they are made unique using \code{\link[base]{make.unique}}.
#'
#' @param resident_traits A `data.frame` containing a species ID column
#'   (specified by `species_col`) **or** species IDs as row names, plus one or
#'   more trait columns (numeric, factor, character supported).
#' @param n_inv Integer; number of invaders to simulate (default `10`).
#' @param species_col `NULL` or character; the species ID column name in
#'   `resident_traits`. If `NULL`, species IDs are taken from row names
#'   (default `"species"`).
#' @param trait_cols `NULL` or character vector; which trait columns to use.
#'   Default: all columns except `species_col` (when present).
#' @param mode Either `"columnwise"` (new trait combinations; default) or
#'   `"rowwise"` (preserve covariance by resampling rows).
#' @param numeric_method For columnwise **numeric** traits: one of
#'   `"bootstrap"` (default; sample from observed values),
#'   `"normal"` (draws from \eqn{N(\bar{x}, s)}; truncated to [min, max] if
#'   `keep_bounds = TRUE`), or `"uniform"` (draws from `Uniform[min, max]`).
#' @param keep_bounds Logical; if `TRUE`, normal or uniform draws are constrained
#'   to the observed `[min, max]`. This does not apply to `"bootstrap"`
#'   (default `TRUE`).
#' @param inv_prefix Character; prefix used to construct invader IDs
#'   (default `"inv"`).
#' @param keep_species_column Logical; when `species_col` is not `NULL`, keep
#'   the species ID column after setting row names (default `TRUE`). Ignored if
#'   `species_col = NULL`.
#' @param seed `NULL` or integer; optional RNG seed for reproducibility.
#'
#' @return A `data.frame` of simulated invaders. Row names are the invader IDs.
#'   If `species_col` is not `NULL` and `keep_species_column = TRUE`, that
#'   column will contain the same IDs as the row names.
#'
#' @examples
#' ## ---------------------------
#' ## Example 1: species IDs in a column
#' ## ---------------------------
#' set.seed(1)
#' residents_col = data.frame(
#'   species = paste0("sp", 1:5),
#'   height  = c(10.2, 9.8, 11.1, 10.5, 9.9),
#'   SLA     = c(15.0, 15.2, 14.7, 15.5, 15.1),
#'   lifeform = factor(c("tree", "shrub", "shrub", "tree", "herb"))
#' )
#'
#' # Columnwise bootstrap for numeric traits; factor sampled from levels.
#' inv1 = simulate_invaders(
#'   residents_col, n_inv = 4, species_col = "species",
#'   mode = "columnwise", numeric_method = "bootstrap",
#'   inv_prefix = "inv", keep_species_column = TRUE, seed = 42
#' )
#' head(inv1)
#'
#' # Rowwise resampling preserves covariance (entire rows).
#' inv2 = simulate_invaders(
#'   residents_col, n_inv = 3, species_col = "species",
#'   mode = "rowwise", inv_prefix = "inv", seed = 123
#' )
#' head(inv2)
#'
#' ## ---------------------------
#' ## Example 2: species IDs in row names (species_col = NULL)
#' ## ---------------------------
#' set.seed(2)
#' residents_rn = data.frame(
#'   height  = c(10.2, 9.8, 11.1, 10.5, 9.9),
#'   SLA     = c(15.0, 15.2, 14.7, 15.5, 15.1),
#'   lifeform = factor(c("tree", "shrub", "shrub", "tree", "herb"))
#' )
#' rownames(residents_rn) = paste0("sp", 1:5)
#'
#' # Columnwise with normal draws truncated to observed min–max
#' inv3 = simulate_invaders(
#'   residents_rn, n_inv = 5, species_col = NULL,
#'   mode = "columnwise", numeric_method = "normal",
#'   keep_bounds = TRUE, inv_prefix = "x", seed = 99
#' )
#' head(inv3)  # row names are the invader IDs; no species column present
#'
#' # Rowwise with species_col = NULL: entire rows resampled; row names replaced
#' inv4 = simulate_invaders(
#'   residents_rn, n_inv = 3, species_col = NULL,
#'   mode = "rowwise", inv_prefix = "new", seed = 77
#' )
#' head(inv4)
#'
#' @seealso [base::sample], [base::make.unique], [stats::rnorm], [stats::runif], [stats::sd]
#' @keywords simulation resampling bootstrap traits ecology
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
