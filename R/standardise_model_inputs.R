# ======================================================================
# standardise_model_inputs.R
# ======================================================================

#' Standardise model inputs (no leakage) for residents and invaders
#'
#' @description
#' Column-wise z-scores **environment** and **resident trait** numerics, then
#' scales **invader trait** numerics **using resident moments only** (to avoid
#' information leakage). Invader factor/character columns are coerced to the
#' resident levels; unseen levels become `NA`. Optionally drops invader-only
#' columns so the resident/invader trait schemas match.
#'
#' @details
#' **What gets standardised and how**
#' - **Environment (`env_df`)**: numeric columns are z-scored (mean 0, sd 1);
#'   non-numeric columns are kept as-is. Zero variance is guarded by setting sd=1.
#' - **Resident traits (`traits_res`)**: numeric columns are z-scored; mixed
#'   types allowed—non-numeric columns are kept.
#' - **Invader traits (`traits_inv`)**: numeric columns are scaled using the
#'   **resident** trait means/SDs only (never computed from invaders). Factor/
#'   character columns are coerced to resident levels; unseen levels become `NA`.
#'   Extra invader columns are dropped (with a note).
#'
#' **Returned objects**
#' - `env_df_z`: environment table with numeric columns standardised (or `NULL`)
#' - `traits_res_glmm`: resident traits with numeric columns standardised
#' - `traits_inv_glmm`: invader traits, scaled **like residents** + factor levels matched (or `NULL`)
#' - `moments`: resident/reference moments used for scaling (`env_*`, `trait_*`)
#' - `info$notes`: human-readable notes on coercions/drops
#'
#' **Where this is used in the workflow**
#' - Called explicitly prior to GLMM fitting and when harmonising invaders,
#'   and implicitly by wrappers such as `prepare_trait_space()` (if available).
#'
#' @param env_df   Optional `data.frame` (sites × environment). Numeric columns
#'   are z-scored; non-numeric are preserved.
#' @param traits_res `data.frame` (residents × traits). Mixed types allowed;
#'   numeric columns are z-scored.
#' @param traits_inv Optional `data.frame` (invaders × traits). Must contain at
#'   least the trait columns present in `traits_res`. Numeric columns are scaled
#'   using **resident** means/SDs; factors are coerced to resident levels.
#' @param drop_extra_invader_cols Logical; if `TRUE`, invader-only columns are
#'   dropped (not used downstream). If `FALSE`, they are still dropped for
#'   alignment but flagged in `info$notes`.
#' @param verbose Logical; print messages about what was standardised/coerced.
#'
#' @return A named `list` with components:
#' \describe{
#'   \item{env_df_z}{Environment table with numeric columns z-scored (or `NULL`).}
#'   \item{traits_res_glmm}{Resident trait table with numeric columns z-scored.}
#'   \item{traits_inv_glmm}{Invader trait table scaled to resident moments and factor levels matched (or `NULL`).}
#'   \item{moments}{`list(env_means, env_sds, trait_means_res, trait_sds_res)` used for scaling.}
#'   \item{info}{`list(notes=character())` with human-readable notes.}
#' }
#'
#' @section Invariants and guards:
#' - Column names and row names are preserved.
#' - Zero-variance numeric columns use `sd=1` so z-scores stay defined.
#' - Invader trait numerics are always scaled by **resident** moments (no leakage).
#' - Invader extra columns are dropped for alignment; missing required columns error.
#'
#' @seealso
#' `prepare_inputs()`, `prepare_trait_space()`, `build_invader_predictors()`,
#' `compute_trait_space()`, `compute_resident_crowding()`
#'
#' @examples
#' # Minimal reproducible example ----------------------------------------------
#' set.seed(1)
#' env_df = data.frame(elev = rnorm(5), temp = rnorm(5), zone = factor(sample(c("A","B"), 5, TRUE)))
#' rownames(env_df) = paste0("s", 1:5)
#'
#' traits_res = data.frame(
#'   size = rlnorm(4), leaf = factor(c("broad","needle","broad","needle")),
#'   stringsAsFactors = FALSE
#' )
#' rownames(traits_res) = paste0("sp", 1:4)
#'
#' traits_inv = data.frame(
#'   size = c(10, 1), leaf = factor(c("broad","unknown"))  # 'unknown' -> NA after coercion
#' )
#' rownames(traits_inv) = c("inv1","inv2")
#'
#' std = standardise_model_inputs(env_df, traits_res, traits_inv, verbose = FALSE)
#' str(std, 1)
#' head(std$traits_res_glmm)
#' head(std$traits_inv_glmm)   # note: 'leaf' for 'unknown' becomes NA
#'
#' # Workflow sketch ------------------------------------------------------------
#' # fit = prepare_inputs(long_df = longDF, ...)         # gives fit$inputs$env_df, $traits_res
#' # inv  = simulate_invaders(fit$inputs$traits_res, n_inv = 10)
#' # std  = standardise_model_inputs(fit$inputs$env_df, fit$inputs$traits_res, inv)
#' # std$traits_res_glmm; std$traits_inv_glmm              # pass to GLMM / trait-space steps
#'
#' @export
standardise_model_inputs = function(env_df = NULL,
                                     traits_res,
                                     traits_inv = NULL,
                                     drop_extra_invader_cols = FALSE,
                                     verbose = TRUE) {
  stopifnot(is.data.frame(traits_res))

  # ---- small internal helpers (safe if utils not loaded) ----------------------
  .std_df_local = function(df, ref_means = NULL, ref_sds = NULL) {
    stopifnot(is.data.frame(df))
    out = df
    num = vapply(df, is.numeric, logical(1))
    if (!any(num)) {
      return(list(data = out, means = numeric(0), sds = numeric(0)))
    }
    X = as.matrix(df[, num, drop = FALSE])
    if (is.null(ref_means) || is.null(ref_sds)) {
      m = colMeans(X, na.rm = TRUE)
      s = apply(X, 2, stats::sd, na.rm = TRUE)
    } else {
      m = ref_means[colnames(X)]
      s = ref_sds  [colnames(X)]
    }
    s[!is.finite(s) | s == 0] = 1
    Z = sweep(sweep(X, 2, m, "-"), 2, s, "/")
    out[, num] = Z
    list(data = out, means = m, sds = s)
  }

  .std_df = if (exists(".standardise_df", mode = "function")) .standardise_df else .std_df_local

  .coerce_invader_factors = function(inv_df, res_df, notes) {
    for (nm in intersect(colnames(res_df), colnames(inv_df))) {
      rcol = res_df[[nm]]
      icol = inv_df[[nm]]
      if (is.factor(rcol) || is.character(rcol)) {
        lev = if (is.factor(rcol)) levels(rcol) else sort(unique(as.character(rcol)))
        if (is.factor(icol)) icol = as.character(icol)
        inv_df[[nm]] = factor(as.character(icol), levels = lev)
        unseen = setdiff(unique(as.character(icol)), lev)
        if (length(unseen)) {
          notes = c(notes, sprintf(
            "Coerced invader levels in '%s' to resident levels; unseen -> NA: %s",
            nm, paste(unseen, collapse = ", ")
          ))
        }
      }
    }
    list(inv_df = inv_df, notes = notes)
  }

  # ---- 1) ENV standardisation (optional) --------------------------------------
  env_df_z = NULL
  env_means = env_sds = numeric(0)
  if (!is.null(env_df)) {
    stopifnot(is.data.frame(env_df))
    env_std = .std_df(as.data.frame(env_df))
    env_df_z = env_std$data
    env_means = env_std$means
    env_sds   = env_std$sds
    if (isTRUE(verbose)) {
      message(sprintf(
        "ENV: z-scored %d numeric column(s); kept %d non-numeric.",
        length(env_means), ncol(env_df) - length(env_means)
      ))
    }
  }

  # ---- 2) RESIDENT trait standardisation --------------------------------------
  tr_res = as.data.frame(traits_res)
  tr_res_std = .std_df(tr_res)
  traits_res_glmm = tr_res_std$data
  trait_means_res = tr_res_std$means
  trait_sds_res   = tr_res_std$sds
  if (isTRUE(verbose)) {
    message(sprintf("TRAITS (residents): z-scored %d numeric column(s).",
                    length(trait_means_res)))
  }

  notes = character()

  # ---- 3) INVADER trait alignment & standardisation (no leakage) --------------
  traits_inv_glmm = NULL
  if (!is.null(traits_inv)) {
    inv = as.data.frame(traits_inv)

    # (a) enforce same columns as residents
    res_cols = colnames(traits_res_glmm)
    inv_cols = colnames(inv)

    extra_inv = setdiff(inv_cols, res_cols)
    miss_inv  = setdiff(res_cols, inv_cols)

    if (length(extra_inv)) {
      notes = c(notes, sprintf(
        "Invader has %d extra trait column(s) not in residents: %s",
        length(extra_inv), paste(extra_inv, collapse = ", ")
      ))
      # Must align to resident schema for downstream steps
      inv = inv[, setdiff(inv_cols, extra_inv), drop = FALSE]
      notes = c(notes, if (drop_extra_invader_cols)
        "Dropped invader-only columns (drop_extra_invader_cols = TRUE)."
        else
          "Dropped invader-only columns to align with residents (required downstream).")
    }

    if (length(miss_inv)) {
      stop("traits_inv is missing required resident trait columns: ",
           paste(miss_inv, collapse = ", "),
           "\nAdd these to traits_inv (even if NA) or re-simulate invaders.")
    }

    inv = inv[, res_cols, drop = FALSE]

    # (b) factor/character coercion to resident levels
    co = .coerce_invader_factors(inv, traits_res, notes)
    inv = co$inv_df
    notes = co$notes

    # (c) numeric scaling using **resident** moments only
    num_mask_res = vapply(traits_res_glmm, is.numeric, logical(1))
    num_cols = names(which(num_mask_res))
    if (length(num_cols)) {
      Xinv = inv[, num_cols, drop = FALSE]
      ref_m = trait_means_res[num_cols]; ref_s = trait_sds_res[num_cols]
      ref_s[!is.finite(ref_s) | ref_s == 0] = 1
      Zin  = sweep(sweep(as.matrix(Xinv), 2, ref_m, "-"), 2, ref_s, "/")
      inv[, num_cols] = Zin
    }

    traits_inv_glmm = inv

    if (isTRUE(verbose)) {
      message(sprintf(
        "TRAITS (invaders): scaled %d numeric column(s) using resident moments; factors coerced.",
        length(num_cols)
      ))
    }
  }

  # ---- 4) Return --------------------------------------------------------------
  list(
    env_df_z        = env_df_z,
    traits_res_glmm = traits_res_glmm,
    traits_inv_glmm = traits_inv_glmm,
    moments = list(
      env_means        = env_means,
      env_sds          = env_sds,
      trait_means_res  = trait_means_res,
      trait_sds_res    = trait_sds_res
    ),
    info = list(
      notes = notes
    )
  )
}
