#' Build standardized invader predictors (r_is_z, C_is_z, S_is_z)
#'
#' Uses resident-side moments and the residents-only model to construct
#' invader-level predictors on the **site-standardised** scale:
#' \eqn{r^{(z)}_{is}}, \eqn{C^{(z)}_{is}}, \eqn{S^{(z)}_{is}}.
#' Handles PCA projection with “dummy” factor expansion so that invader
#' covariates exactly match the **training design** used in
#' [`model_residents()`] / the stored resident GLMM.
#'
#' @section Rationale:
#' Many downstream steps (e.g. invasion fitness, establishment probability)
#' assume invader predictors live on the **same centred/scaled basis** as the
#' resident training data. This function guarantees alignment by: (i) coercing
#' raw trait factor levels to training levels (with `_other_` fallbacks), (ii)
#' rebuilding the design matrix with the training dummy map, (iii) applying
#' the stored centring/scaling, and (iv) projecting through the stored PCA
#' rotations before calling [`build_invader_predictors()`].
#'
#' @param fit `invasimapr_fit` after
#'   [`prepare_trait_space()`] / [`learn_sensitivities()`] and resident model
#'   fitting (i.e., it contains resident GLMM artefacts, PCA objects/metadata,
#'   and resident moments such as `r_mu_s`, `r_sd_s`, `C_mu_s`, `C_sd_s`).
#' @param traits_inv **data.frame** of *raw* invader trait values; columns must
#'   match the resident trait names used at training (numerics and/or factors).
#'
#' @return The input `fit` with an updated `fit$invaders` list containing:
#' \describe{
#'   \item{`traits_inv_raw`}{Raw (unmodified) invader trait table (as supplied).}
#'   \item{`traits_inv_glmm`}{Design-scale trait data passed to the resident model
#'         (post dummy-expansion, standardisation, PCA, and factor alignment).}
#'   \item{`r_is_z`, `C_is_z`, `S_is_z`}{Site-standardised matrices for invaders.}
#'   \item{`df`}{Tidy table used/returned by [`build_invader_predictors()`] for joins.}
#' }
#'
#' @details
#' **Pipeline**
#' 1. **Environment**: if an environment PCA was used at training, rebuild the
#'    environment design matrix with the saved factor map and project using the
#'    stored rotation; otherwise reuse the stored `env_df_z`.
#' 2. **Traits**: coerce invader raw trait factors to training levels (redirect
#'    unseen levels to `_other_` if present); rebuild the dummy matrix and
#'    standardise with stored means/SDs; project via trait PCA if present; then
#'    append any raw factor terms that remained in the fixed-effects formula.
#' 3. **Predictors**: call
#'    \href{https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.html}{`build_invader_predictors()`}
#'    to compute \eqn{r^{(z)}_{is}}, \eqn{C^{(z)}_{is}}, \eqn{S^{(z)}_{is}} using
#'    resident moments, crowding kernels, and similarity structures stored in `fit`.
#'
#' **Assumptions & safeguards**
#' - Requires a trained resident model stored in `fit$residents$fit_r` and PCA
#'   metadata in `fit$model$*_pca*` if PCA was used.
#' - Uses resident moments (`r_mu_s`, `r_sd_s`, `C_mu_s`, `C_sd_s`) and
#'   similarity/crowding information (`W_site`, `gower`).
#' - Site ordering is taken from `fit$meta$sites`; inputs are conformed to it.
#'
#' @seealso
#' - Input assembly: \href{https://b-cubed-eu.github.io/invasimapr/reference/assemble_matrices.html}{`assemble_matrices()`}
#' - Predictor builder: \href{https://b-cubed-eu.github.io/invasimapr/reference/build_invader_predictors.html}{`build_invader_predictors()`}
#' - Residents GLMM: \href{https://b-cubed-eu.github.io/invasimapr/reference/fit_auxiliary_residents_glmm.html}{`fit_auxiliary_residents_glmm()`}
#' - Sensitivities: \href{https://b-cubed-eu.github.io/invasimapr/reference/learn_sensitivities.html}{`learn_sensitivities()`}
#'
#' @examples
#' \dontrun{
#' fit <- prepare_inputs(sites = site_df, residents = resident_df,
#'                       invaders = invader_df, traits = trait_df)
#' fit <- learn_sensitivities(fit)
#'
#' # New invader trait table (raw scale, same columns as residents' traits)
#' inv_traits <- data.frame(height = c(1.3, 0.9), SLA = c(12, 18),
#'                          life_form = c("shrub","forb"),
#'                          row.names = c("spA","spB"))
#'
#' fit <- predict_invaders(fit, inv_traits)
#' str(fit$invaders, 1)
#' dim(fit$invaders$r_is_z)  # |sites| × |invaders|
#' }
#'
#' @export
predict_invaders = function(fit, traits_inv) {

  # --- Preconditions -----------------------------------------------------------
  stopifnot(inherits(fit, "invasimapr_fit"))
  `%||%` = function(a, b) if (!is.null(a)) a else b

  # ---------- helpers paired with model_residents() ----------------------------
  # Align raw factor levels to the training model frame (redirect unseen to `_other_`).
  .coerce_factor_levels = function(df, model_frame) {
    mf = model_frame
    common = intersect(names(df), names(mf))
    for (v in common) {
      if (is.factor(mf[[v]])) {
        lv = levels(mf[[v]])
        x  = as.character(df[[v]])
        if ("_other_" %in% lv) x[!x %in% lv] = "_other_"
        df[[v]] = factor(x, levels = lv, ordered = is.ordered(mf[[v]]))
      }
    }
    df
  }

  # Construct a dummy matrix using the saved factor map; conform columns to train_vars.
  .dummy_matrix_from_info = function(df_raw, vars_info, train_vars) {
    df = as.data.frame(df_raw, stringsAsFactors = FALSE)
    rn = rownames(df)
    for (nm in names(df)) if (is.character(df[[nm]])) df[[nm]] = factor(df[[nm]])

    # numeric block
    num_idx = names(df)[vapply(df, is.numeric, logical(1))]
    X_num = if (length(num_idx)) as.matrix(df[, num_idx, drop = FALSE]) else
      matrix(numeric(0), nrow = nrow(df), ncol = 0, dimnames = list(rn, NULL))

    # factor block via training map
    X_fac = if (!is.null(vars_info$factors) && length(vars_info$factors)) {
      MM_list = lapply(names(vars_info$factors), function(v) {
        if (!v %in% names(df)) df[[v]] = NA
        x  = as.character(df[[v]])
        lv = as.character(vars_info$factors[[v]])
        if ("_other_" %in% lv) x[!x %in% lv] = "_other_"
        f = factor(x, levels = lv)
        mm = stats::model.matrix(~ . - 1, data = data.frame(f = f))
        colnames(mm) = paste0(v, lv)  # matches model_residents() naming
        mm
      })
      do.call(cbind, MM_list)
    } else {
      matrix(numeric(0), nrow = nrow(df), ncol = 0, dimnames = list(rn, NULL))
    }

    # bind and conform to training variable order
    X = if (ncol(X_num) && ncol(X_fac)) cbind(X_num, X_fac) else if (ncol(X_num)) X_num else X_fac
    if (length(train_vars)) {
      miss = setdiff(train_vars, colnames(X))
      if (length(miss)) {
        add = matrix(0, nrow = nrow(X), ncol = length(miss),
                     dimnames = list(rownames(X), miss))
        X = if (ncol(X)) cbind(X, add) else add
      }
      X = X[, train_vars, drop = FALSE]
    } else {
      X = matrix(numeric(0), nrow = nrow(df), ncol = 0, dimnames = list(rn, NULL))
    }
    X
  }

  # Standardise columns like training (centre/scale vectors named by train_vars).
  .standardize_like_training = function(X, center, scale) {
    if (!length(center) && !length(scale)) return(X)
    for (j in seq_along(center)) {
      nm = names(center)[j]
      if (!nm %in% colnames(X)) next
      x = X[, nm]
      x[!is.finite(x)] = center[j]                # ensures 0 after standardisation
      sc = if (is.finite(scale[j]) && scale[j] != 0) scale[j] else 1
      X[, nm] = (x - center[j]) / sc
    }
    X
  }

  # Project a standardised design matrix through a stored prcomp rotation.
  .predict_pcs = function(pca, X_std, pc_prefix) {
    if (is.null(pca) || !inherits(pca, "prcomp"))
      return(data.frame(row.names = rownames(X_std)))
    rot_names = rownames(pca$rotation)
    if (is.null(rot_names)) stop("Stored PCA rotation has no rownames; cannot project.")
    miss = setdiff(rot_names, colnames(X_std))
    if (length(miss)) {
      add = matrix(0, nrow = nrow(X_std), ncol = length(miss),
                   dimnames = list(rownames(X_std), miss))
      X_std = cbind(X_std, add)
    }
    X_std = X_std[, rot_names, drop = FALSE]
    S = X_std %*% pca$rotation
    S = as.matrix(S)
    colnames(S) = paste0(pc_prefix, "PC", seq_len(ncol(S)))
    as.data.frame(S)
  }

  # ---------- Prepare environment PCs for all sites ----------------------------
  if (!is.null(fit$model$env_pca)) {
    env_raw = as.data.frame(fit$inputs$env_df, stringsAsFactors = FALSE)
    rownames(env_raw) = rownames(fit$inputs$env_df)
    X_env     = .dummy_matrix_from_info(env_raw, fit$model$env_pca_info,  fit$model$env_pca_vars)
    X_env_std = .standardize_like_training(X_env, fit$model$env_pca_center, fit$model$env_pca_scale)
    env_pc_df = .predict_pcs(fit$model$env_pca, X_env_std, pc_prefix = "ENV_")
    rownames(env_pc_df) = rownames(env_raw)
    env_df_for_model = env_pc_df
  } else {
    env_df_for_model = as.data.frame(fit$model$env_df_z)
  }

  # ---------- Prepare trait PCs for invaders ----------------------------------
  # Align raw factor levels to training levels if fixed-effects retained raw factors.
  mf     = try(fit$residents$fit_r$frame, silent = TRUE)
  tr_raw = as.data.frame(traits_inv, stringsAsFactors = FALSE)
  tr_raw = .coerce_factor_levels(tr_raw, if (!inherits(mf, "try-error")) mf else tr_raw)

  if (!is.null(fit$model$traits_pca)) {
    X_tr     = .dummy_matrix_from_info(tr_raw, fit$model$traits_pca_info,  fit$model$traits_pca_vars)
    X_tr_std = .standardize_like_training(X_tr, fit$model$traits_pca_center, fit$model$traits_pca_scale)
    tr_pc_df = .predict_pcs(fit$model$traits_pca, X_tr_std, pc_prefix = "TR_")
    rownames(tr_pc_df) = rownames(tr_raw)
    # keep raw factor terms that appear in the model frame
    keep_raw = if (!inherits(mf, "try-error")) intersect(names(tr_raw), names(mf)) else character(0)
    traits_for_model = cbind(tr_pc_df, tr_raw[, keep_raw, drop = FALSE])
  } else {
    # No trait PCA: standardise like residents (helper expected to be available in package)
    traits_for_model = .scale_like(as.data.frame(tr_raw),
                                   ref_means = fit$model$trait_means,
                                   ref_sds   = fit$model$trait_sds)
  }

  # ---------- Build invader predictors via the package builder -----------------
  sites           = fit$meta$sites
  env_for_builder = env_df_for_model[sites, , drop = FALSE]

  inv = build_invader_predictors(
    fit_r            = fit$residents$fit_r,
    env_df_model     = env_for_builder,
    traits_inv_model = traits_for_model,
    sites            = sites,
    inv_ids          = rownames(traits_inv),
    r_mu_s           = fit$residents$r_mu_s %||% stats::setNames(rep(0, length(sites)), sites),
    r_sd_s           = fit$residents$r_sd_s %||% stats::setNames(rep(1, length(sites)), sites),
    W_site           = fit$crowding$W_site,
    gower_all        = fit$traits$gower,
    res_ids          = fit$meta$residents,
    sigma_alpha      = fit$crowding$sigma_alpha,
    C_mu_s           = fit$residents$C_mu_s %||% stats::setNames(rep(0, length(sites)), sites),
    C_sd_s           = fit$residents$C_sd_s %||% stats::setNames(rep(1, length(sites)), sites),
    S_s_z            = fit$residents$S_s_z %||% stats::setNames(rep(0, length(sites)), sites),
    verbose          = TRUE
  )

  # ---------- Assemble outputs for downstream use ------------------------------
  traits_inv_raw  = as.data.frame(traits_inv)
  if (is.null(rownames(traits_inv_raw))) {
    # Fall back to IDs from the builder or tidy df (keeps names consistent)
    traits_inv_raw <- tibble::rownames_to_column(traits_inv_raw, "invader")
    rownames(traits_inv_raw) <- inv$inv_ids %||%
      rownames(inv$df |> dplyr::distinct(invader) |> as.data.frame())
  }
  traits_inv_glmm <- traits_for_model  # store the design-scale traits we actually used

  fit$invaders = list(
    traits_inv_raw  = traits_inv_raw,
    traits_inv_glmm = traits_inv_glmm,
    r_is_z = inv$r_is_z,
    C_is_z = inv$C_is_z,
    S_is_z = inv$S_is_z,
    df     = inv$df
  )

  fit
}
