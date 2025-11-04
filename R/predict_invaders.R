# ======================================================================
# 4) INVADER PREDICTORS (standardised): r_is_z, C_is_z, S_is_z
# ======================================================================

#' @title Build standardized invader predictors
#' @description Uses resident moments and the residents-only model to produce
#'   r_is_z, C_is_z, S_is_z. Handles PCA projection with dummy expansion
#'   to match the training design used in model_residents().
#' @param fit from learn_sensitivities()/prepare_trait_space() + model_residents()
#' @param traits_inv Invader trait table (raw scale; columns match resident trait names)
#' @return updated invasimapr_fit
#' @export
predict_invaders = function(fit, traits_inv) {

  # # build_invader_predictors
  # try({
  #   source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/build_invader_predictors.R")
  # }, silent = TRUE)
  # if (!exists("build_invader_predictors")) {
  #   stop("build_invader_predictors() not found.")
  # }

  stopifnot(inherits(fit, "invasimapr_fit"))

  `%||%` = function(a, b) if (!is.null(a)) a else b

  # ---------- helpers copied/paired with model_residents() ----------
  # Force factor levels to match training where those raw factor terms remain in the formula
  .coerce_factor_levels = function(df, model_frame) {
    # any columns in model frame that are factors and also present in df -> align levels
    mf = model_frame
    common = intersect(names(df), names(mf))
    for (v in common) {
      if (is.factor(mf[[v]])) {
        lv = levels(mf[[v]])
        x  = df[[v]]
        # coerce to character then factor with target levels (drop unseen)
        x = as.character(x)
        # redirect unseen to "_other_" if present
        if ("_other_" %in% lv) {
          x[!x %in% lv] = "_other_"
        }
        df[[v]] = factor(x, levels = lv, ordered = is.ordered(mf[[v]]))
      }
    }
    df
  }

  # Build a dummy matrix for a raw df using the training factor level map and
  # return an X with columns exactly equal to train_vars (order preserved).
  .dummy_matrix_from_info = function(df_raw, vars_info, train_vars) {
    df = as.data.frame(df_raw, stringsAsFactors = FALSE)
    rn = rownames(df)

    # Coerce characters to factors
    for (nm in names(df)) if (is.character(df[[nm]])) df[[nm]] = factor(df[[nm]])

    # numeric block
    num_idx = names(df)[vapply(df, is.numeric, logical(1))]
    X_num = if (length(num_idx)) as.matrix(df[, num_idx, drop = FALSE]) else
      matrix(numeric(0), nrow = nrow(df), ncol = 0, dimnames = list(rn, NULL))

    # dummy block from factor map saved at training
    # vars_info$factors is a named list: var -> levels used in training
    X_fac = NULL
    if (!is.null(vars_info$factors) && length(vars_info$factors)) {
      MM_list = lapply(names(vars_info$factors), function(v) {
        # Ensure the column exists; if missing, fabricate NA
        if (!v %in% names(df)) df[[v]] = NA
        x = as.character(df[[v]])
        lv = as.character(vars_info$factors[[v]])
        # redirect unseen to "_other_" if present; otherwise keep as-is (will become all zeros)
        if ("_other_" %in% lv) {
          x[!x %in% lv] = "_other_"
        }
        f = factor(x, levels = lv)
        mm = stats::model.matrix(~ . - 1, data = data.frame(f = f))
        colnames(mm) = paste0(v, lv)  # same naming as in model_residents()
        mm
      })
      X_fac = do.call(cbind, MM_list)
    } else {
      X_fac = matrix(numeric(0), nrow = nrow(df), ncol = 0, dimnames = list(rn, NULL))
    }

    # Bind numeric + dummy; drop any dup column names if any (shouldn't happen)
    X = if (ncol(X_num) && ncol(X_fac)) cbind(X_num, X_fac) else
      if (ncol(X_num)) X_num else X_fac

    # Now conform to train_vars: add missing columns (all zeros), reorder, drop extras
    if (length(train_vars)) {
      miss = setdiff(train_vars, colnames(X))
      if (length(miss)) {
        add = matrix(0, nrow = nrow(X), ncol = length(miss),
                      dimnames = list(rownames(X), miss))
        X = if (ncol(X)) cbind(X, add) else add
      }
      X = X[, train_vars, drop = FALSE]
    } else {
      # no training vars → empty matrix
      X = matrix(numeric(0), nrow = nrow(df), ncol = 0, dimnames = list(rn, NULL))
    }

    X
  }

  # Standardize columns of X using training center/scale (vectors named by train_vars)
  .standardize_like_training = function(X, center, scale) {
    if (!length(center) && !length(scale)) return(X)
    # fill any NA columns with training center to get 0 after standardization
    for (j in seq_along(center)) {
      nm = names(center)[j]
      if (!nm %in% colnames(X)) next
      x = X[, nm]
      x[!is.finite(x)] = center[j]
      X[, nm] = (x - center[j]) / (if (is.finite(scale[j]) && scale[j] != 0) scale[j] else 1)
    }
    X
  }

  # Project into PCs using prcomp object; assumes X_std has columns == rownames(rotation)
  .predict_pcs = function(pca, X_std, pc_prefix) {
    if (is.null(pca) || !inherits(pca, "prcomp")) {
      return(data.frame(row.names = rownames(X_std)))
    }
    rot_names = rownames(pca$rotation)
    if (is.null(rot_names)) {
      stop("Stored PCA rotation has no rownames; cannot project.")
    }
    # Align columns to rotation rownames
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

  # ---------- prepare PCA frames for env & traits ----------
  # ENV PCs for all sites
  if (!is.null(fit$model$env_pca)) {
    # raw env df used for training (rows = sites)
    env_raw = as.data.frame(fit$inputs$env_df, stringsAsFactors = FALSE)
    rownames(env_raw) = rownames(fit$inputs$env_df)
    # Build dummy matrix with training factor map & vars, then standardize and project
    X_env = .dummy_matrix_from_info(env_raw, fit$model$env_pca_info, fit$model$env_pca_vars)
    X_env_std = .standardize_like_training(X_env, fit$model$env_pca_center, fit$model$env_pca_scale)
    env_pc_df = .predict_pcs(fit$model$env_pca, X_env_std, pc_prefix = "ENV_")
    rownames(env_pc_df) = rownames(env_raw)
    env_df_for_model = env_pc_df
  } else {
    env_df_for_model = as.data.frame(fit$model$env_df_z)
  }

  # TRAIT PCs for invaders + keep raw factor trait columns for the fixed effects
  # Start from user-supplied traits_inv (raw); align raw factor columns to model frame levels
  mf = try(fit$residents$fit_r$frame, silent = TRUE)
  tr_raw = as.data.frame(traits_inv, stringsAsFactors = FALSE)
  tr_raw = .coerce_factor_levels(tr_raw, if (!inherits(mf, "try-error")) mf else tr_raw)

  if (!is.null(fit$model$traits_pca)) {
    X_tr = .dummy_matrix_from_info(tr_raw, fit$model$traits_pca_info, fit$model$traits_pca_vars)
    X_tr_std = .standardize_like_training(X_tr, fit$model$traits_pca_center, fit$model$traits_pca_scale)
    tr_pc_df = .predict_pcs(fit$model$traits_pca, X_tr_std, pc_prefix = "TR_")
    rownames(tr_pc_df) = rownames(tr_raw)
    traits_for_model = cbind(tr_pc_df, tr_raw[, intersect(names(tr_raw), names(mf)), drop = FALSE])
  } else {
    # no PCA: just standardize like residents
    traits_for_model = .scale_like(as.data.frame(tr_raw),
                                    ref_means = fit$model$trait_means,
                                    ref_sds   = fit$model$trait_sds)
  }

  # ---------- call the builder ----------
  sites = fit$meta$sites
  env_for_builder = env_df_for_model[sites, , drop = FALSE]

  inv = build_invader_predictors(
    fit_r            = fit$residents$fit_r,
    env_df_model     = env_for_builder,   # renamed arg (see function below)
    traits_inv_model = traits_for_model,  # renamed arg (see function below)
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

  # inside predict_invaders() right before returning ----

  # keep the raw (pre-PCA / pre-standardisation) traits for summaries
  traits_inv_raw = as.data.frame(traits_inv)
  if (is.null(rownames(traits_inv_raw))) {
    rownames(traits_inv_raw) = inv$inv_ids %||% rownames(inv$df |> dplyr::distinct(invader) |> as.data.frame())
  }

  fit$invaders = list(
    traits_inv_raw  = traits_inv_raw,                 # =- NEW
    traits_inv_glmm = traits_inv_glmm,
    r_is_z = inv$r_is_z,
    C_is_z = inv$C_is_z,
    S_is_z = inv$S_is_z,
    df     = inv$df
  )
  fit

}
