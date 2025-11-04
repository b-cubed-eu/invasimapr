# ======================================================================
# 2) RESIDENT BACKBONE (GLMM) & STANDARDIZED PREDICTORS
# ======================================================================

#' @title Model residents and build standardized predictors (r_js_z, C_js_z, S_js_z)
#' @description
#' Builds a GLMM formula, fits the residents-only model, and constructs
#' standardized matrices. **Per-site z** for both r_js and C_js happens here.
#' If `reduce_strategy = "pca"`, this function expands factors to **dummies**,
#' standardizes all predictors, and runs PCA on the standardized matrix. The
#' PCA objects plus dummy/centering metadata are stored for robust projection.
#' @param fit from prepare_trait_space() (must contain inputs + crowding$C_js (raw))
#' @param family glmmTMB family
#' @param include_env_trait_interactions logical; for build_model_formula()
#' @param saturation_mode one of "evenness_deficit","opportunity_penalty","modelled_dominance"
#' @param robust_r logical; robust row-wise z for r_js
#' @param robust_c logical; robust row-wise z for C_js
#' @param fit_model logical; if FALSE, builds/preflights but does not fit
#' @param max_dense_gb numeric; threshold (GB) for dense design size preflight
#' @param reduce_strategy one of c("auto","none","no_interactions","pca")
#' @param pca_env_k, pca_trait_k integers; number of PCs if PCA reduction is used
#' @param verbose logical; emit messages
#' @return updated invasimapr_fit
#' @import data.table
#' @export
model_residents = function(fit,
                            family = glmmTMB::tweedie(link = "log"),
                            include_env_trait_interactions = TRUE,
                            saturation_mode = c("evenness_deficit",
                                                "opportunity_penalty",
                                                "modelled_dominance"),
                            robust_r = TRUE,
                            robust_c = TRUE,
                            fit_model = TRUE,
                            # preflight knobs:
                            max_dense_gb    = 8,
                            reduce_strategy = c("auto","none","no_interactions","pca"),
                            pca_env_k       = 5,
                            pca_trait_k     = 5,
                            verbose = TRUE) {

  # # bring helpers
  # try({ source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/build_model_formula.R") }, silent = TRUE)
  # if (!exists("build_model_formula")) stop("build_model_formula() not found.")
  # try({ source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/prep_resident_glmm.R") }, silent = TRUE)
  # if (!exists("prep_resident_glmm")) stop("prep_resident_glmm() not found.")
  # try({ source("D:/Methods/R/myR_Packages/b-cubed-versions/invasimapr/R/compute_site_saturation.R") }, silent = TRUE)
  # if (!exists("compute_site_saturation")) stop("compute_site_saturation() not found.")

  stopifnot(inherits(fit, "invasimapr_fit"))
  saturation_mode = match.arg(saturation_mode)
  reduce_strategy = match.arg(reduce_strategy)

  .std = function(df) .standardise_df(as.data.frame(df))
  .msg = function(...) if (isTRUE(verbose)) message(sprintf(...))

  # ---------- Standardize model inputs (base scale, kept for back-compat) ----------
  env_std = .std(fit$inputs$env_df)
  tr_std  = .std(fit$inputs$traits_res)

  fit$model = list(
    env_df_z        = env_std$data,
    env_means       = env_std$means,
    env_sds         = env_std$sds,
    traits_res_glmm = tr_std$data,
    trait_means     = tr_std$means,
    trait_sds       = tr_std$sds
  )

  # Ensure character → factor (so dummy-expansion is stable)
  for (nm in names(fit$model$env_df_z)) {
    if (is.character(fit$model$env_df_z[[nm]])) fit$model$env_df_z[[nm]] = factor(fit$model$env_df_z[[nm]])
  }
  for (nm in names(fit$model$traits_res_glmm)) {
    if (is.character(fit$model$traits_res_glmm[[nm]])) fit$model$traits_res_glmm[[nm]] = factor(fit$model$traits_res_glmm[[nm]])
  }

  # ---------- helpers for preflight sizing ----------
  strip_random = function(fml) {
    f = paste0(deparse(fml), collapse = "")
    rhs = sub("^[^~]+~", "", f)
    rhs = gsub("\\+?\\s*\\([^()]*\\|[^()]*\\)", "", rhs)  # remove (..|..), (..||..)
    rhs = gsub("\\s+\\+\\s+$", "", rhs)
    stats::as.formula(paste0("abundance ~ ", rhs))
  }
  estimate_dense_gb = function(fml_fixed, dat) {
    trm = stats::terms(fml_fixed, data = dat)
    X   = Matrix::sparse.model.matrix(trm, dat)
    as.numeric(nrow(X)) * as.numeric(ncol(X)) * 8 / 1e9
  }

  # ---------- NEW: dummy-expand + PCA on standardized matrix ----------
  # Build a numeric matrix from mixed df:
  #  - numeric kept as-is
  #  - factor/ordered/character -> one-hot dummies (no intercept), full rank
  # Then standardize each column (store center/scale) and run prcomp on the standardized X.
  dummy_pca = function(df, k, pc_prefix) {
    df = as.data.frame(df, stringsAsFactors = FALSE)
    rn = rownames(df)

    # Coerce character -> factor early
    for (nm in names(df)) if (is.character(df[[nm]])) df[[nm]] = factor(df[[nm]])

    num_idx = vapply(df, is.numeric, logical(1))
    fac_idx = vapply(df, is.factor,  logical(1)) | vapply(df, is.ordered, logical(1))

    X = NULL
    vars_info = list()

    # numeric block
    if (any(num_idx)) {
      X_num = as.matrix(df[, num_idx, drop = FALSE])
      colnames(X_num) = make.unique(colnames(X_num))
      X = X_num
      vars_info$numeric = colnames(X_num)
    } else {
      X = matrix(numeric(0), nrow = nrow(df), ncol = 0, dimnames = list(rn, NULL))
      vars_info$numeric = character(0)
    }

    # dummy block for each factor column (full one-hot, no intercept)
    if (any(fac_idx)) {
      fac_names = names(df)[fac_idx]
      vars_info$factors = lapply(fac_names, function(v) levels(df[[v]]))
      names(vars_info$factors) = fac_names

      # Build one model.matrix per factor to preserve clean column names
      MM_list = lapply(fac_names, function(v) {
        # use no intercept; this gives one column per level
        mm = stats::model.matrix(~ . - 1, data = df[v], contrasts.arg = stats::contrasts(df[[v]], contrasts = FALSE))
        # prefix with variable name to avoid collisions: v=habitat -> habitatGrassland, etc.
        colnames(mm) = paste0(v, levels(df[[v]]))
        mm
      })

      X_fac = do.call(cbind, MM_list)
      if (ncol(X_fac) > 0L) {
        # Bind preserving row order
        if (ncol(X) == 0L) {
          X = X_fac
        } else {
          X = cbind(X, X_fac)
        }
        vars_info$dummy = colnames(X_fac)
      } else {
        vars_info$dummy = character(0)
      }
    } else {
      vars_info$factors = list()
      vars_info$dummy   = character(0)
    }

    # Drop zero-variance columns (constant dummy or numeric)
    if (ncol(X) > 0L) {
      sds = apply(X, 2, stats::sd, na.rm = TRUE)
      keep = is.finite(sds) & sds > 0
      if (any(!keep)) {
        X = X[, keep, drop = FALSE]
      }
    }

    # Standardize (store center/scale used for PCA)
    if (ncol(X) > 0L) {
      center = apply(X, 2, stats::median, na.rm = TRUE)  # robust-ish center
      scalev = apply(X, 2, stats::sd,     na.rm = TRUE)
      scalev[!is.finite(scalev) | scalev == 0] = 1
      Xs = scale(X, center = center, scale = scalev)
    } else {
      center = numeric(0); scalev = numeric(0)
      Xs = X
    }

    # PCA on standardized matrix (no extra centering/scaling in prcomp)
    if (ncol(Xs) > 0L) {
      pcs_obj = stats::prcomp(Xs, center = FALSE, scale. = FALSE)
      kk = min(k, ncol(pcs_obj$x))
      pcs_df = as.data.frame(pcs_obj$x[, seq_len(kk), drop = FALSE])
      colnames(pcs_df) = paste0(pc_prefix, "PC", seq_len(kk))
      rownames(pcs_df) = rn
      # ensure rotation rownames record the training variables (after dummy+std)
      # prcomp sets rownames(rotation) to colnames(Xs) automatically
    } else {
      pcs_obj = NULL
      pcs_df  = data.frame(row.names = rn)
    }

    list(
      pcs_df   = pcs_df,
      pca      = pcs_obj,
      train_vars = colnames(Xs),
      center   = center,
      scale    = scalev,
      vars_info = vars_info
    )
  }

  # ---------- choose frames & interactions per reduce_strategy ----------
  env_use    = fit$model$env_df_z
  traits_use = fit$model$traits_res_glmm
  incl_int   = include_env_trait_interactions
  path_used  = "original"

  if (reduce_strategy == "no_interactions") {
    incl_int  = FALSE
    path_used = "no_interactions"
  } else if (reduce_strategy == "pca") {
    envc = dummy_pca(env_use,    pca_env_k,   pc_prefix = "ENV_")
    trc  = dummy_pca(traits_use, pca_trait_k, pc_prefix = "TR_")

    env_use    = envc$pcs_df
    traits_use = trc$pcs_df

    fit$model$env_pca       = envc$pca
    fit$model$traits_pca    = trc$pca
    fit$model$env_pca_vars  = envc$train_vars         # dummy+numeric names used in PCA
    fit$model$traits_pca_vars = trc$train_vars
    fit$model$env_pca_center  = envc$center
    fit$model$env_pca_scale   = envc$scale
    fit$model$traits_pca_center = trc$center
    fit$model$traits_pca_scale  = trc$scale
    fit$model$env_pca_info      = envc$vars_info       # factor levels per original var
    fit$model$traits_pca_info   = trc$vars_info

    path_used = "pca"
  }

  # ---------- build initial formula ----------
  fml0 = build_model_formula(
    response = "abundance",
    env_df   = env_use,
    trait_df = traits_use,
    include_env_trait_interactions = incl_int,
    random_intercepts = c("site","species"),
    random_slopes     = NULL,
    backend = "glmmTMB"
  )

  # ---------- preflight (without fitting) ----------
  preflight_once = function(fml, env_df, trait_df) {
    rg = prep_resident_glmm(
      comm_res         = fit$inputs$comm_res,
      env_df_z         = env_df,
      traits_res_glmm  = trait_df,
      fml              = fml,
      family           = family,
      fit_model        = FALSE
    )
    gb = tryCatch(estimate_dense_gb(strip_random(fml), rg$dat_r), error = function(e) Inf)
    list(rg = rg, gb = gb)
  }

  attempt = preflight_once(fml0, env_use, traits_use)
  .msg("Design estimate (%s): %.2f GB", path_used, attempt$gb)

  if (reduce_strategy == "auto" && is.finite(attempt$gb) && attempt$gb > max_dense_gb) {
    .msg("Too large (%.2f GB > %.2f GB). Trying no_interactions …", attempt$gb, max_dense_gb)
    f_no = build_model_formula(
      response = "abundance", env_df = env_use, trait_df = traits_use,
      include_env_trait_interactions = FALSE,
      random_intercepts = c("site","species"),
      random_slopes     = NULL, backend = "glmmTMB"
    )
    attempt2 = preflight_once(f_no, env_use, traits_use)
    .msg("Design estimate (no_interactions): %.2f GB", attempt2$gb)

    if (is.finite(attempt2$gb) && attempt2$gb <= max_dense_gb) {
      attempt   = attempt2
      fml0      = f_no
      path_used = "auto->no_interactions"
    } else {
      .msg("Still large. Trying PCA compression …")
      # Run PCA with dummy-expansion here
      envc = dummy_pca(fit$model$env_df_z,        pca_env_k,   pc_prefix = "ENV_")
      trc  = dummy_pca(fit$model$traits_res_glmm, pca_trait_k, pc_prefix = "TR_")
      env_use2    = envc$pcs_df
      traits_use2 = trc$pcs_df

      fit$model$env_pca         = envc$pca
      fit$model$traits_pca      = trc$pca
      fit$model$env_pca_vars    = envc$train_vars
      fit$model$traits_pca_vars = trc$train_vars
      fit$model$env_pca_center  = envc$center
      fit$model$env_pca_scale   = envc$scale
      fit$model$traits_pca_center = trc$center
      fit$model$traits_pca_scale  = trc$scale
      fit$model$env_pca_info      = envc$vars_info
      fit$model$traits_pca_info   = trc$vars_info

      f_pca = build_model_formula(
        response = "abundance",
        env_df   = env_use2,
        trait_df = traits_use2,
        include_env_trait_interactions = include_env_trait_interactions,
        random_intercepts = c("site","species"),
        random_slopes     = NULL,
        backend = "glmmTMB"
      )
      attempt3 = preflight_once(f_pca, env_use2, traits_use2)
      .msg("Design estimate (pca): %.2f GB", attempt3$gb)

      if (is.finite(attempt3$gb) && attempt3$gb <= max_dense_gb) {
        attempt   = attempt3
        fml0      = f_pca
        env_use   = env_use2
        traits_use= traits_use2
        path_used = "auto->pca"
      } else {
        .msg("Still large. Forcing PCA + no_interactions …")
        f_pca_no = build_model_formula(
          response = "abundance",
          env_df   = env_use2,
          trait_df = traits_use2,
          include_env_trait_interactions = FALSE,
          random_intercepts = c("site","species"),
          random_slopes     = NULL,
          backend = "glmmTMB"
        )
        attempt4 = preflight_once(f_pca_no, env_use2, traits_use2)
        .msg("Design estimate (pca+no_interactions): %.2f GB", attempt4$gb)

        if (is.finite(attempt4$gb) && attempt4$gb <= max_dense_gb) {
          attempt   = attempt4
          fml0      = f_pca_no
          env_use   = env_use2
          traits_use= traits_use2
          path_used = "auto->pca+no_interactions"
        } else {
          stop(sprintf(
            "Fixed-effects design remains too large (%.2f GB > %.2f GB) even after reductions. ",
            attempt4$gb, max_dense_gb
          ), "Reduce predictors further or partition the data.")
        }
      }
    }
  }

  # ---------- fit ----------
  if (!isTRUE(fit_model)) {
    fit$residents = list(
      fit_r    = NULL,
      dat_r    = attempt$rg$dat_r,
      grid_res = NULL,
      fml      = fml0,
      r_js     = NULL,
      mu_js    = NULL,
      messages = c(sprintf("preflight_gb=%.2f", attempt$gb), sprintf("path=%s", path_used))
    )
    return(fit)
  }

  fit_r = glmmTMB::glmmTMB(formula = fml0, data = attempt$rg$dat_r, family = family)

  # ---------- predictions ----------
  sites   = levels(attempt$rg$dat_r$site)
  res_ids = levels(attempt$rg$dat_r$species)

  grid_res = dplyr::distinct(attempt$rg$dat_r[, c("site","species")])
  grid_res = dplyr::left_join(grid_res, tibble::rownames_to_column(env_use,    "site"),    by = "site")
  grid_res = dplyr::left_join(grid_res, tibble::rownames_to_column(traits_use, "species"), by = "species")

  grid_res$eta_r = stats::predict(fit_r, newdata = grid_res, type = "link",
                                   re.form = NA, allow.new.levels = TRUE)

  r_js = matrix(NA_real_, nrow = length(sites), ncol = length(res_ids),
                 dimnames = list(sites, res_ids))
  ii = match(grid_res$site, sites); jj = match(grid_res$species, res_ids)
  r_js[cbind(ii, jj)] = grid_res$eta_r
  mu_js = exp(r_js)

  # ---------- z-standardisations & saturation ----------
  r_std = .row_z(r_js,              robust = robust_r)
  C_std = .row_z(fit$crowding$C_js, robust = robust_c)

  sat = compute_site_saturation(
    mode     = saturation_mode,
    comm_res = fit$inputs$comm_res,
    res_ids  = colnames(fit$inputs$comm_res),
    mu_js    = if (identical(saturation_mode, "modelled_dominance")) mu_js else NULL
  )

  fit$residents = list(
    fit_r    = fit_r,
    dat_r    = attempt$rg$dat_r,
    grid_res = grid_res,
    fml      = fml0,
    r_js     = r_js,
    mu_js    = mu_js,
    r_js_z   = r_std$z,
    r_mu_s   = r_std$mu,
    r_sd_s   = r_std$sd,
    C_js_z   = C_std$z,
    C_mu_s   = C_std$mu,
    C_sd_s   = C_std$sd,
    S_s      = sat$S_s,
    S_s_z    = sat$S_s_z,
    S_js_z   = sat$S_js_z,
    messages = c(sprintf("preflight_gb=%.2f", attempt$gb), sprintf("path=%s", path_used))
  )

  fit$model$pca_used = list(
    env      = !is.null(fit$model$env_pca),
    traits   = !is.null(fit$model$traits_pca),
    k_env    = if (!is.null(fit$model$env_pca)) pca_env_k else 0,
    k_traits = if (!is.null(fit$model$traits_pca)) pca_trait_k else 0
  )

  fit
}
