#' Model residents and build standardized predictors (\eqn{r_{js}^{(z)}, C_{js}^{(z)}, S_{js}^{(z)}})
#'
#' Fits the **residents-only** GLMM and constructs site-standardised predictors
#' used to learn sensitivities and to project invaders. Specifically:
#' (i) builds a design using environment and resident traits,
#' (ii) optionally **dummy-expands** factors and performs **PCA compression**
#' with stored centering/scale metadata for reproducible projection,
#' (iii) fits a `glmmTMB` model on resident abundance,
#' (iv) predicts linear predictors \eqn{r_{js}} then applies **row-wise
#' z-scoring** to obtain \eqn{r_{js}^{(z)}}, and (v) row-z-standardises
#' **resident crowding** \eqn{C_{js}} and computes site saturation summaries
#' \eqn{S_s}, \eqn{S_{s}^{(z)}}, \eqn{S_{js}^{(z)}}.
#'
#' @section Reduction strategies:
#' Use `reduce_strategy` to control the fixed-effects size/complexity.
#' - `"auto"`: preflights a dense design estimate; if memory budget exceeded,
#'   falls back to `"no_interactions"`, then `"pca"`, then `"pca+no_interactions"`.
#' - `"no_interactions"`: drops all env × trait interactions.
#' - `"pca"`: dummy-expand → standardise → `prcomp()` on the design blocks
#'   (environment and traits) and use the first `pca_env_k`/`pca_trait_k` PCs.
#' - `"none"`: keep full model as requested by `include_env_trait_interactions`.
#'
#' @param fit An `invasimapr_fit` returned by
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/prepare_trait_space.html}{`prepare_trait_space()`};
#'   must contain `inputs$comm_res`, `inputs$env_df`, `inputs$traits_res`,
#'   and `crowding$C_js` (raw, not z-scored).
#' @param family A `glmmTMB` family (default
#'   `glmmTMB::tweedie(link = "log")`).
#' @param include_env_trait_interactions Logical; whether to include
#'   environment × trait interactions in the fixed effects (when not reduced).
#' @param saturation_mode One of
#'   `c("evenness_deficit","opportunity_penalty","modelled_dominance")`;
#'   forwarded to
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_site_saturation.html}{`compute_site_saturation()`}.
#' @param robust_r,robust_c Logical; use robust row-wise z-scoring for
#'   \eqn{r_{js}} and \eqn{C_{js}} respectively.
#' @param fit_model Logical; if `FALSE`, preflight only (build design/estimate
#'   memory footprint) and return without fitting.
#' @param max_dense_gb Numeric; maximum allowed dense fixed-effects design size
#'   (in GB) for preflight in `"auto"` mode.
#' @param reduce_strategy One of `c("auto","none","no_interactions","pca")`.
#' @param pca_env_k,pca_trait_k Integers; number of PCs to retain for the environment
#'   and trait blocks under PCA reduction.
#' @param verbose Logical; emit progress/messages.
#'
#' @return The input `invasimapr_fit` with updated components:
#' \describe{
#'   \item{`model`}{Lists of standardised inputs (`env_df_z`, `traits_res_glmm`),
#'     means/SDs, and—if PCA used—`*_pca`, `*_pca_vars`, `*_pca_center`,
#'     `*_pca_scale`, and `*_pca_info` needed for invader projection.}
#'   \item{`residents`}{GLMM fit (`fit_r`), model frame (`dat_r`), grid for predictions,
#'     raw linear predictor matrix `r_js`, mean scale `mu_js`, and row-z outputs
#'     `r_js_z` with per-site means/SDs; similarly `C_js_z` and saturation summaries
#'     (`S_s`, `S_s_z`, `S_js_z`).}
#'   \item{`model$pca_used`}{Flags and counts recording whether PCA was used and
#'     the retained dimensionality.}
#' }
#'
#' @details
#' **Design construction**
#' The base frames `env_df_z` and `traits_res_glmm` are created by robust
#' standardisation. Characters are coerced to factors to ensure stable
#' dummy-expansion. In `"pca"` reduction, each block undergoes:
#' `model.matrix` one-hot encoding (no intercept) → per-column standardisation →
#' `prcomp()` with **stored** centring/scales and rotation rownames. This metadata
#' supports lossless, future projection of invaders (see
#' \href{https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.html}{`predict_invaders()`}).
#'
#' **Preflight memory check**
#' A rough dense design GB estimate is computed on the fixed-effects portion
#' (random terms stripped). `"auto"` reduction iteratively lowers complexity to
#' meet `max_dense_gb`.
#'
#' **Z-standardisation**
#' Row-wise z-scoring uses `.row_z()` (robust if requested) and stores site-wise
#' means/SDs for later use in mapping probabilities and fitness.
#'
#' @seealso
#' - Formula builder:
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/build_model_formula.html}{`build_model_formula()`}
#' - GLMM preparation:
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/prep_resident_glmm.html}{`prep_resident_glmm()`}
#' - Site saturation:
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/compute_site_saturation.html}{`compute_site_saturation()`}
#' - Trait preparation:
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/prepare_trait_space.html}{`prepare_trait_space()`}
#' - Invader projection:
#'   \href{https://b-cubed-eu.github.io/invasimapr/reference/predict_invaders.html}{`predict_invaders()`}
#'
#' @examples
#' \dontrun{
#' fit <- prepare_inputs(sites = site_df, residents = resident_df,
#'                       invaders = invader_df, traits = trait_df)
#' fit <- prepare_trait_space(fit, traits_inv)
#' fit <- model_residents(fit, reduce_strategy = "auto", verbose = TRUE)
#'
#' # Inspect z-standardised predictors for residents
#' dim(fit$residents$r_js_z); dim(fit$residents$C_js_z)
#' fit$residents$messages
#' }
#'
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

  # --- Preconditions ------------------------------------------------------------
  stopifnot(inherits(fit, "invasimapr_fit"))
  saturation_mode = match.arg(saturation_mode)
  reduce_strategy = match.arg(reduce_strategy)

  .std = function(df) .standardise_df(as.data.frame(df))  # robust standardiser (package-internal)
  .msg = function(...) if (isTRUE(verbose)) message(sprintf(...))

  # --- Standardise base model inputs (stored for back-compat) ------------------
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

  # Ensure character → factor before any dummy-expansion
  for (nm in names(fit$model$env_df_z)) {
    if (is.character(fit$model$env_df_z[[nm]])) fit$model$env_df_z[[nm]] = factor(fit$model$env_df_z[[nm]])
  }
  for (nm in names(fit$model$traits_res_glmm)) {
    if (is.character(fit$model$traits_res_glmm[[nm]])) fit$model$traits_res_glmm[[nm]] = factor(fit$model$traits_res_glmm[[nm]])
  }

  # --- Helpers for preflight sizing --------------------------------------------
  strip_random = function(fml) {
    # Remove random terms (..|..) /(..||..) to estimate fixed-effects size only
    f = paste0(deparse(fml), collapse = "")
    rhs = sub("^[^~]+~", "", f)
    rhs = gsub("\\+?\\s*\\([^()]*\\|[^()]*\\)", "", rhs)
    rhs = gsub("\\s+\\+\\s+$", "", rhs)
    stats::as.formula(paste0("abundance ~ ", rhs))
  }
  estimate_dense_gb = function(fml_fixed, dat) {
    trm = stats::terms(fml_fixed, data = dat)
    X   = Matrix::sparse.model.matrix(trm, dat)
    as.numeric(nrow(X)) * as.numeric(ncol(X)) * 8 / 1e9
  }

  # --- Dummy-expand + PCA on standardised matrix (block-wise) ------------------
  # Creates reproducible PCA objects with stored center/scale and training vars.
  dummy_pca = function(df, k, pc_prefix) {
    df = as.data.frame(df, stringsAsFactors = FALSE)
    rn = rownames(df)
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

    # full one-hot dummies (no intercept) for factor columns
    if (any(fac_idx)) {
      fac_names = names(df)[fac_idx]
      vars_info$factors = lapply(fac_names, function(v) levels(df[[v]]))
      names(vars_info$factors) = fac_names

      MM_list = lapply(fac_names, function(v) {
        mm = stats::model.matrix(~ . - 1, data = df[v],
                                 contrasts.arg = stats::contrasts(df[[v]], contrasts = FALSE))
        colnames(mm) = paste0(v, levels(df[[v]]))
        mm
      })
      X_fac = do.call(cbind, MM_list)
      X = if (ncol(X) == 0L) X_fac else cbind(X, X_fac)
      vars_info$dummy = colnames(X_fac)
    } else {
      vars_info$factors = list(); vars_info$dummy = character(0)
    }

    # drop zero-variance columns (constant dummies/numerics)
    if (ncol(X) > 0L) {
      sds = apply(X, 2, stats::sd, na.rm = TRUE)
      keep = is.finite(sds) & sds > 0
      if (any(!keep)) X = X[, keep, drop = FALSE]
    }

    # per-column standardisation (robust center) for PCA
    if (ncol(X) > 0L) {
      center = apply(X, 2, stats::median, na.rm = TRUE)
      scalev = apply(X, 2, stats::sd,     na.rm = TRUE)
      scalev[!is.finite(scalev) | scalev == 0] = 1
      Xs = scale(X, center = center, scale = scalev)
    } else {
      center = numeric(0); scalev = numeric(0); Xs = X
    }

    if (ncol(Xs) > 0L) {
      pcs_obj = stats::prcomp(Xs, center = FALSE, scale. = FALSE)
      kk = min(k, ncol(pcs_obj$x))
      pcs_df = as.data.frame(pcs_obj$x[, seq_len(kk), drop = FALSE])
      colnames(pcs_df) = paste0(pc_prefix, "PC", seq_len(kk))
      rownames(pcs_df) = rn
    } else {
      pcs_obj = NULL; pcs_df = data.frame(row.names = rn)
    }

    list(
      pcs_df     = pcs_df,
      pca        = pcs_obj,
      train_vars = colnames(Xs),
      center     = center,
      scale      = scalev,
      vars_info  = vars_info
    )
  }

  # --- Choose frames & interactions per reduction strategy ---------------------
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

    fit$model$env_pca          = envc$pca
    fit$model$traits_pca       = trc$pca
    fit$model$env_pca_vars     = envc$train_vars
    fit$model$traits_pca_vars  = trc$train_vars
    fit$model$env_pca_center   = envc$center
    fit$model$env_pca_scale    = envc$scale
    fit$model$traits_pca_center= trc$center
    fit$model$traits_pca_scale = trc$scale
    fit$model$env_pca_info     = envc$vars_info
    fit$model$traits_pca_info  = trc$vars_info

    path_used = "pca"
  }

  # --- Build initial formula ---------------------------------------------------
  fml0 = build_model_formula(
    response = "abundance",
    env_df   = env_use,
    trait_df = traits_use,
    include_env_trait_interactions = incl_int,
    random_intercepts = c("site","species"),
    random_slopes     = NULL,
    backend = "glmmTMB"
  )

  # --- Preflight (without fitting) --------------------------------------------
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

  # Automatic fallback path when dense design is too large
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
      attempt   = attempt2; fml0 = f_no; path_used = "auto->no_interactions"
    } else {
      .msg("Still large. Trying PCA compression …")
      envc = dummy_pca(fit$model$env_df_z,        pca_env_k,   pc_prefix = "ENV_")
      trc  = dummy_pca(fit$model$traits_res_glmm, pca_trait_k, pc_prefix = "TR_")
      env_use2    = envc$pcs_df
      traits_use2 = trc$pcs_df

      fit$model$env_pca           = envc$pca
      fit$model$traits_pca        = trc$pca
      fit$model$env_pca_vars      = envc$train_vars
      fit$model$traits_pca_vars   = trc$train_vars
      fit$model$env_pca_center    = envc$center
      fit$model$env_pca_scale     = envc$scale
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
        attempt   = attempt3; fml0 = f_pca
        env_use   = env_use2; traits_use = traits_use2
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
          attempt   = attempt4; fml0 = f_pca_no
          env_use   = env_use2; traits_use = traits_use2
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

  # --- Early return if only preflighting ---------------------------------------
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

  # --- Fit residents GLMM ------------------------------------------------------
  fit_r = glmmTMB::glmmTMB(formula = fml0, data = attempt$rg$dat_r, family = family)

  # --- Predict resident linear predictor on full site × species grid -----------
  sites   = levels(attempt$rg$dat_r$site)
  res_ids = levels(attempt$rg$dat_r$species)

  grid_res = dplyr::distinct(attempt$rg$dat_r[, c("site","species")])
  grid_res = dplyr::left_join(grid_res, tibble::rownames_to_column(env_use,    "site"),    by = "site")
  grid_res = dplyr::left_join(grid_res, tibble::rownames_to_column(traits_use, "species"), by = "species")

  grid_res$eta_r = stats::predict(
    fit_r, newdata = grid_res, type = "link",
    re.form = NA, allow.new.levels = TRUE
  )

  r_js = matrix(NA_real_, nrow = length(sites), ncol = length(res_ids),
                dimnames = list(sites, res_ids))
  ii = match(grid_res$site, sites); jj = match(grid_res$species, res_ids)
  r_js[cbind(ii, jj)] = grid_res$eta_r
  mu_js = exp(r_js)

  # --- Row-wise z-standardisations & saturation summaries ----------------------
  r_std = .row_z(r_js,              robust = robust_r)             # r_js → r_js_z
  C_std = .row_z(fit$crowding$C_js, robust = robust_c)             # C_js → C_js_z

  sat = compute_site_saturation(
    mode     = saturation_mode,
    comm_res = fit$inputs$comm_res,
    res_ids  = colnames(fit$inputs$comm_res),
    mu_js    = if (identical(saturation_mode, "modelled_dominance")) mu_js else NULL
  )

  # --- Pack resident outputs ---------------------------------------------------
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

  # Record PCA usage for downstream invader projection
  fit$model$pca_used = list(
    env      = !is.null(fit$model$env_pca),
    traits   = !is.null(fit$model$traits_pca),
    k_env    = if (!is.null(fit$model$env_pca)) pca_env_k else 0,
    k_traits = if (!is.null(fit$model$traits_pca)) pca_trait_k else 0
  )

  fit
}
