#' Build standardized invader predictors r_is_z, C_is_z, S_is_z
#' (expects env/traits frames already projected to model columns)
#'
#' @param fit_r Fitted residents-only glmmTMB model.
#' @param env_df_model Data frame for environment predictors used by the model
#'        (rows = sites). This is either ENV_PC* (if PCA used) or env_df_z.
#' @param traits_inv_model Data frame for invader traits used by the model
#'        (rows = invaders). This includes TR_PC* (if PCA used) and any raw
#'        factor traits present in the model formula.
#' @param ... (unchanged other args)
#' @export
build_invader_predictors = function(fit_r, env_df_model, traits_inv_model,
                                     sites, inv_ids,
                                     r_mu_s, r_sd_s,
                                     W_site, gower_all, res_ids, sigma_alpha,
                                     C_mu_s, C_sd_s, S_s_z,
                                     verbose = TRUE) {

  ## ---------- Shapes & basic checks ----------
  W_site   = if (is.matrix(W_site)) W_site else as.matrix(W_site)
  stopifnot(identical(rownames(W_site), sites),
            identical(colnames(W_site), res_ids))

  # Ensure rownames
  if (is.null(rownames(env_df_model))) stop("env_df_model must have rownames = sites.")
  if (is.null(rownames(traits_inv_model))) stop("traits_inv_model must have rownames = invader IDs.")

  ## ---------- Gower -> kernel -> C_is (S×I) ----------
  D_all = if (inherits(gower_all, "dist")) as.matrix(gower_all) else as.matrix(gower_all)
  if (is.null(colnames(D_all)) && !is.null(rownames(D_all))) colnames(D_all) = rownames(D_all)

  miss_r = setdiff(res_ids, rownames(D_all))
  miss_i = setdiff(inv_ids, rownames(D_all))
  if (length(miss_r) || length(miss_i)) {
    stop("gower_all is missing IDs: ",
         if (length(miss_r)) paste0("[residents] ", paste(miss_r, collapse = ", ")),
         if (length(miss_i) && length(miss_r)) "; ",
         if (length(miss_i)) paste0("[invaders] ", paste(miss_i, collapse = ", ")))
  }

  # Bandwidth default if needed
  if (is.null(sigma_alpha) || !is.finite(sigma_alpha[1]) || sigma_alpha[1] <= 0) {
    D_rr = D_all[res_ids, res_ids, drop = FALSE]
    med_d = suppressWarnings(stats::median(D_rr[upper.tri(D_rr)], na.rm = TRUE))
    if (!is.finite(med_d) || is.na(med_d)) med_d = 1
    sigma_alpha = max(.Machine$double.eps, med_d)
    if (isTRUE(verbose)) message("Using data-driven sigma_alpha = ", signif(sigma_alpha, 4))
  }
  sig = as.numeric(sigma_alpha[1])

  D_inv_res = D_all[inv_ids, res_ids, drop = FALSE]
  K_inv_res = exp(-(D_inv_res^2) / (2 * sig^2))   # I×J
  C_is_raw  = W_site %*% t(K_inv_res)              # S×I
  dimnames(C_is_raw) = list(sites, inv_ids)

  # Standardize C by site moments
  C_sd_s[!is.finite(C_sd_s) | C_sd_s == 0] = 1
  C_is_z = sweep(sweep(C_is_raw, 1, C_mu_s, `-`), 1, C_sd_s, `/`)

  # Broadcast site-only S
  S_is_z = matrix(S_s_z, nrow = length(sites), ncol = length(inv_ids),
                   dimnames = list(sites, inv_ids))

  ## ---------- Abiotic suitability r_is_z via stats::predict(fit_r) ----------
  # Cross join sites × invaders
  nd = tidyr::expand_grid(site = sites, invader = inv_ids)

  # Left joins (env by site, traits by invader)
  env_df = as.data.frame(env_df_model)
  env_df$site = rownames(env_df_model)
  trt_df = as.data.frame(traits_inv_model)
  trt_df$invader = rownames(traits_inv_model)

  # Prefer many-to-one semantics when available (dplyr >= 1.1)
  if ("relationship" %in% names(formals(dplyr::left_join))) {
    nd = dplyr::left_join(nd, env_df, by = "site",    relationship = "many-to-one")
    nd = dplyr::left_join(nd, trt_df, by = "invader", relationship = "many-to-one")
  } else {
    nd = dplyr::left_join(nd, env_df, by = "site")
    nd = dplyr::left_join(nd, trt_df, by = "invader")
  }

  # Coerce newdata types to match model frame (factor levels, ordered, etc.)
  coerce_like_model = function(newdata, fit) {
    mf = try(stats::model.frame(fit), silent = TRUE)
    if (inherits(mf, "try-error")) return(newdata)
    # Ensure 'site' and 'species' exist with valid levels for random effects
    for (grp in c("site","species")) {
      if (grp %in% names(mf) && !grp %in% names(newdata)) {
        lv = levels(mf[[grp]])
        newdata[[grp]] = factor(rep(lv[1], nrow(newdata)), levels = lv)
      }
    }
    common = intersect(names(mf), names(newdata))
    for (v in common) {
      trn = mf[[v]]
      if (is.factor(trn)) {
        new_lv = levels(trn)
        # redirect unseen to "_other_" if present
        xx = as.character(newdata[[v]])
        if ("_other_" %in% new_lv) xx[!xx %in% new_lv] = "_other_"
        newdata[[v]] = factor(xx, levels = new_lv, ordered = is.ordered(trn))
      } else if (is.numeric(trn)) {
        newdata[[v]] = suppressWarnings(as.numeric(newdata[[v]]))
      } else if (is.integer(trn)) {
        newdata[[v]] = suppressWarnings(as.integer(newdata[[v]]))
      }
    }
    newdata
  }
  nd = coerce_like_model(nd, fit_r)

  # Predict on LINK scale; ignore random effects
  eta = try({
    as.numeric(stats::predict(fit_r, newdata = nd, type = "link", re.form = ~0))
  }, silent = TRUE)

  R_link = matrix(NA_real_, nrow = length(sites), ncol = length(inv_ids),
                   dimnames = list(sites, inv_ids))
  r_is_z = matrix(0, nrow = length(sites), ncol = length(inv_ids),
                   dimnames = list(sites, inv_ids))

  if (!inherits(eta, "try-error") && length(eta) == nrow(nd) && all(is.finite(eta))) {
    R_link[,] = matrix(eta, nrow = length(sites), ncol = length(inv_ids),
                        byrow = FALSE, dimnames = list(sites, inv_ids))
    r_sd_s[!is.finite(r_sd_s) | r_sd_s == 0] = 1
    r_is_z = sweep(sweep(R_link, 1, r_mu_s, `-`), 1, r_sd_s, `/`)
  } else if (isTRUE(verbose)) {
    pred_err = if (inherits(eta, "try-error")) as.character(eta) else "non-finite predictions"
    message("r_is_z could not be built from fit_r; returning zeros.\n",
            "Hints: (a) ensure ENV_PC* and TR_PC* are present when PCA used;\n",
            "       (b) raw trait factor columns must match model levels;\n",
            "       (c) env/trait numeric columns must be numeric.\n",
            "stats::predict error: ", pred_err)
  }

  ## ---------- Tidy long data frame (site × invader) ----------
  pick_mat = function(M) as.numeric(M[cbind(match(nd$site, rownames(M)),
                                             match(nd$invader, colnames(M)))])
  tidy_df = data.frame(
    site    = nd$site,
    invader = nd$invader,
    r_link  = pick_mat(R_link),
    r_z     = pick_mat(r_is_z),
    C_z     = pick_mat(C_is_z),
    S_z     = pick_mat(S_is_z),
    stringsAsFactors = FALSE
  )

  list(
    r_is_z      = r_is_z,
    C_is_z      = C_is_z,
    S_is_z      = S_is_z,
    df          = tidy_df,
    kernels     = list(K_inv_res = K_inv_res, D_inv_res = D_inv_res),
    inputs_used = list(sites = sites, inv_ids = inv_ids, res_ids = res_ids, sigma_alpha = sig)
  )
}
