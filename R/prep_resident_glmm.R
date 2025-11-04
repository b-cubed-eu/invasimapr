# ======================================================================
# Fast residents long table + GLMM fit (no-copy joins, optional data.table)
# ======================================================================

#' Build residents-by-sites long table and fit the residents-only GLMM
#'
#' @description
#' Creates a long data set with one row per (site, resident) pair, attaches
#' standardized environment and resident traits, then fits a GLMM for residents.
#' Also returns fixed-effect predictions on the link scale as a site×resident
#' matrix (r_js) and the response-scale mean (mu_js).
#'
#' @param comm_res Numeric matrix (sites × residents) of observed counts or abundance.
#'   Row names are site IDs; column names are resident species IDs.
#' @param env_df_z Data frame of standardized site-level environment.
#'   Row names are site IDs and must match `rownames(comm_res)`.
#' @param traits_res_glmm Data frame of standardized resident traits used in the GLMM.
#'   Row names are resident IDs and must match `colnames(comm_res)`.
#' @param fml A model formula for `glmmTMB::glmmTMB` (e.g., from `build_model_formula()`).
#' @param family A `glmmTMB` family object. Default is `glmmTMB::tweedie(link="log")`.
#' @param seed Integer. Random seed.
#' @param fit_model logical; if FALSE, return dat_r/fml without fitting.
#'
#' @return list(dat_r, fit_r, grid_res, r_js, mu_js, sites, res_ids)
#' @import data.table
#' @importFrom glmmTMB glmmTMB tweedie
#' @export
prep_resident_glmm = function(comm_res,
                               env_df_z,
                               traits_res_glmm,
                               fml,
                               family = glmmTMB::tweedie(link = "log"),
                               seed = 123,
                               fit_model = TRUE) {

  # --- validate & align --------------------------------------------------------
  if (inherits(comm_res, "Matrix")) comm_res = as.matrix(comm_res)
  if (!is.matrix(comm_res))        comm_res = as.matrix(comm_res)
  storage.mode(comm_res) = "double"

  if (is.null(rownames(comm_res)) || is.null(colnames(comm_res))) {
    stop("comm_res must have rownames (sites) and colnames (resident IDs).")
  }
  sites   = rownames(comm_res)
  res_ids = colnames(comm_res)

  if (is.null(rownames(env_df_z))) stop("env_df_z must have rownames = site IDs.")
  if (!identical(sites, rownames(env_df_z))) {
    if (setequal(sites, rownames(env_df_z))) {
      env_df_z = env_df_z[sites, , drop = FALSE]
    } else {
      stop("Row names of env_df_z do not match sites in comm_res.")
    }
  }

  if (is.null(rownames(traits_res_glmm))) {
    stop("traits_res_glmm must have rownames = resident IDs.")
  }
  if (!all(res_ids %in% rownames(traits_res_glmm))) {
    miss = setdiff(res_ids, rownames(traits_res_glmm))
    stop("traits_res_glmm is missing rows for: ", paste(miss, collapse = ", "))
  }
  traits_res_glmm = traits_res_glmm[res_ids, , drop = FALSE]

  # --- ultra-fast long table construction -------------------------------------
  .use_dt = requireNamespace("data.table", quietly = TRUE)

  if (.use_dt) {
    # data.table path (very fast, very memory efficient)
    dt = data.table::CJ(site = sites, species = res_ids, sorted = FALSE)
    dt[, abundance := as.vector(comm_res[cbind(match(site, sites), match(species, res_ids))])]

    env_tbl_dt    = data.table::data.table(site = rownames(env_df_z),    env_df_z,    check.names = FALSE)
    traits_tbl_dt = data.table::data.table(species = rownames(traits_res_glmm), traits_res_glmm, check.names = FALSE)

    dt = env_tbl_dt[dt, on = "site"]
    dt = traits_tbl_dt[dt, on = "species"]

    dat_r = as.data.frame(dt, check.names = FALSE)
    dat_r$site    = factor(dat_r$site,    levels = sites)
    dat_r$species = factor(dat_r$species, levels = res_ids)

  } else {
    # Base-R "indexed fill" path (no large copies, still fast)
    dat_r = tidyr::expand_grid(site = sites, species = res_ids)
    ri = match(dat_r$site,    rownames(comm_res))
    cj = match(dat_r$species, colnames(comm_res))
    dat_r$abundance = comm_res[cbind(ri, cj)]

    env_tbl    = data.frame(site    = rownames(env_df_z),    env_df_z,    row.names = NULL, check.names = FALSE)
    traits_tbl = data.frame(species = rownames(traits_res_glmm), traits_res_glmm, row.names = NULL, check.names = FALSE)

    env_cols   = setdiff(names(env_tbl), "site")
    trait_cols = setdiff(names(traits_tbl), "species")
    # preallocate
    for (nm in env_cols)   dat_r[[nm]] = NA_real_
    for (nm in trait_cols) dat_r[[nm]] = NA_real_

    idx_env = match(dat_r$site,    env_tbl$site)
    idx_tr  = match(dat_r$species, traits_tbl$species)
    for (nm in env_cols)   dat_r[[nm]] = env_tbl[[nm]][idx_env]
    for (nm in trait_cols) dat_r[[nm]] = traits_tbl[[nm]][idx_tr]

    dat_r$site    = factor(dat_r$site,    levels = sites)
    dat_r$species = factor(dat_r$species, levels = res_ids)
  }

  # --- formula sanity ----------------------------------------------------------
  need = setdiff(all.vars(fml), names(dat_r))
  if (length(need)) stop("Your formula references variables not found in the residents dataset: ",
                         paste(need, collapse = ", "))

  if (!isTRUE(fit_model)) {
    return(list(dat_r = dat_r, fml = fml, sites = sites, res_ids = res_ids))
  }

  # --- fit GLMM ---------------------------------------------------------------
  set.seed(seed)
  fit_r = glmmTMB::glmmTMB(formula = fml, data = dat_r, family = family)

  # --- predictions (fixed effects only) ---------------------------------------
  grid_res = unique(dat_r[, c("site","species")])
  grid_res = dplyr::left_join(grid_res, tibble::rownames_to_column(env_df_z, "site"), by = "site")
  grid_res = dplyr::left_join(grid_res, tibble::rownames_to_column(traits_res_glmm, "species"), by = "species")
  grid_res = grid_res[, c("site","species", setdiff(names(grid_res), c("site","species"))), drop = FALSE]

  grid_res$eta_r = stats::predict(fit_r, newdata = grid_res, type = "link",
                                   re.form = NA, allow.new.levels = TRUE)

  r_js = matrix(NA_real_, nrow = length(sites), ncol = length(res_ids),
                 dimnames = list(sites, res_ids))
  ii = match(grid_res$site, sites); jj = match(grid_res$species, res_ids)
  r_js[cbind(ii, jj)] = grid_res$eta_r
  mu_js = exp(r_js)

  list(dat_r = dat_r, fit_r = fit_r, grid_res = grid_res, r_js = r_js, mu_js = mu_js,
       sites = sites, res_ids = res_ids)
}
