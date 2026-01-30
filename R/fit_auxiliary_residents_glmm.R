#' Fit an auxiliary residents-only GLMM on standardized predictors
#'
#' @title Auxiliary GLMM for trait-varying and site-varying slopes
#'
#' @description
#' \code{fit_auxiliary_residents_glmm()} assembles a long table of resident
#' site by species cells with standardized predictors
#' \eqn{r^{(z)}_{js}}, \eqn{C^{(z)}_{js}}, and \eqn{S^{(z)}_{js}}, together with
#' two-dimensional trait scores (\code{tr1}, \code{tr2}). It then fits a Gaussian
#' generalized linear mixed model to estimate how slopes on abiotic suitability,
#' niche crowding, and site saturation vary across the trait plane.
#'
#' Optional site-level random slopes can be included for \code{r_z} and
#' \code{C_z}. Random slopes for \code{S_z} are intentionally excluded because
#' \code{S_z} is site-only and has no within-site variation across residents.
#'
#' This function is intentionally forgiving: non-matrix inputs are coerced using
#' \code{as.matrix()}, missing dimension names are repaired when possible, and
#' inputs are aligned to \code{comm_res}.
#'
#' @param comm_res Numeric matrix of site by resident abundances. Row names are
#'   site identifiers and column names are resident identifiers.
#' @param r_js_z Numeric matrix of standardized abiotic suitability with the same
#'   dimensions and names as \code{comm_res}.
#' @param C_js_z Numeric matrix of standardized niche crowding with the same
#'   dimensions and names as \code{comm_res}.
#' @param S_js_z Numeric matrix of standardized site saturation broadcast across
#'   residents, with the same dimensions and names as \code{comm_res}.
#' @param Q_res Data frame of resident trait scores. Must contain numeric columns
#'   \code{tr1} and \code{tr2}, with row names matching
#'   \code{colnames(comm_res)}.
#' @param use_site_random_slopes Logical; if \code{TRUE}, include random slopes
#'   for \code{r_z} and \code{C_z} at the site level.
#' @param family GLM family for the conditional model. Default is
#'   \code{gaussian()} applied to \code{log1p(abundance)}.
#' @param control Control object passed to \code{glmmTMB::glmmTMB()}.
#' @param na_action How to handle rows with missing values. Either
#'   \code{"drop"} (default) or \code{"error"}.
#' @param verbose Logical; if \code{TRUE}, print a short model summary line.
#'
#' @details
#' The fixed-effects component of the model is
#' \deqn{\log(1 + \mathrm{abundance}) \sim (r_z + C_z + S_z) \times (tr1 + tr2),}
#' with random intercepts for species and site always included.
#'
#' When \code{use_site_random_slopes = TRUE}, random slopes
#' \eqn{(0 + r_z | site)} and \eqn{(0 + C_z | site)} are added. Random slopes
#' for \code{S_z} are omitted because \code{S_z} is constant within each site
#' and therefore not identifiable.
#'
#' @section Role in invasion fitness:
#' This auxiliary model estimates trait-dependent slope systems that can be
#' mapped to invader-level parameters such as
#' \eqn{\theta_i}, \eqn{\alpha_i}, and \eqn{\beta_i}, and optionally to
#' site-varying counterparts via random slopes. These parameters enter directly
#' into the linear invasion-fitness decomposition
#' \deqn{\lambda_{is} = \Gamma_{is} r^{(z)}_{is} -
#'       \alpha_{is} C^{(z)}_{is} -
#'       \beta_i S^{(z)}_{is},}
#' linking trait position to abiotic suitability, niche crowding, and resident
#' competition.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{fit}: the fitted \code{glmmTMB} model
#'   \item \code{data}: the long-format data frame used for fitting
#'   \item \code{formula}: the model formula
#'   \item \code{args}: resolved arguments used in the fit
#' }
#'
#' @importFrom stats complete.cases
#' @importFrom utils head
#' @importFrom stats gaussian
#'
#' @examples
#' \dontrun{
#' # Toy dimensions
#' set.seed(1)
#' S = 8  # sites
#' J = 6  # residents
#' sites   = paste0("s", 1:S)
#' res_ids = paste0("sp", 1:J)
#'
#' # Site x resident abundance
#' comm_res = matrix(rpois(S*J, lambda = 2), S, J,
#'                    dimnames = list(sites, res_ids))
#'
#' # Standardized predictors (row-wise z by site for r and C; site-only for S)
#' r_raw = matrix(rnorm(S*J), S, J, dimnames = dimnames(comm_res))
#' C_raw = matrix(rnorm(S*J), S, J, dimnames = dimnames(comm_res))
#' S_raw = matrix(rep(scale(rnorm(S))[,1], each = J), S, J,
#'                 dimnames = dimnames(comm_res))
#' # Simple per-site z for demo
#' r_js_z = t(scale(t(r_raw))); r_js_z[!is.finite(r_js_z)] = 0
#' C_js_z = t(scale(t(C_raw))); C_js_z[!is.finite(C_js_z)] = 0
#' S_js_z = S_raw
#'
#' # Resident trait-plane scores (PCoA on Gower in real workflow)
#' Q_res = data.frame(tr1 = rnorm(J), tr2 = rnorm(J))
#' rownames(Q_res) = res_ids
#'
#' aux = fit_auxiliary_residents_glmm(
#'   comm_res   = comm_res,
#'   r_js_z     = r_js_z,
#'   C_js_z     = C_js_z,
#'   S_js_z     = S_js_z,
#'   Q_res      = Q_res,
#'   use_site_random_slopes = TRUE
#' )
#'
#' summary(aux$fit)
#' head(aux$data)
#' aux$formula
#' }
#'
#' @importFrom stats complete.cases
#' @importFrom utils head
#' @export
fit_auxiliary_residents_glmm = function(
    comm_res,
    r_js_z,
    C_js_z,
    S_js_z,
    Q_res,
    use_site_random_slopes = TRUE,
    family  = gaussian(),
    control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1e5, eval.max = 1e5)),
    na_action = c("drop", "error"),
    verbose = TRUE
){
  na_action = match.arg(na_action)

  ## ---- Robust coercion to base matrices --------------------------------------
  to_base_matrix = function(x, name) {
    # Accept Matrix, data.frame, tibble, data.table, etc.
    if (inherits(x, "Matrix")) x = as.matrix(x)
    else if (!is.matrix(x))    x = as.matrix(x)

    if (!is.matrix(x)) stop(sprintf("%s could not be coerced to a matrix.", name))

    # Ensure 2D numeric with finite dims
    if (length(dim(x)) != 2L || any(dim(x) == 0L))
      stop(sprintf("%s must be a non-empty 2D matrix; got dim = (%s).",
                   name, paste(dim(x), collapse = "x")))
    suppressWarnings(storage.mode(x) == "double")
    x
  }

  comm_res = to_base_matrix(comm_res, "comm_res")
  r_js_z   = to_base_matrix(r_js_z,   "r_js_z")
  C_js_z   = to_base_matrix(C_js_z,   "C_js_z")
  S_js_z   = to_base_matrix(S_js_z,   "S_js_z")

  ## ---- Repair/propagate dimnames --------------------------------------------
  if (is.null(rownames(comm_res))) rownames(comm_res) = paste0("s", seq_len(nrow(comm_res)))
  if (is.null(colnames(comm_res))) colnames(comm_res) = paste0("sp", seq_len(ncol(comm_res)))

  propagate_dimnames = function(A, ref) {
    if (is.null(rownames(A)) && nrow(A) == nrow(ref)) rownames(A) = rownames(ref)
    if (is.null(colnames(A)) && ncol(A) == ncol(ref)) colnames(A) = colnames(ref)
    A
  }
  r_js_z = propagate_dimnames(r_js_z, comm_res)
  C_js_z = propagate_dimnames(C_js_z, comm_res)
  S_js_z = propagate_dimnames(S_js_z, comm_res)

  same_dimnames = function(A, B) {
    identical(dim(A), dim(B)) &&
      identical(rownames(A), rownames(B)) &&
      identical(colnames(A), colnames(B))
  }
  if (!same_dimnames(r_js_z, comm_res)) stop("r_js_z must align with comm_res (rows/cols and names).")
  if (!same_dimnames(C_js_z, comm_res)) stop("C_js_z must align with comm_res (rows/cols and names).")
  if (!same_dimnames(S_js_z, comm_res)) stop("S_js_z must align with comm_res (rows/cols and names).")

  ## ---- Q_res checks -----------------------------------------------------------
  stopifnot(is.data.frame(Q_res))
  if (!all(c("tr1","tr2") %in% colnames(Q_res)))
    stop("Q_res must contain columns 'tr1' and 'tr2'.")
  if (is.null(rownames(Q_res)) && "species" %in% names(Q_res)) {
    rn = as.character(Q_res$species)
    if (anyDuplicated(rn)) stop("Q_res$species contains duplicate IDs; cannot set as rownames.")
    rownames(Q_res) = rn
  }

  sites   = rownames(comm_res)
  res_ids = colnames(comm_res)
  if (!all(res_ids %in% rownames(Q_res))) {
    missing = setdiff(res_ids, rownames(Q_res))
    stop("Q_res is missing trait rows for residents: ",
         paste(utils::head(missing, 10), collapse = ", "),
         if (length(missing) > 10) " ...")
  }

  ## ---- Build long table -------------------------------------------------------
  dat2 = tidyr::expand_grid(site = sites, species = res_ids)

  idx = cbind(match(dat2$site, sites), match(dat2$species, res_ids))
  dat2$abundance = comm_res[idx]
  dat2$r_z       = r_js_z[idx]
  dat2$C_z       = C_js_z[idx]
  dat2$S_z       = S_js_z[idx]

  Q_res2 = Q_res[res_ids, , drop = FALSE]
  Q_res2 = tibble::rownames_to_column(Q_res2, "species")

  dat2 = dat2 |>
    dplyr::left_join(Q_res2, by = "species") |>
    dplyr::mutate(
      site    = factor(site,    levels = sites),
      species = factor(species, levels = res_ids)
    )

  needed = c("abundance","r_z","C_z","S_z","tr1","tr2")
  bad = !stats::complete.cases(dat2[, needed])
  if (any(bad)) {
    msg = sprintf("Found %d rows with missing values in required columns.", sum(bad))
    if (na_action == "error") stop(msg)
    if (isTRUE(verbose)) message(msg, " Dropping those rows for the fit.")
    dat2 = dat2[!bad, , drop = FALSE]
  }

  ## ---- Formula & fit ----------------------------------------------------------
  fixed = "(r_z + C_z + S_z) * (tr1 + tr2)"
  rand  = c("(1 | species)", "(1 | site)")
  if (isTRUE(use_site_random_slopes)) {
    rand = c(rand, "(0 + r_z || site)", "(0 + C_z || site)")
  }
  rhs = paste(c(fixed, rand), collapse = " + ")
  fml = stats::as.formula(paste0("log1p(abundance) ~ ", rhs))

  fit = glmmTMB::glmmTMB(
    formula = fml,
    data    = dat2,
    family  = family,
    control = glmmTMB::glmmTMBControl(
      optCtrl = list(iter.max = 1e5, eval.max = 1e5),
      profile = TRUE) # sometimes stabilizes variance comps
  )

  if (isTRUE(verbose)) {
    fe = try(names(glmmTMB::fixef(fit)$cond), silent = TRUE)
    message("Auxiliary GLMM fitted: ",
            if (!inherits(fe, "try-error")) paste0(length(fe), " fixed terms; ") else "",
            "random: ", paste(rand, collapse = " + "))
  }

  list(
    fit     = fit,
    data    = dat2,
    formula = fml,
    args    = list(
      use_site_random_slopes = use_site_random_slopes,
      family  = family,
      control = control,
      na_action = na_action
    )
  )
}

