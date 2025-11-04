#' Build a GLMM-ready model formula from trait and environment terms
#'
#' @title Flexible formula constructor for residents-only trait–environment models
#'
#' @description
#' `build_model_formula()` assembles a right-hand side (RHS) for a GLMM (or LM/GLM)
#' from environment terms, trait terms, and (optionally) all pairwise
#' environment-by-trait interactions. It also appends random-effects structures,
#' such as `(1 | site) + (1 | species)` and optional zero-correlation random
#' slopes like `(0 + r_z || site)`.
#'
#' You can pass the terms directly as character vectors, or let the function
#' derive them from `env_df` and/or `trait_df` column names.
#'
#' @param response Character scalar. Name of the response on the LHS
#'   (default `"abundance"`).
#' @param env_terms Optional character vector of environment term names to include
#'   as fixed effects. If `NULL`, they can be derived from `env_df`.
#' @param trait_terms Optional character vector of trait term names to include
#'   as fixed effects. If `NULL`, they can be derived from `trait_df`.
#' @param env_df Optional data frame containing environment predictors; used only
#'   to infer `env_terms` when `env_terms` is `NULL`.
#' @param trait_df Optional data frame containing trait predictors; used only
#'   to infer `trait_terms` when `trait_terms` is `NULL`.
#' @param include_intercept Logical. Include the fixed-effect intercept?
#'   If `FALSE`, the intercept is removed via `0` (equivalent to `-1`).
#'   Default `TRUE`.
#' @param include_env_main Logical. Include environment main effects?
#'   Default `TRUE`.
#' @param include_trait_main Logical. Include trait main effects?
#'   Default `TRUE`.
#' @param include_env_trait_interactions Logical. Include all pairwise
#'   environment-by-trait interactions? Implemented as
#'   `(E1 + E2 + ...):(T1 + T2 + ...)`. Default `TRUE`.
#' @param extra_fixed Optional character vector of additional fixed-effect terms
#'   to append verbatim (e.g., `"poly(temp,2)"`, `"I(pH^2)"`).
#' @param random_intercepts Character vector of grouping factors for random
#'   intercepts (e.g., `c("site","species")`). Use `NULL` to omit.
#'   Default `c("site","species")`.
#' @param random_slopes Named list of the form
#'   `list(site = c("r_z","C_z"), species = "r_z")` to add zero-correlation
#'   slopes `(0 + term || group)`. Use `NULL` (default) for none.
#' @param backend Character flag used only for messaging; both **lme4**
#'   and **glmmTMB** accept the same syntax here. Default `c("glmmTMB","lme4")`.
#' @param verbose Logical. If `TRUE`, prints the assembled formula string.
#'
#' @return An object of class `formula`, e.g.:
#'   \code{abundance ~ env1 + env2 + tr1 + tr2 + (env1 + env2):(tr1 + tr2) + (1 | site) + (1 | species) + (0 + r_z || site)}
#'
#' @examples
#' \dontrun{
#' # Toy data
#' set.seed(1)
#' env_df_z        = data.frame(env1 = rnorm(10), env2 = rnorm(10))
#' traits_res_glmm = data.frame(tr1  = rnorm(5),  tr2  = rnorm(5))
#'
#' fml = build_model_formula(
#'   response   = "abundance",
#'   env_df     = env_df_z,
#'   trait_df   = traits_res_glmm,
#'   random_intercepts = c("site","species"),
#'   random_slopes     = list(site = c("r_z","C_z"))
#' )
#' fml
#'
#' # Fit a GLMM (example; requires your long residents×sites table `dat_r`)
#' # library(glmmTMB)
#' # fit = glmmTMB::glmmTMB(fml, family = glmmTMB::tweedie(link="log"), data = dat_r)
#' # summary(fit)
#' }
#'
#' @export
build_model_formula = function(
    response = "abundance",
    env_terms = NULL,
    trait_terms = NULL,
    env_df = NULL,
    trait_df = NULL,
    include_intercept = TRUE,
    include_env_main = TRUE,
    include_trait_main = TRUE,
    include_env_trait_interactions = TRUE,
    extra_fixed = NULL,
    random_intercepts = c("site", "species"),
    random_slopes = NULL,
    backend = c("glmmTMB","lme4"),
    verbose = FALSE
){
  backend = match.arg(backend)

  # After `backend = match.arg(backend)`
  if (!is.character(response) || length(response) != 1) {
    stop(
      "First argument 'response' must be a single character string (e.g., 'abundance').\n",
      "Did you mean to call with named arguments like:\n",
      "  build_model_formula(env_df = env_df_z, trait_df = traits_res_glmm)"
    )
  }

  # ---- Helpers ----------------------------------------------------------------
  .as_char = function(x){
    if (is.null(x)) return(character(0))
    if (is.data.frame(x)) return(colnames(x))
    if (is.matrix(x))     return(colnames(x))
    if (is.character(x))  return(unique(stats::na.omit(x)))
    stop("Terms must be character vectors or data frames with named columns.")
  }

  # Derive terms from data frames if not explicitly provided
  if (is.null(env_terms))   env_terms   = .as_char(env_df)
  if (is.null(trait_terms)) trait_terms = .as_char(trait_df)

  env_terms   = .as_char(env_terms)
  trait_terms = .as_char(trait_terms)

  # Build fixed-effect pieces
  fixed_parts = character(0)

  # Intercept control: "0" removes the intercept; otherwise nothing is needed
  if (!isTRUE(include_intercept)) fixed_parts = c(fixed_parts, "0")

  # Main effects
  if (isTRUE(include_env_main)   && length(env_terms))   fixed_parts = c(fixed_parts, env_terms)
  if (isTRUE(include_trait_main) && length(trait_terms)) fixed_parts = c(fixed_parts, trait_terms)

  # All pairwise E:T interactions
  if (isTRUE(include_env_trait_interactions) &&
      length(env_terms) && length(trait_terms)) {
    intx = sprintf("(%s):(%s)",
                    paste(env_terms,   collapse = " + "),
                    paste(trait_terms, collapse = " + "))
    fixed_parts = c(fixed_parts, intx)
  }

  # Any additional user-supplied fixed terms
  if (length(extra_fixed)) fixed_parts = c(fixed_parts, extra_fixed)

  # If nothing made it in, default to "1" (intercept-only) to keep a valid formula
  if (!length(fixed_parts)) fixed_parts = "1"

  # Random intercepts: (1 | group)
  re_parts = character(0)
  if (length(random_intercepts)) {
    stopifnot(is.character(random_intercepts))
    re_parts = c(re_parts, sprintf("(1 | %s)", random_intercepts))
  }

  # Random slopes (zero-correlation): (0 + term || group)
  # Expect a named list: list(site = c("r_z","C_z"), species = "r_z")
  if (!is.null(random_slopes)) {
    stopifnot(is.list(random_slopes))
    for (grp in names(random_slopes)) {
      terms_here = as.character(random_slopes[[grp]])
      terms_here = terms_here[nzchar(terms_here)]
      if (!length(terms_here)) next
      re_parts = c(re_parts, sprintf("(0 + %s || %s)", terms_here, grp))
    }
  }

  # Concatenate all RHS parts
  rhs = paste(c(fixed_parts, re_parts), collapse = " + ")

  # Assemble the full formula
  fml_str = paste0(response, " ~ ", rhs)
  if (isTRUE(verbose)) message("Model formula: ", fml_str)

  stats::as.formula(fml_str)
}
