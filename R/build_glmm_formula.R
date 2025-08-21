#' Build a trait-environment GLMM formula safely and flexibly
#'
#' @description
#' Constructs a model formula for trait-environment analyses in a single step.
#' The function (i) auto-detects trait and environment columns from a long-format
#' table, (ii) assembles fixed effects for all traits and all environment variables,
#' (iii) optionally includes all pairwise \eqn{trait \times environment} interactions,
#' and (iv) appends user-specified random-effects terms. The returned object is a
#' standard \code{formula} suitable for \pkg{glmmTMB}, \pkg{lme4}, etc.
#'
#' @param data \code{data.frame}. Long-format observations (e.g., species-by-site),
#'   including the response, species ID, site ID, trait columns, and environment columns.
#' @param response \code{character} (default \code{"count"}).
#'   Name of the response variable (e.g., count/abundance).
#' @param species_col \code{character} (default \code{"species"}).
#'   Column name identifying species.
#' @param site_col \code{character} (default \code{"site_id"}).
#'   Column name identifying sites.
#' @param trait_cols \code{NULL} or \code{character} vector.
#'   If \code{NULL} (default), traits are auto-detected using name prefixes
#'   \verb{^trait_}, \verb{^t_}, or \verb{^trt_}. If not found, falls back to
#'   "everything not excluded" (see \code{env_exclude}). Pass explicit names
#'   for full control.
#' @param env_cols \code{NULL} or \code{character} vector.
#'   If \code{NULL} (default), environment variables are auto-detected using name
#'   prefixes \verb{^env_}, \verb{^e_}, \verb{^clim_}, \verb{^soil_}. If not found,
#'   falls back to "everything not in traits and not excluded".
#' @param env_exclude \code{character} vector.
#'   Columns to exclude from environment auto-detection. Defaults to
#'   \code{c("site_id","x","y","count","species")}. Adjust to your schema.
#' @param include_interactions \code{logical} (default \code{TRUE}).
#'   If \code{TRUE}, adds a single block term \code{(traits):(envs)} which expands
#'   to all pairwise \eqn{trait \times environment} interactions.
#' @param random_effects \code{character} vector.
#'   Random-effect terms to append to the RHS (e.g., \code{"(1 | species)"}).
#'   Use \code{character(0)} to omit random effects. Default adds random intercepts
#'   for species and site: \code{c("(1 | species)", "(1 | site_id)")}.
#'
#' @return A \code{formula} with fixed effects (traits + envs [+ interactions])
#'   and any requested random effects, e.g.:
#'   \preformatted{
#'   count ~ trait_cont1 + ... + trait_cat + env1 + ... + envK +
#'           (trait_cont1 + ... + trait_cat):(env1 + ... + envK) +
#'           (1 | species) + (1 | site_id)
#'   }
#'
#' @details
#' \strong{Auto-detection:}
#' - Traits: first tries prefixes \verb{^trait_}, \verb{^t_}, \verb{^trt_}. If none match,
#'   uses all columns not in \code{env_exclude}, not \code{response}, not \code{species_col},
#'   and not \code{site_col}.
#' - Environment: first tries prefixes \verb{^env_}, \verb{^e_}, \verb{^clim_}, \verb{^soil_}.
#'   If none match, uses remaining non-excluded columns not already assigned as traits.
#'
#' \strong{Interactions:}
#' When \code{include_interactions = TRUE}, a single block term
#' \code{(t1 + t2 + ...):(e1 + e2 + ...)} is inserted; model-fitting packages will
#' expand it to all pairwise interactions. Disable with \code{FALSE} if the design
#' is too large or you prefer targeted interactions.
#'
#' \strong{Random effects:}
#' Supplied verbatim (e.g., random intercepts/slopes). For example,
#' \code{c("(1 | species)", "(1 | site_id)")} or \code{c("(1 + key_trait | species)")}.
#'
#' @examples
#' # Minimal reproducible toy example -----------------------------------------
#' set.seed(1)
#' n <- 100
#' longDF <- data.frame(
#'   site_id = factor(sample(paste0("s", 1:10), n, TRUE)),
#'   species = factor(sample(paste0("sp", 1:15), n, TRUE)),
#'   x = runif(n), y = runif(n),
#'   count = rpois(n, lambda = 3),
#'   # traits
#'   trait_cont1 = rnorm(n),
#'   trait_cont2 = rnorm(n),
#'   trait_cat = factor(sample(letters[1:3], n, TRUE)),
#'   # environments
#'   env1 = scale(rnorm(n))[, 1],
#'   env2 = scale(runif(n))[, 1]
#' )
#'
#' # Build a full formula with all trait x environment interactions and default REs
#' fml <- build_glmm_formula(longDF)
#' fml
#'
#' # Example fit (uncomment if glmmTMB is available)
#' # mod = glmmTMB::glmmTMB(fml, data = longDF, family = glmmTMB::tweedie(link = "log"))
#' # summary(mod)
#'
#' # Targeted columns & no interactions
#' fml2 <- build_glmm_formula(
#'   data = longDF,
#'   trait_cols = c("trait_cont1", "trait_cont2", "trait_cat"),
#'   env_cols = c("env1", "env2"),
#'   include_interactions = FALSE,
#'   random_effects = character(0)
#' )
#' fml2
#'
#' @seealso \code{\link[glmmTMB]{glmmTMB}}, \code{\link[lme4]{lmer}}, \code{\link[stats]{as.formula}}
#' @importFrom stats as.formula
#' @export
build_glmm_formula <- function(
    data,
    response = "count",
    species_col = "species",
    site_col = "site_id",
    # If NULL, auto-detect below. Otherwise pass character vectors of names.
    trait_cols = NULL,
    env_cols = NULL,
    # Columns to exclude from env auto-detection (IDs, coords, response).
    env_exclude = c("site_id", "x", "y", "count", "species"),
    include_interactions = TRUE,
    # Random effects as character terms; set to character(0) to omit.
    random_effects = c("(1 | species)", "(1 | site_id)")) {
  # --- Basic checks ----------------------------------------------------------
  stopifnot(is.data.frame(data))
  req <- c(response, species_col, site_col)
  missing_req <- setdiff(req, names(data))
  if (length(missing_req)) {
    stop(
      "Missing required columns in 'data': ",
      paste(missing_req, collapse = ", ")
    )
  }

  # --- Auto-detect TRAIT columns --------------------------------------------
  # Heuristic 1: prefer common trait prefixes; adjust to your naming scheme.
  if (is.null(trait_cols)) {
    trait_guess <- grep("^(trait_|t_|trt_)", names(data), value = TRUE)
    # Heuristic 2 (fallback): everything not excluded and not IDs/response.
    if (!length(trait_guess)) {
      trait_guess <- setdiff(
        names(data),
        unique(c(env_exclude, response, species_col, site_col))
      )
    }
    trait_cols <- trait_guess
  }

  # --- Auto-detect ENV columns ----------------------------------------------
  # Heuristic 1: prefer common env prefixes; adjust as needed.
  if (is.null(env_cols)) {
    env_guess <- grep("^(env_|e_|clim_|soil_|bio_)", names(data), value = TRUE)
    # Heuristic 2 (fallback): remaining non-excluded, non-trait columns.
    if (!length(env_guess)) {
      env_guess <- setdiff(
        names(data),
        unique(c(trait_cols, env_exclude, response, species_col, site_col))
      )
    }
    env_cols <- env_guess
  }

  # --- Final sanity checks on detected sets ---------------------------------
  trait_cols <- intersect(trait_cols, names(data))
  env_cols <- intersect(env_cols, names(data))
  if (!length(trait_cols)) stop("No trait columns detected; supply 'trait_cols'.")
  if (!length(env_cols)) stop("No environment columns detected; supply 'env_cols'.")

  # --- Build FIXED-EFFECT terms ---------------------------------------------
  # Concatenate trait and env main-effect terms.
  trait_terms <- paste(trait_cols, collapse = " + ")
  env_terms <- paste(env_cols, collapse = " + ")

  # Interactions as a single block (t1 + t2 + ...):(e1 + e2 + ...)
  inter_terms <- if (isTRUE(include_interactions)) {
    paste0("(", trait_terms, "):(", env_terms, ")")
  } else {
    NULL
  }

  # --- RANDOM-EFFECT terms ---------------------------------------------------
  # Supplied verbatim; allows intercepts/slopes. Omit with character(0).
  re_terms <- character(0)
  if (length(random_effects)) re_terms <- random_effects

  # --- Assemble RHS and full formula ----------------------------------------
  rhs_parts <- c(trait_terms, env_terms, inter_terms, re_terms)
  rhs_parts <- rhs_parts[!is.na(rhs_parts) & nzchar(rhs_parts)]
  rhs <- paste(rhs_parts, collapse = " + ")

  fml_chr <- paste(response, "~", rhs)

  # --- Return standard formula ----------------------------------------------
  stats::as.formula(fml_chr)
}
