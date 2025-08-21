#' Predict site-level responses for residents and simulated invaders
#'
#' @description
#' Builds a complete **site x species** prediction grid by crossing a species-trait table
#' (residents + simulated invaders) with a site-environment table, then calls the fitted
#' model's `predict()` to obtain expected responses (e.g., abundance). The function is
#' robust to factor level issues and mirrors a typical `expand_grid()` + left-join workflow.
#'
#' @details
#' **Model classes supported**
#' \itemize{
#'   \item \strong{glmmTMB}: uses `type` and `re.form`. Population-level predictions set
#'         `re.form = ~ 0` (default, `include_random = FALSE`); conditional predictions
#'         use `re.form = NULL` with `allow.new.levels = TRUE`.
#'   \item \strong{lme4::merMod} (`lmer`/`glmer`): population-level uses `re.form = NA`;
#'         conditional uses `re.form = NULL` with `allow.new.levels = TRUE`.
#'   \item \strong{mgcv::gam}: passes `type` through (e.g., `"response"`).
#'   \item \strong{stats::glm} / \strong{stats::lm}: standard `predict()`; `type` for GLM, none for LM.
#' }
#'
#' **Why population-level predictions by default?** For novel invaders, random-effect
#' levels (e.g., species/site intercepts) are unknown. Excluding random effects yields
#' fixed-effects expectations driven by traits, environments, and their interactions -
#' appropriate for invasion screening and ranking.
#'
#' **Augmentation inputs.** `site_aug` and/or `species_aug` let you add extra predictors
#' (e.g., `obs_sum`, `spp_rich`) prior to prediction, replicating the common
#' `left_join()` pattern used in analysis notebooks.
#'
#' @param model Fitted model object (\pkg{glmmTMB}, \pkg{lme4}, \pkg{mgcv}, or base \pkg{stats}).
#' @param species_traits \code{data.frame}. Species traits (one row per species) containing
#'   \code{species_col} and all trait variables referenced in the model's fixed effects.
#' @param site_env \code{data.frame}. Site predictors (one row per site) containing
#'   \code{site_col} and all environmental variables referenced in the model's fixed effects.
#' @param species_col \code{character}. Species ID column in \code{species_traits}. Default \code{"species"}.
#' @param site_col \code{character}. Site ID column in \code{site_env}. Default \code{"site_id"}.
#' @param response_type \code{character}. Prediction scale for \code{predict()} where applicable
#'   (e.g., \code{"response"} or \code{"link"}). Default \code{"response"}.
#' @param include_random \code{logical}. If \code{FALSE} (default), compute population-level
#'   predictions (exclude random effects). If \code{TRUE}, compute conditional predictions
#'   (include random effects; allow new levels where supported).
#' @param site_aug \code{NULL} or \code{data.frame}. Optional site-level augmentation table
#'   (must contain \code{site_col}) to be left-joined into \code{site_env} before prediction.
#' @param species_aug \code{NULL} or \code{data.frame}. Optional species-level augmentation table
#'   (must contain \code{species_col}) to be left-joined into \code{species_traits}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{newdata}: the site x species table used in \code{predict()}.
#'   \item \code{predictions}: long table with columns \code{site_col}, \code{species_col}, and \code{pred}.
#'   \item \code{prediction_matrix}: wide matrix (sites x species) of predicted values
#'         (rows ordered by site ID, columns by species ID).
#' }
#'
#' @importFrom stats formula predict
#'
#' @examples
#' \dontrun{
#' # Residents + invaders (traits), sites (env), fitted Tweedie glmmTMB model:
#' all_traits <- rbind(spp_trait, inv_traits)
#' out <- predict_invader_response(
#'   model          = mod,
#'   species_traits = all_traits,
#'   site_env       = site_env,
#'   species_col    = "species",
#'   site_col       = "site_id",
#'   response_type  = "response",
#'   include_random = FALSE, # population-level
#'   site_aug       = dplyr::select(spp_rich_obs, site_id, obs_sum, spp_rich)
#' )
#' head(out$predictions)
#' dim(out$prediction_matrix)
#' }
#'
#' @seealso \code{\link[glmmTMB]{predict.glmmTMB}}, \code{\link[lme4]{predict.merMod}},
#'   \code{\link[mgcv]{predict.gam}}, \code{\link[stats]{predict}}
#' @export
predict_invader_response <- function(
    model,
    species_traits, # residents + invaders (must include species_col)
    site_env, # site-level predictors (must include site_col)
    species_col = "species",
    site_col = "site_id",
    response_type = "response",
    include_random = FALSE, # population-level by default
    site_aug = NULL, # optional extra site-level columns (e.g., obs_sum, spp_rich)
    species_aug = NULL # optional extra species-level columns
    ) {
  # ---- Basic checks ----------------------------------------------------------
  stopifnot(is.data.frame(species_traits), is.data.frame(site_env))
  if (!species_col %in% names(species_traits)) stop("species_col not found in species_traits.")
  if (!site_col %in% names(site_env)) stop("site_col not found in site_env.")

  # ---- Optional augments (mirror typical dplyr joins) ------------------------
  if (!is.null(site_aug)) {
    if (!is.data.frame(site_aug) || !site_col %in% names(site_aug)) {
      stop("site_aug must be a data.frame containing column: ", site_col)
    }
    site_env <- dplyr::left_join(site_env, site_aug, by = site_col)
  }
  if (!is.null(species_aug)) {
    if (!is.data.frame(species_aug) || !species_col %in% names(species_aug)) {
      stop("species_aug must be a data.frame containing column: ", species_col)
    }
    species_traits <- dplyr::left_join(species_traits, species_aug, by = species_col)
  }

  # ---- Construct site x species grid robustly --------------------------------
  sp_vals <- unique(species_traits[[species_col]])
  site_vals <- unique(site_env[[site_col]])

  # Use tidyr::expand_grid if available; else replicate manually
  if (requireNamespace("tidyr", quietly = TRUE)) {
    grid <- tidyr::expand_grid(
      !!species_col := sp_vals,
      !!site_col := site_vals
    )
  } else {
    grid <- data.frame(
      species_tmp = rep(sp_vals, each = length(site_vals)),
      site_tmp = rep(site_vals, times = length(sp_vals)),
      stringsAsFactors = FALSE
    )
    names(grid) <- c(species_col, site_col)
  }

  # Coerce join key types to match sources (avoids failed joins)
  if (inherits(species_traits[[species_col]], "factor")) {
    grid[[species_col]] <- factor(grid[[species_col]], levels = levels(species_traits[[species_col]]))
  } else {
    storage.mode(grid[[species_col]]) <- storage.mode(species_traits[[species_col]])
  }
  if (inherits(site_env[[site_col]], "factor")) {
    grid[[site_col]] <- factor(grid[[site_col]], levels = levels(site_env[[site_col]]))
  } else {
    storage.mode(grid[[site_col]]) <- storage.mode(site_env[[site_col]])
  }

  # ---- Join traits and site covariates (same logic as your working code) -----
  newdata <- dplyr::left_join(grid, species_traits, by = species_col) %>%
    dplyr::left_join(site_env, by = site_col)

  # ---- Fixed-effects terms & factor-level harmonisation ----------------------
  # Try to recover the training frame (handles glmmTMB/lme4/glm/lm/gam)
  train_df <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (is.null(train_df) && inherits(model, "glmmTMB")) {
    train_df <- tryCatch(model$frame, error = function(e) NULL)
  }

  # Extract fixed-effects formula per model class
  fe_formula <- tryCatch(
    {
      if (inherits(model, "glmmTMB")) {
        stats::formula(model)$cond
      } else if (inherits(model, "merMod")) {
        lme4::formula(model, fixed.only = TRUE)
      } else {
        stats::formula(model)
      }
    },
    error = function(e) stats::formula(model)
  )

  # Expand interaction labels to underlying variables
  rhs_labels <- tryCatch(attr(stats::terms(fe_formula), "term.labels"), error = function(e) NULL)
  rhs_all <- if (length(rhs_labels)) {
    unique(all.vars(stats::terms(stats::as.formula(paste("~", paste(rhs_labels, collapse = "+"))))))
  } else {
    unique(all.vars(fe_formula))
  }

  # Coerce non-ID fixed-effect factors in newdata to training levels (prevents MM mismatch)
  if (!is.null(train_df)) {
    for (v in intersect(rhs_all, intersect(names(train_df), names(newdata)))) {
      if (v %in% c(site_col, species_col)) next # never touch ID columns
      trv <- train_df[[v]]
      if (is.factor(trv)) {
        newdata[[v]] <- if (is.ordered(trv)) {
          factor(newdata[[v]], levels = levels(trv), ordered = TRUE)
        } else {
          factor(newdata[[v]], levels = levels(trv))
        }
      }
    }
  }

  # ---- Pre-check: ensure all fixed-effect variables are present --------------
  needed <- setdiff(rhs_all, c("(Intercept)"))
  missing_vars <- setdiff(needed, names(newdata))
  if (length(missing_vars)) {
    stop(
      "Missing fixed-effect variable(s) in newdata: ",
      paste(missing_vars, collapse = ", "),
      ". Add them via site_env/species_traits or site_aug/species_aug."
    )
  }

  # ---- Predict (with glmmTMB-safe fallback) ----------------------------------
  if (inherits(model, "glmmTMB")) {
    # Try population-level (or conditional) first; if it errors, fall back to allow.new.levels only
    pred <- try(
      predict(
        model,
        newdata = newdata, type = response_type,
        re.form = if (isTRUE(include_random)) NULL else ~0,
        allow.new.levels = TRUE
      ),
      silent = TRUE
    )
    if (inherits(pred, "try-error")) {
      pred <- predict(model, newdata = newdata, type = response_type, allow.new.levels = TRUE)
    }
  } else if (inherits(model, "merMod")) {
    pred <- if (isTRUE(include_random)) {
      stats::predictMerMod(model,
        newdata = newdata, type = response_type,
        re.form = NULL, allow.new.levels = TRUE
      )
    } else {
      stats::predictMerMod(model, newdata = newdata, type = response_type, re.form = NA)
    }
  } else if (inherits(model, "gam")) {
    pred <- mgcv::predict.gam(model, newdata = newdata, type = response_type)
  } else if (inherits(model, "glm")) {
    pred <- stats::predict(model, newdata = newdata, type = response_type)
  } else if (inherits(model, "lm")) {
    pred <- stats::predict(model, newdata = newdata) # no 'type'
  } else {
    stop(
      "Unsupported model class: ", paste(class(model), collapse = ", "),
      ". Supported: glmmTMB, lme4::merMod, mgcv::gam, stats::glm, stats::lm."
    )
  }

  # ---- Tidy outputs; NA-safe matrix construction ----------------------------
  predictions <- data.frame(
    site = newdata[[site_col]],
    species = newdata[[species_col]],
    pred = as.numeric(pred),
    check.names = FALSE
  )
  names(predictions)[1:2] <- c(site_col, species_col)

  # Drop rows with missing IDs (prevents NA indices during matrix fill)
  keep <- stats::complete.cases(predictions[, c(site_col, species_col)])
  if (!all(keep)) {
    warning(sum(!keep), " prediction rows dropped due to missing site/species IDs.")
    predictions <- predictions[keep, , drop = FALSE]
  }

  # Construct wide matrix (sites x species)
  site_levels <- sort(unique(as.character(predictions[[site_col]])))
  species_levels <- sort(unique(as.character(predictions[[species_col]])))
  prediction_matrix <- matrix(
    NA_real_,
    nrow = length(site_levels),
    ncol = length(species_levels),
    dimnames = list(site_levels, species_levels)
  )
  idx <- cbind(
    match(as.character(predictions[[site_col]]), site_levels),
    match(as.character(predictions[[species_col]]), species_levels)
  )
  prediction_matrix[idx] <- predictions$pred

  # ---- Return ---------------------------------------------------------------
  list(
    newdata           = newdata,
    predictions       = predictions[, c(site_col, species_col, "pred")],
    prediction_matrix = prediction_matrix
  )
}
