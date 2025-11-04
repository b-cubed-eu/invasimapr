#' Assemble community matrices for invasion-fitness workflows
#'
#' @description
#' `assemble_matrices()` creates the standard inputs used throughout the
#' invasion-fitness pipeline from either a single “long” table or a set of
#' already-formed tables. It returns aligned objects:
#' - `site_df` (site × x,y),
#' - `env_df` (site × environment; **numeric matrix**),
#' - `comm_res` (site × resident abundances; numeric matrix),
#' - `pa_res` (site × resident presence/absence; integer matrix),
#' - `traits_res` (resident × traits; data frame, keeps mixed types).
#'
#' @param long_df Optional long data frame with at least `site`, `x`, `y`,
#'   `species`, `count`.
#' @param site_df Optional data frame with columns `site`, `x`, `y` (used when
#'   `long_df` is `NULL`).
#' @param env_df Optional site × environment **numeric** data frame/matrix
#'   (row names must be site IDs; used when `long_df` is `NULL`).
#' @param comm_res Optional site × resident numeric matrix/data frame (wide) or
#'   a long table with `site`, `species`, `count` (see `comm_long`).
#' @param traits_res Optional resident × traits data frame (row names = species).
#' @param comm_long How to interpret a separately supplied `comm_res`:
#'   one of `c("auto","long","wide")`; `"auto"` detects a long table and pivots.
#' @param site_col,x_col,y_col,species_col,count_col Column names to use when
#'   building from `long_df`. Defaults: `"site"`, `"x"`, `"y"`, `"species"`, `"count"`.
#' @param env_cols,env_prefix,trait_cols,trait_prefix,drop_empty_sites,
#'   drop_empty_species,return_diversity,make_plots See details in original docs.
#'
#' @importFrom dplyr %>% distinct select arrange group_by summarise mutate left_join relocate
#' @importFrom tidyr pivot_wider
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom stats setNames
#' @export
assemble_matrices = function(
    long_df = NULL,
    site_df = NULL,
    env_df = NULL,
    comm_res = NULL,
    traits_res = NULL,
    comm_long = c("auto", "long", "wide"),
    site_col = "site",
    x_col = "x",
    y_col = "y",
    species_col = "species",
    count_col = "count",
    env_cols = NULL,
    env_prefix = "^env",
    trait_cols = NULL,
    trait_prefix = "^trait",
    drop_empty_sites = TRUE,
    drop_empty_species = TRUE,
    return_diversity = TRUE,
    make_plots = FALSE
) {
  `%||%` = function(a, b) if (!is.null(a)) a else b
  comm_long = match.arg(comm_long)

  # Helper: standardize a separately supplied comm_res (long -> wide if needed)
  .standardize_comm = function(comm_res, site_col, species_col, count_col, comm_long){
    if (is.null(comm_res)) return(NULL)
    if (is.matrix(comm_res) && is.numeric(comm_res)) return(comm_res)
    df = as.data.frame(comm_res)

    treat_as_long = switch(
      comm_long,
      "long" = TRUE,
      "wide" = FALSE,
      "auto" = all(c(site_col, species_col, count_col) %in% names(df))
    )

    if (treat_as_long) {
      if (!all(c(site_col, species_col, count_col) %in% names(df))) {
        stop("`comm_res` declared as long but missing one of: ",
             paste(c(site_col, species_col, count_col), collapse = ", "))
      }
      tmp = df |>
        dplyr::select(dplyr::all_of(c(site_col, species_col, count_col))) |>
        dplyr::group_by(.data[[site_col]], .data[[species_col]]) |>
        dplyr::summarise(count = sum(.data[[count_col]]), .groups = "drop") |>
        tidyr::pivot_wider(names_from = dplyr::all_of(species_col),
                           values_from = "count",
                           values_fill = 0) |>
        dplyr::arrange(.data[[site_col]]) |>
        as.data.frame()
      rownames(tmp) = as.character(tmp[[site_col]])
      tmp[[site_col]] = NULL
      m = as.matrix(tmp); storage.mode(m) = "double"
      return(m)
    }

    # assume wide; allow a 'site' column; drop non-numeric data cols
    if (site_col %in% names(df)) {
      rn = as.character(df[[site_col]])
      num_idx = vapply(df, is.numeric, logical(1))
      num_idx[match(site_col, names(df))] = FALSE
      m = as.matrix(df[, num_idx, drop = FALSE]); rownames(m) = rn
      storage.mode(m) = "double"; return(m)
    } else {
      num_idx = vapply(df, is.numeric, logical(1))
      m = as.matrix(df[, num_idx, drop = FALSE])
      storage.mode(m) = "double"; return(m)
    }
  }

  # --------------------------- 0) Basic checks --------------------------------
  if (is.null(long_df) && is.null(site_df) && is.null(env_df) &&
      is.null(comm_res) && is.null(traits_res)) {
    stop("Provide either `long_df` or any of the individual tables (site_df/env_df/comm_res/traits_res).")
  }

  if (!is.null(long_df)) long_df = as.data.frame(long_df)

  # --------------------------- 1) Build from long_df --------------------------
  if (!is.null(long_df)) {
    need_cols = c(site_col, x_col, y_col, species_col, count_col)
    if (!all(need_cols %in% names(long_df))) {
      stop("`long_df` must contain: ", paste(need_cols, collapse = ", "), ".")
    }

    # 1a) Site coordinates (site_df)
    if (is.null(site_df)) {
      site_df = long_df |>
        dplyr::distinct(.data[[site_col]], .keep_all = TRUE) |>
        dplyr::select(dplyr::all_of(c(site_col, x_col, y_col))) |>
        dplyr::arrange(.data[[site_col]])
      names(site_df) = c("site","x","y")
      rownames(site_df) = as.character(site_df$site)
      site_df$site = as.character(site_df$site)
    }

    # 1b) Environment (env_df)
    if (is.null(env_df)) {
      env_cols = env_cols %||% grep(env_prefix, names(long_df), value = TRUE)
      if (length(env_cols)) {
        env_df = long_df |>
          dplyr::distinct(.data[[site_col]], .keep_all = TRUE) |>
          dplyr::select(dplyr::all_of(c(site_col, env_cols))) |>
          dplyr::arrange(.data[[site_col]])
        rownames(env_df) = as.character(env_df[[site_col]])
        env_df[[site_col]] = NULL
        # keep original behavior: coerce to numeric where possible
        for (j in seq_along(env_df)) {
          if (!is.numeric(env_df[[j]])) {
            suppressWarnings({
              x = as.numeric(env_df[[j]])
              if (sum(!is.na(x)) >= sum(!is.na(env_df[[j]]))) env_df[[j]] = x
            })
          }
        }
      } else env_df = NULL
    }

    # 1c) Community matrix (comm_res) + presence/absence (pa_res)
    if (is.null(comm_res)) {
      tmp = long_df |>
        dplyr::select(dplyr::all_of(c(site_col, species_col, count_col))) |>
        dplyr::group_by(.data[[site_col]], .data[[species_col]]) |>
        dplyr::summarise(count = sum(.data[[count_col]]), .groups = "drop")
      comm_res = tmp |>
        tidyr::pivot_wider(names_from = dplyr::all_of(species_col),
                           values_from = "count",
                           values_fill = 0) |>
        dplyr::arrange(.data[[site_col]]) |>
        as.data.frame()
      rownames(comm_res) = as.character(comm_res[[site_col]])
      comm_res[[site_col]] = NULL
    }

    # 1d) Traits (traits_res)
    if (is.null(traits_res)) {
      trait_cols = trait_cols %||% grep(trait_prefix, names(long_df), value = TRUE)
      if (length(trait_cols)) {
        traits_res = long_df |>
          dplyr::distinct(.data[[species_col]], .keep_all = TRUE) |>
          dplyr::select(dplyr::all_of(c(species_col, trait_cols))) |>
          dplyr::arrange(.data[[species_col]])
        rownames(traits_res) = as.character(traits_res[[species_col]])
        traits_res[[species_col]] = NULL
      } else traits_res = NULL
    }
  }

  # make all character trait columns factors (Gower needs factors, not chars)
  if (!is.null(traits_res)) {
    for (j in seq_along(traits_res)) {
      if (is.character(traits_res[[j]])) traits_res[[j]] = factor(traits_res[[j]])
    }
  }

  # ---- NEW: standardize separate comm_res (long -> wide if needed) ------------
  if (is.null(long_df) && !is.null(comm_res)) {
    comm_res = .standardize_comm(comm_res, site_col, species_col, count_col, comm_long)
  }

  # ---- Standardize separately supplied traits_res (use species_cola as rownames) ----
  if (is.null(long_df) && !is.null(traits_res)) {
    traits_res = as.data.frame(traits_res, stringsAsFactors = FALSE)
    if (species_col %in% names(traits_res)) {
      rn = trimws(as.character(traits_res[[species_col]]))
      traits_res[[species_col]] = NULL
      rownames(traits_res) = rn
    } else if (is.null(rownames(traits_res))) {
      stop("When supplying `traits_res` separately, provide a '", species_col,
           "' column or row names corresponding to species IDs.")
    }
  }

  # ---- NEW: rename coords to x/y when supplied separately ---------------------
  if (is.null(long_df)) {
    if (!is.null(site_df)) {
      if (!"x" %in% names(site_df) && x_col %in% names(site_df))
        names(site_df)[names(site_df) == x_col] = "x"
      if (!"y" %in% names(site_df) && y_col %in% names(site_df))
        names(site_df)[names(site_df) == y_col] = "y"
    }
    if (!is.null(env_df)) {
      if (!"x" %in% names(env_df) && x_col %in% names(env_df))
        names(env_df)[names(env_df) == x_col] = "x"
      if (!"y" %in% names(env_df) && y_col %in% names(env_df))
        names(env_df)[names(env_df) == y_col] = "y"
    }
  }

  # Presence/absence from comm_res (always derivable here)
  pa_res = if (!is.null(comm_res)) {
    as.data.frame((comm_res > 0) * 1L, check.names = FALSE)
  } else NULL

  # ---- Normalize site_df / env_df IDs and coords (ensure 'site' matches rownames) ----
  if (!is.null(site_df)) {
    site_df = as.data.frame(site_df, stringsAsFactors = FALSE)

    # Set row names from site_col when present; otherwise keep existing
    if (site_col %in% names(site_df)) {
      rownames(site_df) = as.character(site_df[[site_col]])
    }

    # Rename coordinates to x/y if user-specified names differ
    if (!"x" %in% names(site_df) && x_col %in% names(site_df)) {
      names(site_df)[names(site_df) == x_col] = "x"
    }
    if (!"y" %in% names(site_df) && y_col %in% names(site_df)) {
      names(site_df)[names(site_df) == y_col] = "y"
    }

    # Ensure canonical 'site' column exists and matches rownames (as character)
    site_df$site = as.character(rownames(site_df))
  }

  if (!is.null(env_df)) {
    env_df = as.data.frame(env_df, stringsAsFactors = FALSE)

    # Set row names from site_col when present; otherwise keep existing
    if (site_col %in% names(env_df)) {
      rownames(env_df) = as.character(env_df[[site_col]])
    }

    # Rename coordinates to x/y if user-specified names differ
    if (!"x" %in% names(env_df) && x_col %in% names(env_df)) {
      names(env_df)[names(env_df) == x_col] = "x"
    }
    if (!"y" %in% names(env_df) && y_col %in% names(env_df)) {
      names(env_df)[names(env_df) == y_col] = "y"
    }

    # Ensure canonical 'site' column exists and matches rownames (as character)
    # env_df$site = as.character(rownames(env_df))
  }

  # --------------------------- 2) Drop empties (opt) --------------------------
  if (!is.null(comm_res)) {
    if (isTRUE(drop_empty_species)) {
      keep_cols = colSums(comm_res, na.rm = TRUE) > 0
      if (any(!keep_cols)) {
        comm_res = comm_res[, keep_cols, drop = FALSE]
        if (!is.null(pa_res)) pa_res = pa_res[, keep_cols, drop = FALSE]
        if (!is.null(traits_res)) {
          traits_res = traits_res[intersect(rownames(traits_res), colnames(comm_res)), , drop = FALSE]
        }
      }
    }
    if (isTRUE(drop_empty_sites)) {
      keep_rows = rowSums(comm_res, na.rm = TRUE) > 0
      if (any(!keep_rows)) {
        comm_res = comm_res[keep_rows, , drop = FALSE]
        if (!is.null(pa_res))  pa_res  = pa_res [keep_rows, , drop = FALSE]
        if (!is.null(site_df)) site_df = site_df[rownames(comm_res), , drop = FALSE]
        if (!is.null(env_df))  env_df  = env_df [rownames(comm_res), , drop = FALSE]
      }
    }
  }

  # --------------------------- 3) Align by common sets ------------------------
  site_lists = list()
  if (!is.null(site_df)) site_lists = c(site_lists, list(rownames(site_df)))
  if (!is.null(env_df))  site_lists = c(site_lists, list(rownames(env_df)))
  if (!is.null(comm_res))site_lists = c(site_lists, list(rownames(comm_res)))
  if (!length(site_lists)) stop("No site-bearing objects to align.")
  common_sites = Reduce(intersect, site_lists)
  if (!length(common_sites)) {
    stop("No common sites among the provided inputs. ",
         "Check that site IDs match across tables (same case/prefix, no extra whitespace).")
  }
  if (!is.null(site_df)) site_df = site_df[common_sites, , drop = FALSE]
  if (!is.null(env_df))  env_df  = env_df [common_sites, , drop = FALSE]
  if (!is.null(comm_res)) {
    comm_res = comm_res[common_sites, , drop = FALSE]
    if (!is.null(pa_res)) pa_res = pa_res[common_sites, , drop = FALSE]
  }

  # Residents/species alignment
  if (!is.null(comm_res) && !is.null(traits_res)) {
    common_residents = intersect(colnames(comm_res), rownames(traits_res))
    if (!length(common_residents)) {
      stop("No overlap between species in `comm_res` and `traits_res`. ",
           "Check species IDs / spelling / case.")
    }
    comm_res   = comm_res[, common_residents, drop = FALSE]
    pa_res     = pa_res[,   common_residents, drop = FALSE]
    traits_res = traits_res[common_residents, , drop = FALSE]
  } else if (!is.null(comm_res)) {
    common_residents = colnames(comm_res)
  } else if (!is.null(traits_res)) {
    common_residents = rownames(traits_res)
  } else {
    common_residents = character(0)
  }

  # --------------------------- 4) Diversity (optional) ------------------------
  diversity = NULL
  if (isTRUE(return_diversity) && !is.null(comm_res)) {
    comm_res_df = comm_res |>
      as.data.frame() |>
      tibble::rownames_to_column("site")

    spp_mat = as.matrix(comm_res_df[, setdiff(names(comm_res_df), "site"), drop = FALSE])

    diversity = tibble::tibble(
      site     = comm_res_df$site,
      spp_rich = rowSums(spp_mat > 0, na.rm = TRUE),
      obs_sum  = rowSums(spp_mat,      na.rm = TRUE),
      obs_mean = rowMeans(spp_mat,     na.rm = TRUE)
    )

    if (!is.null(site_df)) {
      has_site_col = "site" %in% names(site_df)
      coord_tbl = if (has_site_col) {
        site_df[, intersect(c("site","x","y"), names(site_df)), drop = FALSE]
      } else {
        tibble::rownames_to_column(as.data.frame(site_df), "site") |>
          dplyr::select(dplyr::all_of(intersect(c("site","x","y"), names(site_df))))
      }
      if (all(c("site","x","y") %in% names(coord_tbl))) {
        diversity = diversity |>
          dplyr::mutate(site = as.character(site)) |>
          dplyr::left_join(coord_tbl, by = "site") |>
          dplyr::relocate(site, x, y)
      }
    }
  }

  # (Optional) quick plot
  plots = list()
  if (isTRUE(make_plots) && !is.null(diversity) && requireNamespace("ggplot2", quietly = TRUE)) {
    plots$richness = ggplot2::ggplot(diversity, ggplot2::aes(x, y, fill = sqrt(spp_rich))) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(name = "√(Richness)") +
      ggplot2::labs(title = "Site-level species richness", x = "x", y = "y") +
      ggplot2::theme_minimal()
  }

  # site_df$site = row.names(site_df)
  # env_df$site = row.names(env_df)

  # Just before the return()
  if (!is.null(env_df) && "site" %in% names(env_df)) {
    env_df[["site"]] = NULL
  }

  # --------------------------- 5) Return --------------------------------------
  list(
    site_df          = site_df,
    env_df           = env_df,
    comm_res         = comm_res,
    pa_res           = pa_res,
    traits_res       = traits_res,
    diversity        = diversity,
    sites            = rownames(site_df),
    residents        = common_residents,
    n_sites          = nrow(site_df),
    n_env            = if (!is.null(env_df)) ncol(env_df) else 0L,
    n_residents      = length(common_residents),
    n_traits         = if (!is.null(traits_res)) ncol(traits_res) else 0L,
    plots            = plots
  )
}
