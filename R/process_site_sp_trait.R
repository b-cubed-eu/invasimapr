#' Process site-level species–trait data
#'
#' Harmonises raw TRY trait records with site–by–species occurrence tables
#' and returns matched, analysis-ready objects for functional‐diversity work.
#' The workflow:
#' \enumerate{
#'   \item Imports or validates a TRY‐formatted data frame (or file).
#'   \item Fuzzy-matches the species in \code{grid_obs} to accepted TRY names
#'         (Jaro–Winkler distance \code{<= 0.05} via
#'         \code{fuzzyjoin::stringdist_inner_join()}).
#'   \item Aggregates multiple numeric observations per species–trait pair to a
#'         single value, then reshapes to a wide species × trait matrix.
#'   \item Filters the site × species matrix (\code{grid_sp}) to retain only
#'         species with trait information.
#'   \item Returns a named list containing the filtered site and trait tables,
#'         spatial coordinates, a trait-name lookup, and (optionally) diagnostic
#'         objects.
#' }
#'
#' @param try_data \code{character} or \code{data.frame}.  Either the path to a
#'   tab-delimited TRY export or an in-memory data frame with the columns
#'   \code{AccSpeciesName}, \code{TraitID}, \code{TraitName}, \code{OrigValueStr}.
#' @param grid_obs \code{data.frame}.  Long table (\code{grid_id}, \code{species})
#'   listing every species observation per grid cell.
#' @param grid_sp \code{data.frame}.  Wide site × species matrix with coordinates;
#'   must contain \code{grid_id}, \code{x}, \code{y} and the columns named in
#'   \code{sp_cols}.
#' @param sp_cols \code{character}.  The species columns in \code{grid_sp} that
#'   should be checked for matching trait data.
#' @param appendix \code{logical(1)}.  If \code{TRUE}, also returns
#'   \code{sbt} (the unfiltered species × trait matrix) in addition to the core
#'   outputs.  Defaults to \code{FALSE}.
#'
#' @return A list with components
#' \describe{
#'   \item{\code{sp_traits}}{Matched species × trait data frame (wide).}
#'   \item{\code{site_sp}}{Filtered site × species matrix (wide).}
#'   \item{\code{site_xy}}{Data frame of \code{x}, \code{y} coordinates for the
#'         sites in \code{site_sp}.}
#'   \item{\code{traitname}}{Lookup table linking synthetic trait IDs
#'         (\code{trait_1}, …) to original TRY trait names.}
#'   \item{\code{species_mapping}}{Table mapping original species names to
#'         matched TRY names.}
#'   \item{\code{sbt}}{(Only if \code{appendix = TRUE}) the full
#'         species × trait matrix before filtering.}
#' }
#'
#' @details
#' Non-numeric \code{OrigValueStr} entries are coerced to \code{NA}.
#' Where multiple trait measurements exist for a species, the first numeric
#' record is kept.  Users can adjust this behaviour manually after inspection of
#' the returned \code{sbt} object.
#'
#' @seealso \code{\link[fuzzyjoin]{stringdist_inner_join}},
#'   \code{\link[rtry]{rtry_import}}
#'
#' @importFrom rtry rtry_import
#' @importFrom fuzzyjoin stringdist_inner_join
#' @importFrom dplyr select distinct group_by summarise filter mutate across
#'   arrange inner_join pull
#' @importFrom tidyr drop_na pivot_wider pivot_longer
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
#'
#' @examples
#' \dontrun{
#' # --- toy data ------------------------------------------------------------
#' toy_try <- data.frame(
#'   AccSpeciesName = c("Acacia karroo", "Panicum maximum"),
#'   TraitID        = c(23, 23),
#'   TraitName      = c("Specific leaf area", "Specific leaf area"),
#'   OrigValueStr   = c("14.2", "28.1")
#' )
#'
#' grid_obs <- data.frame(
#'   grid_id = c(1, 1, 2, 2),
#'   species = c("Acacia karroo", "Panicum maximum",
#'               "Acacia karoo",  "Panicum maximum")
#' )
#'
#' grid_sp <- data.frame(
#'   grid_id = c(1, 2),
#'   x = c(30.1, 30.2),
#'   y = c(-24.8, -24.7),
#'   `Acacia karroo`   = c(1, 0),
#'   `Panicum maximum` = c(1, 1)
#' )
#'
#' res <- process_site_sp_trait(
#'   try_data = toy_try,
#'   grid_obs = grid_obs,
#'   grid_sp  = grid_sp,
#'   sp_cols  = c("Acacia karroo", "Panicum maximum")
#' )
#' str(res)
#' }
process_site_sp_trait = function(try_data,
                                  grid_obs, # Unfiltered site-by-observations table (long format)
                                  grid_sp, # Unfiltered site-by-species table
                                  sp_cols,
                                  appendix = FALSE) {
  # Ensure required packages are loaded
  required_packages = c("rtry", "cli", "dplyr", "tidyr", "tibble", "fuzzyjoin")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required but not installed.")
  }

  # Step 1: Process TRY data
  if ("character" %in% class(try_data)) {
    trydata = rtry::rtry_import(
      input = try_data,
      separator = "\t",
      encoding = "UTF-8",
      quote = "",
      showOverview = TRUE
    )
  } else if ("data.frame" %in% class(try_data)) {
    if (any(!c("AccSpeciesName", "TraitID", "TraitName", "OrigValueStr") %in% colnames(try_data))) {
      requiredcol = c("AccSpeciesName", "TraitID", "TraitName", "OrigValueStr")
      missingcol = requiredcol[!c("AccSpeciesName", "TraitID", "TraitName", "OrigValueStr") %in% colnames(try_data)]
      cli::cli_abort(c("{missingcol} is/are not in the {.var try_data} column "))
    }
    trydata = try_data
  } else {
    cli::cli_abort(c("{.var try_data} is not a file path or dataframe"))
  }

  # Step 2: Extract unique species from grid_obs and perform fuzzy matching
  species_list = sort(unique(grid_obs$species))
  merged_data = fuzzyjoin::stringdist_inner_join(
    grid_obs %>% dplyr::select(species),
    trydata %>% dplyr::select(AccSpeciesName),
    by = c("species" = "AccSpeciesName"),
    method = "jw",
    max_dist = 0.05,
    distance_col = "match_level"
  )
  species_mapping = merged_data %>%
    dplyr::select(species, AccSpeciesName) %>%
    dplyr::distinct()

  # Step 3: Process trait data
  SpeciesbyTrait = trydata %>%
    tidyr::drop_na(TraitID) %>%
    dplyr::select(AccSpeciesName, TraitID, OrigValueStr) %>%
    dplyr::mutate(
      OrigValueStr = iconv(OrigValueStr, from = "UTF-8", to = "UTF-8", sub = NA),
      OrigValueStr = suppressWarnings(as.numeric(OrigValueStr))
    ) %>%
    dplyr::group_by(AccSpeciesName, TraitID) %>%
    dplyr::summarise(across(OrigValueStr, first), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = TraitID, values_from = OrigValueStr) %>%
    dplyr::filter(AccSpeciesName %in% species_mapping$AccSpeciesName) %>%
    tibble::column_to_rownames(var = "AccSpeciesName") %>%
    dplyr::select_if(~ !all(is.na(.)))

  # Step 4: Extract trait names
  trait = colnames(SpeciesbyTrait)
  traitname = trydata %>%
    tidyr::drop_na(TraitID) %>%
    dplyr::select(TraitID, TraitName) %>%
    dplyr::group_by(TraitID) %>%
    dplyr::summarise(across(TraitName, first), .groups = "drop") %>%
    tibble::column_to_rownames("TraitID")
  traitname = c(traitname[trait, ])
  traitname = data.frame('TraitID' = paste0("trait_", 1:ncol(SpeciesbyTrait)),
                          'TraitName' = traitname)
  names(SpeciesbyTrait) = paste0("trait_", 1:ncol(SpeciesbyTrait))

  # Step 5: Filter traits based on grid_sp and sp_cols
  grid_long = grid_sp %>%
    tidyr::pivot_longer(
      cols = all_of(sp_cols),
      names_to = "species",
      values_to = "value"
    ) %>%
    dplyr::filter(value > 0)

  try_wide = SpeciesbyTrait %>%
    tibble::rownames_to_column("species")

  sp_traits = grid_long %>%
    dplyr::distinct(species) %>%
    dplyr::arrange(species) %>%
    dplyr::inner_join(try_wide, by = "species")

  matched_species = sp_traits$species

  site_sp = grid_long %>%
    dplyr::filter(species %in% matched_species) %>%
    tidyr::pivot_wider(
      names_from = species,
      values_from = value,
      values_fill = 0
    )

  site_xy = grid_sp %>%
    dplyr::filter(grid_id %in% site_sp$grid_id) %>%
    dplyr::select(grid_id, x, y)

  # Return combined output
  output = list(sp_traits = sp_traits, site_sp = site_sp, site_xy = site_xy)

  if (appendix) {
    output = c(output, list("sbt" = SpeciesbyTrait, "traitname" = traitname, "species_mapping" = species_mapping))
  } else {
    output = c(output, list("traitname" = traitname, "species_mapping" = species_mapping))
  }

  return(output)
}
