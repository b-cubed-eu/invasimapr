#' Format Data Frame
#'
#' This function formats biodiversity data into long or wide format for further analysis.
#'
#' @param data Input data frame containing biodiversity data.
#' @param format Character specifying the format of the data. One of:
#'   - "long": Data is in long format with `species_col` and `value_col`.
#'   - "wide": Data is in wide format with species columns starting with `sp_`.
#'   If NULL, the function will attempt to detect the format automatically.
#' @param x_col Column name for x-coordinates (longitude). If NULL, will be detected.
#' @param y_col Column name for y-coordinates (latitude). If NULL, will be detected.
#' @param site_id_col Column name for site identifiers. If NULL, will be created.
#' @param species_col Column name for species names (required for long format). If NULL, will be detected.
#' @param value_col Column name for species values (e.g., presence/absence or abundance; required for long format). If NULL, will be detected.
#' @param sp_col_range Numeric range specifying species columns in wide format. If NULL, will be detected.
#' @param extra_cols Character vector of additional columns to retain in the output (default: common metadata columns).
#'
#' @return A list with the following elements:
#'   - `site_obs`: A data frame of site-level observations in long format.
#'   - `site_sp`: A data frame of site-level data in wide format.
#'
#' @export
#'
#' @examples
#' # Example for long format data
#' long_data <- data.frame(
#'   x = c(10, 10, 20),
#'   y = c(50, 50, 60),
#'   species = c("sp1", "sp2", "sp1"),
#'   pa = c(1, 1, 0),
#'   coordinateUncertaintyInMeters = c(100, 100, 200)
#' )
#' result <- format_df(data = long_data, format = "long")
#' print(result$site_obs)
#' print(result$site_sp)
#'
#' # Example for wide format data
#' wide_data <- data.frame(
#'   x = c(10, 10, 20),
#'   y = c(50, 50, 60),
#'   species1 = c(1, 1, 0),
#'   species2 = c(2, 2, 1),
#'   recordedBy = c("A", "B", "A"),
#'   eventDate = c("2024-01-01", "2024-01-02", "2024-01-03")
#' )
#' w_result <- format_df(data = wide_data, format = "wide", sp_col_range = 3:4)
#' print(w_result$site_sp)

format_df <- function(data,
                      format = NULL, # Options: "long" or "wide"
                      x_col = NULL, # If not set looks for c("x", "lon", "longitude", "x_coord", "decimalLongitude")
                      y_col = NULL, # ... looks for c("y", "lat", "latitude", "y_coord", "decimalLatitude")
                      site_id_col = NULL, # ... looks for c("site_id", "grid_id", "site", "sample", "id", "plot")
                      species_col = NULL, # ... looks for c("species", "sp_name", "scientific", "verbatimScientificName")
                      value_col = NULL, # ... looks for c("pa", "presence", "abund", "abundance", "count", "individualCount","sum_cnt")
                      sp_col_range = NULL,
                      extra_cols = c("coordinateUncertaintyInMeters",
                                     "recordedBy", "eventDate",
                                     "date", "day", "month", "year")) {
  # Ensure required packages are loaded
  stopifnot(requireNamespace("dplyr", quietly = TRUE),
            requireNamespace("data.table", quietly = TRUE))

  # Default column name mappings
  alt_x_cols <- c("x", "lon", "longitude", "x_coord", "decimalLongitude")
  alt_y_cols <- c("y", "lat", "latitude", "y_coord", "decimalLatitude")
  alt_site_id_cols <- c("site_id","grid_id", "site", "sample", "id", "plot")
  alt_species_cols <- c("species", "sp_name", "scientific", "verbatimScientificName")
  alt_value_cols <- c("pa", "presence", "abund", "abundance", "count", "individualCount","sum_cnt")

  # Helper function to find a column
  find_col <- function(cols, data_names) {
    match <- intersect(tolower(data_names), tolower(cols))
    if (length(match) == 0) return(NULL)
    names(data)[tolower(names(data)) %in% match]
  }

  # Assign columns if not provided
  x_col <- x_col %||% find_col(alt_x_cols, names(data))
  y_col <- y_col %||% find_col(alt_y_cols, names(data))
  site_id_col <- site_id_col %||% find_col(alt_site_id_cols, names(data))
  species_col <- species_col %||% find_col(alt_species_cols, names(data))
  value_col <- value_col %||% find_col(alt_value_cols, names(data))

  # Validate essential columns
  stopifnot(!is.null(x_col), !is.null(y_col))

  # Handle missing `site_id_col`
  if (is.null(site_id_col)) {
    data <- data %>%
      dplyr::group_by(dplyr::across(all_of(c(x_col, y_col)))) %>%
      dplyr::mutate(site_id = paste0("site_", dplyr::cur_group_id())) %>%
      dplyr::ungroup()
    site_id_col <- "site_id"
  }

  # Rename `coordinateUncertaintyInMeters` to `coordUnM` if it exists
  if ("coordinateUncertaintyInMeters" %in% names(data)) {
    data <- data %>%
      dplyr::rename(coordUnM = coordinateUncertaintyInMeters)
    extra_cols <- gsub("coordinateUncertaintyInMeters", "coordUnM", extra_cols)
  }

  # Detect format if not specified
  format <- format %||% ifelse(!is.null(species_col) && !is.null(value_col), "long",
                               ifelse(any(grepl("^sp_", names(data), ignore.case = TRUE)), "wide", NULL))
  stopifnot(!is.null(format))

  # Process long format
  if (format == "long") {
    stopifnot(!is.null(species_col))

    # Rename and select relevant columns
    if ("species" %in% names(data) && species_col != "species") {
      data <- data %>%
        dplyr::rename(species_temp = "species")
    }

    data <- data %>%
      dplyr::rename(site_id = all_of(site_id_col),
                    x = all_of(x_col),
                    y = all_of(y_col),
                    species = all_of(species_col))

    # Handle `value` column
    data <- data %>%
      dplyr::mutate(value = ifelse(is.null(value_col), 1, as.numeric(data[[value_col]]))) %>%
      dplyr::filter(species != "" & !is.na(species))

    site_obs <- data %>%
      dplyr::select(site_id, x, y, species, value, dplyr::any_of(extra_cols))

    # Use data.table for faster pivoting
    dt <- data.table::as.data.table(site_obs)

    # Summarize duplicates before pivoting
    dt <- dt[, .(value = sum(value)), by = .(site_id, x, y, species)]

    # Convert to wide format
    site_sp <- data.table::dcast(dt, site_id + x + y ~ species, value.var = "value", fill = 0)
    site_sp <- as.data.frame(site_sp)
    row.names(site_sp) <- site_sp$site_id

    return(list(site_obs = site_obs, site_sp = site_sp))
  }

  # Process wide format
  if (format == "wide") {
    sp_cols <- if (!is.null(sp_col_range)) names(data)[sp_col_range] else grep("^sp_", names(data), value = TRUE)
    sp_cols <- setdiff(sp_cols, c(site_id_col, x_col, y_col, extra_cols))
    stopifnot(length(sp_cols) > 0)

    # Rename and select relevant columns
    data <- data %>%
      dplyr::rename(site_id = all_of(site_id_col),
                    x = all_of(x_col),
                    y = all_of(y_col))

    # Summarize species columns and group by extra columns
    site_sp <- data %>%
      dplyr::group_by(site_id, x, y, dplyr::across(dplyr::any_of(extra_cols))) %>%
      dplyr::summarize(dplyr::across(all_of(sp_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
    site_sp <- as.data.frame(site_sp)
    row.names(site_sp) <- site_sp$site_id

    return(list(site_sp = site_sp))
  }

  stop("Invalid format specified. Use 'long' or 'wide'.")
}
