#' Generate Spatial Grid
#'
#' This function generates a spatial grid over a geographic extent and assigns grid IDs to points.
#' It also summarizes specified columns within each grid cell.
#'
#' @param data A data frame containing point data with x and y coordinates.
#' @param x_col Character. Column name for x-coordinates (default: "x").
#' @param y_col Character. Column name for y-coordinates (default: "y").
#' @param grid_size Numeric. Size of the grid cells in degrees (default: 0.5).
#' @param sum_col_range Numeric vector. Range of columns to summarize within each grid cell.
#'
#' @return A list containing:
#'   - `grid`: Raster object of the generated grid.
#'   - `grid_sf`: sf object of the grid as polygons.
#'   - `grid_sp`: Data frame summarizing grid cell contents.
#'
#' @import sf
#' @import terra
#' @import dplyr
#' @export
#'
#' @examples
#' # Example usage:
#' data <- data.frame(x = runif(100, -10, 10), y = runif(100, -10, 10), species1 = rpois(100, 5), species2 = rpois(100, 3))
#' grid_result <- generate_grid(data, x_col = "x", y_col = "y", grid_size = 1, sum_col_range = 3:4)
#' print(grid_result$grid_sp)
#' plot(grid_result$grid)

generate_grid <- function(data,
                          x_col = "x",       # Column name for x-coordinates
                          y_col = "y",       # Column name for y-coordinates
                          grid_size = 0.5,    # Grid size in degrees
                          sum_col_range = NULL,  # Columns to summarize (for wide format)
                          species_col = NULL,    # Column name for species (for long format)
                          crs_epsg = 4326       # Set projection *For now only geographic '4326'
) {

  # Ensure required packages are loaded
  required_packages <- c("sf", "terra", "dplyr", "tidyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required but not installed.")
  }

  # Validate inputs
  if (!x_col %in% names(data) || !y_col %in% names(data)) {
    stop("Specified x or y columns do not exist in the input data frame.")
  }
  if (is.null(species_col) && is.null(sum_col_range)) {
    stop("Please specify either `species_col` for long format or `sum_col_range` for wide format.")
  }

  # Convert data to an sf object with WGS84 CRS
  points_sf <- st_as_sf(data, coords = c(x_col, y_col), crs = crs_epsg)

  # Determine and align the geographic extent to the grid size
  bounds <- st_bbox(points_sf)
  x_min <- floor(bounds["xmin"] / grid_size) * grid_size
  x_max <- ceiling(bounds["xmax"] / grid_size) * grid_size
  y_min <- floor(bounds["ymin"] / grid_size) * grid_size
  y_max <- ceiling(bounds["ymax"] / grid_size) * grid_size

  # Extend the extent for a buffer zone of 2 grid cells
  x_min <- x_min - 2 * grid_size
  x_max <- x_max + 2 * grid_size
  y_min <- y_min - 2 * grid_size
  y_max <- y_max + 2 * grid_size

  # Create a raster grid with unique grid IDs
  grid <- terra::rast(
    xmin = x_min, xmax = x_max,
    ymin = y_min, ymax = y_max,
    resolution = grid_size, crs = paste0("EPSG:", crs_epsg)
  )
  grid[] <- 1:terra::ncell(grid)  # Assign unique IDs
  names(grid) <- "grid_id"

  # Extract grid IDs for points
  points <- sf::st_coordinates(points_sf)  # Get coordinates as matrix
  cell_ids <- terra::extract(grid, points)$grid_id
  data$grid_id <- as.character(cell_ids)  # Add grid IDs to data

  # Convert raster to sf polygons for visualization
  grid_sf <- grid %>%
    terra::as.polygons(dissolve = FALSE) %>%
    sf::st_as_sf() %>%
    dplyr::rename(grid_id = grid_id)

  # Generate grid cell centers
  centers <- terra::xyFromCell(grid, 1:terra::ncell(grid)) %>%
    as.data.frame() %>%
    dplyr::rename(x = x, y = y)
  centers$grid_id <- as.character(1:terra::ncell(grid))

  # Summarize data within each grid cell
  if (!is.null(species_col)) {
    # Long format: Summarize species occurrences
    grid_sp <- data %>%
      dplyr::group_by(grid_id) %>%
      dplyr::summarize(obs_sum = n_distinct(!!sym(species_col)), .groups = "drop") %>%
      dplyr::left_join(centers, by = "grid_id")

    # Create wide format for species occurrences
    grid_sp <- data %>%
      dplyr::group_by(grid_id, !!sym(species_col)) %>%
      dplyr::summarize(obs_count = n(), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = !!sym(species_col), values_from = obs_count, values_fill = 0) %>%
      dplyr::left_join(centers, by = "grid_id") %>%
      dplyr::mutate(obs_sum = rowSums(dplyr::select(., -c(grid_id, x, y)), na.rm = TRUE)) %>%
      dplyr::select(grid_id, x, y, obs_sum, everything())

    # Create the grid_obs table
    grid_obs <- data %>%
      dplyr::group_by(grid_id, !!sym(species_col)) %>%
      dplyr::summarize(obs_count = n(), .groups = "drop") %>%
      dplyr::left_join(centers, by = "grid_id") %>%
      dplyr::mutate(obs_sum = grid_sp$obs_sum[match(grid_id, grid_sp$grid_id)]) %>%
      dplyr::select(grid_id, x, y, species = !!sym(species_col), obs_sum, obs_count)
  } else if (!is.null(sum_col_range)) {
    # Wide format: Summarize species columns
    sum_cols <- names(data)[sum_col_range]
    grid_sp <- data %>%
      dplyr::group_by(grid_id) %>%
      dplyr::summarize(across(all_of(sum_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
      dplyr::left_join(centers, by = "grid_id")

    # Calculate total species sum
    grid_sp <- grid_sp %>%
      dplyr::mutate(obs_sum = rowSums(dplyr::select(., all_of(sum_cols)), na.rm = TRUE)) %>%
      dplyr::select(grid_id, x, y, obs_sum, everything())

    # Create the grid_obs table
    grid_obs <- data %>%
      tidyr::pivot_longer(cols = all_of(sum_col_range), names_to = "species", values_to = "presence") %>%
      dplyr::filter(presence > 0) %>%
      dplyr::group_by(grid_id, species) %>%
      dplyr::summarize(obs_count = n(), .groups = "drop") %>%
      dplyr::left_join(centers, by = "grid_id") %>%
      dplyr::mutate(obs_sum = grid_sp$obs_sum[match(grid_id, grid_sp$grid_id)]) %>%
      dplyr::select(grid_id, x, y, species, obs_sum, obs_count)
  }

  # Convert summarized data to a data frame
  grid_sp <- as.data.frame(grid_sp)
  grid_obs <- as.data.frame(grid_obs)
  row.names(grid_sp) <- grid_sp$grid_id

  # Return grid, grid polygons, summarized data, and site observations
  return(list(grid = grid, grid_sf = grid_sf, grid_sp = grid_sp, grid_obs = grid_obs))
}


