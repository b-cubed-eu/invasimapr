#' Scrape and Analyze Wikipedia & Trait Data for a Species
#'
#' Given a binomial species name, this function retrieves optional metadata
#' from Wikipedia (taxonomic summary, taxonomy, image, color palette)
#' and joins relevant plant/trait data from a TRY-style or user-provided trait table.
#' Fuzzy matching is used for both TRY and local tables to handle minor spelling or naming mismatches.
#'
#' @param species Character. Species name (binomial, e.g. "Acacia karroo").
#' @param remove_bg Logical. Remove green/white backgrounds from Wikipedia image? (default: TRUE)
#' @param do_palette,do_taxonomy,do_summary,do_image Logical. Control which metadata to scrape (default: TRUE for all).
#' @param bg_thresh Integer. Brightness threshold for white background removal (default: 80).
#' @param green_delta Integer. How much greener is "green" than R/B? (default: 20).
#' @param n_palette Integer. Number of colors to extract for palette (default: 5).
#' @param preview Logical. Show image after processing? (default: TRUE)
#' @param save_folder Character or NULL. If non-NULL, will save processed PNG image here.
#' @param use_try Logical. If TRUE, join plant traits using a TRY-format database/table (default: FALSE).
#' @param try_data Character (path) or data.frame. Path to TRY file, or data frame containing trait data.
#' @param trait_species_col Name of species column in TRY trait table (default: "AccSpeciesName").
#' @param local_trait_df Optional. Data.frame of local trait data (can be any species-trait table).
#' @param local_species_col Name of species column in local trait table (default: "species").
#' @param max_dist Numeric. Maximum distance for fuzzy join (Levenshtein/Jaro-Winkler; default: 1).
#'
#' @return A tibble (one row) with columns: species, [optional metadata], and all trait columns found.
#' @details
#' \itemize{
#'   \item For TRY tables, `TraitName` is used for wide trait columns. For local tables, all columns except the species column are returned.
#'   \item Fuzzy matching is used to allow for spelling or formatting mismatches.
#'   \item Image-based color palette extraction uses simple k-means clustering; backgrounds can be removed using a color threshold.
#'   \item Requires: dplyr, purrr, tibble, optionally fuzzyjoin, rvest, httr, stringr, jsonlite, magick, abind.
#'   \item You can control which metadata are scraped for speed.
#' }
#'
#' @importFrom rlang .data
#' @importFrom dplyr first
#'
#' @examples
#' \dontrun{
#' # Example using TRY table:
#' get_trait_data("Acacia karroo",
#'                use_try = TRUE,
#'                try_data = try_traits,
#'                trait_species_col = "SpeciesName")
#'
#' # Example using local trait table:
#' get_trait_data("Acraea horta", local_trait_df = traits, local_species_col = "species")
#'
#' # Scrape only metadata (no traits):
#' get_trait_data("Acacia karroo", use_try = FALSE)
#' }
#' @export
get_trait_data <- function(
    species,
    remove_bg = FALSE,
    do_palette = TRUE,
    do_taxonomy = TRUE,
    do_summary = TRUE,
    do_image = TRUE,
    bg_thresh = 80,
    green_delta = 20,
    n_palette = 5,
    preview = FALSE,
    save_folder = NULL,
    use_try = FALSE,
    try_data = NULL,
    trait_species_col = "AccSpeciesName",
    local_trait_df = NULL,
    local_species_col = "species",
    max_dist = 1) {
  requireNamespace("dplyr")
  requireNamespace("purrr")
  requireNamespace("tibble")
  # ensure no connection leakage
  on.exit(closeAllConnections(), add = TRUE)

  # ---- Helper: Clean string encoding ----
  clean_col <- function(x) iconv(as.character(x), from = "", to = "UTF-8", sub = "byte")

  # ---- Wikipedia summary
  # get_wikipedia_summary <- function(species) {
  #   page_title <- gsub(" ", "_", species)
  #   url <- paste0("https://en.wikipedia.org/api/rest_v1/page/summary/", URLencode(page_title, reserved = TRUE))
  #   resp <- tryCatch(httr::GET(url), error = function(e) NULL)
  #   if (is.null(resp) || resp$status_code != 200) return(NA)
  #   cont <- jsonlite::fromJSON(rawToChar(resp$content))
  #   if (!is.null(cont$extract)) return(cont$extract)
  #   NA
  # }
  get_wikipedia_summary <- function(species) {
    page_title <- gsub(" ", "_", species)
    url <- paste0("https://en.wikipedia.org/api/rest_v1/page/summary/", URLencode(page_title, reserved = TRUE))
    resp <- tryCatch(httr::GET(url), error = function(e) NULL)
    if (is.null(resp) || resp$status_code != 200) {
      return(NA)
    }
    cont <- tryCatch(jsonlite::fromJSON(rawToChar(resp$content)), error = function(e) NA)
    # Fix: check if cont is a list and has $extract
    if (is.list(cont) && !is.null(cont$extract)) {
      return(cont$extract)
    } else {
      return(NA)
    }
  }
  # ---- Wikipedia taxonomy via infobox table
  get_wikipedia_taxonomy <- function(species) {
    url <- paste0("https://en.wikipedia.org/wiki/", stringr::str_replace_all(species, " ", "_"))
    page <- tryCatch(rvest::read_html(url), error = function(e) NULL)
    if (is.null(page)) {
      return(setNames(rep(NA_character_, 5), c("Kingdom", "Phylum", "Class", "Order", "Family")))
    }
    infobox <- page %>% rvest::html_node("table.infobox")
    if (inherits(infobox, "xml_missing") || is.null(infobox)) {
      return(setNames(rep(NA_character_, 5), c("Kingdom", "Phylum", "Class", "Order", "Family")))
    }
    tbl <- infobox %>%
      rvest::html_table(fill = TRUE) %>%
      tibble::as_tibble(.name_repair = "minimal") %>%
      purrr::set_names(c("trait", "value")) %>%
      dplyr::filter(trait != "" & value != "")
    tbl <- tbl %>% dplyr::mutate(trait = stringr::str_remove(trait, ":$"), value = stringr::str_squish(value))
    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family")
    vals <- setNames(rep(NA_character_, length(ranks)), ranks)
    for (rk in ranks) {
      v <- tbl$value[tolower(tbl$trait) == tolower(rk)]
      if (length(v) > 0) vals[rk] <- v[1]
    }
    vals
  }
  # ---- Wikipedia image url
  get_wikipedia_image <- function(species) {
    url <- paste0("https://en.wikipedia.org/wiki/", gsub(" ", "_", species))
    html <- tryCatch(rvest::read_html(url), error = function(e) NULL)
    if (is.null(html)) {
      return(NA)
    }
    img_url <- html %>%
      rvest::html_node(".infobox img") %>%
      rvest::html_attr("src")
    if (is.na(img_url) || is.null(img_url)) {
      return(NA)
    }
    if (!startsWith(img_url, "http")) img_url <- paste0("https:", img_url)
    img_url
  }
  # ---- Background removal (simple color threshold)
  remove_background_color <- function(image_url, bg_thresh = 80, green_delta = 20, out_file = NULL) {
    img <- tryCatch(magick::image_read(image_url), error = function(e) NULL)
    if (is.null(img)) {
      return(NA)
    }
    img <- magick::image_scale(img, "200")
    img_data <- magick::image_data(img, channels = "rgb")
    r <- as.integer(img_data[1, , ])
    g <- as.integer(img_data[2, , ])
    b <- as.integer(img_data[3, , ])
    mask <- !(
      (g > r + green_delta & g > b + green_delta & g > bg_thresh) |
        (r > bg_thresh & g > bg_thresh & b > bg_thresh)
    )
    alpha <- as.raw(ifelse(mask, 255, 0))
    img_data_alpha <- abind::abind(img_data, array(alpha, dim = dim(img_data)[2:3]), along = 1)
    img_rgba <- magick::image_read(img_data_alpha)
    if (!is.null(out_file)) magick::image_write(img_rgba, out_file)
    img_rgba
  }
  # ---- Palette extraction
  get_image_palette_local <- function(img, n = 5) {
    img <- magick::image_scale(img, "100")
    img_data <- magick::image_data(img, channels = "rgba")
    alpha <- as.integer(img_data[4, , ])
    keep <- which(alpha > 0)
    if (length(keep) < n) {
      return(NA)
    }
    r <- as.integer(img_data[1, , ])[keep]
    g <- as.integer(img_data[2, , ])[keep]
    b <- as.integer(img_data[3, , ])[keep]
    df <- data.frame(r = r, g = g, b = b)
    if (nrow(unique(df)) < n) {
      return(NA)
    }
    km <- kmeans(df, centers = n)
    rgb <- round(km$centers)
    rgb_hex <- apply(rgb, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    paste(rgb_hex, collapse = ", ")
  }

  # ---- TRY trait table (wide, trait names as columns) ----
  try_traits <- NULL
  if (use_try && !is.null(try_data)) {
    if (missing(trait_species_col) || is.null(trait_species_col) || trait_species_col == "" || !(trait_species_col %in% names(try_data))) {
      stop("trait_species_col must be a valid column name in try_data.")
    }
    trytab <- try_data %>%
      dplyr::filter(!is.na(.data[[trait_species_col]]) & .data[[trait_species_col]] != "")
    trait_cols <- intersect(c(trait_species_col, "TraitName", "OrigValueStr"), names(trytab))
    trytab <- as.data.frame(trytab)[, trait_cols, drop = FALSE]
    suppressWarnings({
      trait_long <- fuzzyjoin::stringdist_left_join(
        tibble::tibble(species = species),
        trytab,
        by = c("species" = trait_species_col),
        max_dist = max_dist,
        method = "jw"
      )
    })
    if (nrow(trait_long) > 0) {
      trait_wide <- trait_long %>%
        dplyr::filter(!is.na(OrigValueStr)) %>%
        dplyr::mutate(
          OrigValueStr = iconv(OrigValueStr, from = "", to = "UTF-8", sub = NA),
          OrigValueStr = suppressWarnings(as.numeric(OrigValueStr))
        ) %>%
        dplyr::filter(!is.na(TraitName), TraitName != "", !is.na(OrigValueStr)) %>%
        dplyr::group_by(TraitName) %>%
        dplyr::summarise(value = first(OrigValueStr), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = TraitName, values_from = value)
      if (nrow(trait_wide) > 0) try_traits <- trait_wide
    }
  }

  # ---- Local trait table (any columns, wide) ----
  local_traits <- NULL
  if (!is.null(local_trait_df) && is.data.frame(local_trait_df)) {
    if (!(local_species_col %in% names(local_trait_df))) {
      stop("local_species_col must be a valid column name in local_trait_df.")
    }
    localtab <- local_trait_df %>%
      dplyr::filter(!is.na(.data[[local_species_col]]) & .data[[local_species_col]] != "")
    trait_cols <- setdiff(names(localtab), local_species_col)
    if (length(trait_cols) == 0) stop("No trait columns found in local_trait_df.")
    suppressWarnings({
      trait_match <- fuzzyjoin::stringdist_left_join(
        tibble::tibble(species = species),
        localtab,
        by = setNames(local_species_col, "species"),
        max_dist = max_dist,
        ignore_case = TRUE,
        distance_col = "distance" # name of the Levenshtein-distance column
      )
    })
    if (nrow(trait_match) > 0) {
      trait_match <- trait_match %>% dplyr::select(dplyr::all_of(trait_cols))
      local_traits <- trait_match[1, , drop = FALSE]
    }
  }

  # ---- Metadata scrape: conditional on user opts ----
  summary <- if (do_summary) get_wikipedia_summary(species) else NA
  taxonomy <- if (do_taxonomy) get_wikipedia_taxonomy(species) else setNames(rep(NA, 5), c("Kingdom", "Phylum", "Class", "Order", "Family"))
  img_url <- if (do_image) get_wikipedia_image(species) else NA
  img_rgba <- if (do_image && !is.na(img_url) && remove_bg) {
    remove_background_color(img_url, bg_thresh, green_delta)
  } else if (do_image && !is.na(img_url)) {
    magick::image_read(img_url)
  } else {
    NULL
  }
  if (preview && inherits(img_rgba, "magick-image")) print(img_rgba)
  if (!is.null(save_folder) && inherits(img_rgba, "magick-image")) {
    if (!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    fname <- paste0(gsub(" ", "_", species), "_bg", if (remove_bg) "_removed" else "", ".png")
    outfile <- file.path(save_folder, fname)
    magick::image_write(img_rgba, outfile)
  }
  palette <- if (do_palette && inherits(img_rgba, "magick-image")) get_image_palette_local(img_rgba, n = n_palette) else NA

  # ---- Build final output row ----
  output <- tibble::tibble(
    species = species,
    summary = summary,
    Kingdom = taxonomy["Kingdom"],
    Phylum = taxonomy["Phylum"],
    Class = taxonomy["Class"],
    Order = taxonomy["Order"],
    Family = taxonomy["Family"],
    img_url = img_url,
    palette = palette
  )
  if (!is.null(try_traits)) output <- dplyr::bind_cols(output, try_traits)
  if (!is.null(local_traits)) output <- dplyr::bind_cols(output, local_traits)
  output
}
