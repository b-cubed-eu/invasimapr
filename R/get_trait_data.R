#' Scrape and Analyze Wikipedia & Trait Data for a Species
#'
#' Given a binomial species name, this function retrieves optional metadata
#' from Wikipedia (taxonomic summary, taxonomy, image, and a color palette)
#' and joins relevant plant/trait data from a TRY-style or user-provided trait table.
#' Fuzzy matching is used for both TRY and local tables to handle minor spelling
#' or naming mismatches.
#'
#' When `remove_bg = TRUE`, the infobox image background is removed using the
#' remove.bg API via an internal helper (`remove_bg_and_save()`), the resulting
#' PNG is re-read with **magick**, and the palette is extracted from that
#' processed image. Set the environment variable `REMOVE_BG_API_KEY` to a valid
#' remove.bg API key before calling.
#'
#' @param species Character. Species name (binomial, e.g. `"Acacia karroo"`).
#' @param remove_bg Logical. If `TRUE`, call remove.bg via `remove_bg_and_save()`
#'   to remove the background from the Wikipedia infobox image and use the
#'   processed PNG for preview/palette. Default: `FALSE`.
#' @param do_palette,do_taxonomy,do_summary,do_image Logical. Control which
#'   metadata to scrape. All default to `TRUE`.
#' @param bg_thresh Integer. Deprecated/ignored when `remove_bg = TRUE`.
#'   Kept for signature compatibility. Default: `80`.
#' @param green_delta Integer. Deprecated/ignored when `remove_bg = TRUE`.
#'   Kept for signature compatibility. Default: `20`.
#' @param n_palette Integer. Number of colors to extract for the palette.
#'   Default: `5`.
#' @param preview Logical. Print the processed image (magick) in the console.
#'   Default: `FALSE`.
#' @param save_folder Character or `NULL`. If non-`NULL`, write the PNG used
#'   for palette/preview into this folder. When `remove_bg = TRUE`, the
#'   background-removed PNG is written here; otherwise the original image is
#'   written. If `NULL` and `remove_bg = TRUE`, a temporary directory is used.
#' @param use_try Logical. If `TRUE`, join plant traits using a TRY-format
#'   database/table. Default: `FALSE`.
#' @param try_data Character (path) or `data.frame`. Path to a TRY file, or a
#'   data frame containing TRY-style trait data (must include `trait_species_col`,
#'   `TraitName`, and `OrigValueStr`).
#' @param trait_species_col Name of the species column in the TRY trait table.
#'   Default: `"AccSpeciesName"`.
#' @param local_trait_df Optional. `data.frame` of local trait data (any wide
#'   table). All columns except the species column are returned.
#' @param local_species_col Name of the species column in the local trait table.
#'   Default: `"species"`.
#' @param max_dist Numeric. Maximum distance for fuzzy matching (Jaro–Winkler
#'   via `fuzzyjoin::stringdist_left_join`). Default: `1`.
#'
#' @return A tibble (one row) with columns: `species`, optional metadata
#'   (`summary`, taxonomy ranks, `img_url`, `palette`, `image_file`), and all
#'   available trait columns found via TRY/local joins. `image_file` contains
#'   the normalized path to the PNG used for palette/preview (or `NA` if none).
#'
#' @details
#' \itemize{
#'   \item **Wikipedia:** summary via REST API; taxonomy parsed from the infobox.
#'   \item **Image:** first infobox image is used; when `remove_bg = TRUE`
#'         the function calls the remove.bg API. Set `Sys.setenv(REMOVE_BG_API_KEY = "…")`.
#'   \item **Palette:** simple k-means on non-transparent pixels of the (processed) PNG.
#'   \item **Traits (TRY):** wide table produced from `TraitName` and numeric `OrigValueStr`.
#'   \item **Traits (local):** returns all columns except the species column.
#'   \item **Dependencies:** `dplyr`, `purrr`, `tibble`, `fuzzyjoin`, `rvest`,
#'         `httr`, `stringr`, `jsonlite`, `magick`, `abind`.
#' }
#'
#' @importFrom rlang .data
#' @importFrom dplyr first
#'
#' @examples
#' \dontrun{
#' # Using TRY table
#' get_trait_data("Acacia karroo",
#'                use_try = TRUE,
#'                try_data = try_traits,
#'                trait_species_col = "SpeciesName")
#'
#' # Using a local trait table
#' get_trait_data("Acraea horta", local_trait_df = traits, local_species_col = "species")
#'
#' # Metadata only, with background removal and saving to a folder
#' Sys.setenv(REMOVE_BG_API_KEY = "<your-removebg-key>")
#' get_trait_data("Acacia karroo", use_try = FALSE, remove_bg = TRUE, save_folder = "out/")
#' }
#' @export
get_trait_data = function(
    species,
    remove_bg = FALSE,
    do_palette = TRUE,
    do_taxonomy = TRUE,
    do_summary = TRUE,
    do_image = TRUE,
    bg_thresh = 80,        # kept for backward-compat signature; unused when remove_bg=TRUE
    green_delta = 20,      # kept for backward-compat signature; unused when remove_bg=TRUE
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
  on.exit(closeAllConnections(), add = TRUE)

  # ---- internal helper: null-coalescing ----
  `%||%` <- get("%||%", envir = parent.frame(), inherits = TRUE)
  if (is.null(`%||%`)) `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- internal helper: remove.bg client (in-scope) ----
  remove_bg_and_save <- function(image_url, outfile,
                                 api_key = Sys.getenv("REMOVE_BG_API_KEY")) {
    if (!nzchar(api_key)) stop("Set REMOVE_BG_API_KEY environment variable")
    # Allow directory as outfile; build filename from image_url
    if (isTRUE(file.info(outfile)$isdir)) {
      base <- tools::file_path_sans_ext(basename(image_url))
      outfile <- file.path(outfile, paste0(base, "_nobg.png"))
    }
    if (!grepl("\\.png$", outfile, ignore.case = TRUE))
      outfile <- paste0(outfile, ".png")

    res <- try(httr::POST(
      url   = "https://api.remove.bg/v1.0/removebg",
      httr::add_headers(`X-Api-Key` = api_key),
      body  = list(image_url = image_url, size = "auto"),
      encode = "multipart"
    ), silent = TRUE)

    if (inherits(res, "try-error") || is.null(res)) {
      message("Request failed before reaching server.")
      return(NA_character_)
    }
    if (httr::http_error(res) || httr::status_code(res) != 200L) {
      txt <- tryCatch(httr::content(res, "text", encoding = "UTF-8"),
                      error = function(e) "<no body>")
      message("API error ", httr::status_code(res), ": ", txt)
      return(NA_character_)
    }

    imgdata <- httr::content(res, "raw")
    if (length(imgdata) < 1000L) {
      txt <- tryCatch(httr::content(res, "text", encoding = "UTF-8"),
                      error = function(e) "<no body>")
      message("Unexpected small payload: ", txt)
      return(NA_character_)
    }

    dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
    con <- file(outfile, "wb"); on.exit(close(con), add = TRUE)
    writeBin(imgdata, con)
    normalizePath(outfile, mustWork = FALSE)
  }

  # ---- Wikipedia summary ----
  clean_col = function(x) iconv(as.character(x), from = "", to = "UTF-8", sub = "byte")
  get_wikipedia_summary = function(species) {
    page_title = gsub(" ", "_", species)
    url = paste0("https://en.wikipedia.org/api/rest_v1/page/summary/", URLencode(page_title, reserved = TRUE))
    resp = tryCatch(httr::GET(url), error = function(e) NULL)
    if (is.null(resp) || resp$status_code != 200) return(NA)
    cont = tryCatch(jsonlite::fromJSON(rawToChar(resp$content)), error = function(e) NA)
    if (is.list(cont) && !is.null(cont$extract)) cont$extract else NA
  }

  # ---- Wikipedia taxonomy via infobox table ----
  get_wikipedia_taxonomy = function(species) {
    url = paste0("https://en.wikipedia.org/wiki/", stringr::str_replace_all(species, " ", "_"))
    page = tryCatch(rvest::read_html(url), error = function(e) NULL)
    if (is.null(page)) return(setNames(rep(NA_character_, 5), c("Kingdom", "Phylum", "Class", "Order", "Family")))
    infobox = page %>% rvest::html_node("table.infobox")
    if (inherits(infobox, "xml_missing") || is.null(infobox))
      return(setNames(rep(NA_character_, 5), c("Kingdom", "Phylum", "Class", "Order", "Family")))
    tbl = infobox %>%
      rvest::html_table(fill = TRUE) %>%
      tibble::as_tibble(.name_repair = "minimal") %>%
      purrr::set_names(c("trait", "value")) %>%
      dplyr::filter(trait != "" & value != "") %>%
      dplyr::mutate(trait = stringr::str_remove(trait, ":$"), value = stringr::str_squish(value))
    ranks = c("Kingdom", "Phylum", "Class", "Order", "Family")
    vals = setNames(rep(NA_character_, length(ranks)), ranks)
    for (rk in ranks) {
      v = tbl$value[tolower(tbl$trait) == tolower(rk)]
      if (length(v) > 0) vals[rk] = v[1]
    }
    vals
  }

  # ---- Wikipedia image url ----
  get_wikipedia_image = function(species) {
    url = paste0("https://en.wikipedia.org/wiki/", gsub(" ", "_", species))
    html = tryCatch(rvest::read_html(url), error = function(e) NULL)
    if (is.null(html)) return(NA)
    img_url = html %>% rvest::html_node(".infobox img") %>% rvest::html_attr("src")
    if (is.na(img_url) || is.null(img_url)) return(NA)
    if (!startsWith(img_url, "http")) img_url = paste0("https:", img_url)
    img_url
  }

  # ---- Palette extraction ----
  get_image_palette_local = function(img, n = 5) {
    img = magick::image_scale(img, "100")
    img_data = magick::image_data(img, channels = "rgba")
    alpha = as.integer(img_data[4, , ])
    keep = which(alpha > 0)
    if (length(keep) < n) return(NA)
    r = as.integer(img_data[1, , ])[keep]
    g = as.integer(img_data[2, , ])[keep]
    b = as.integer(img_data[3, , ])[keep]
    df = data.frame(r = r, g = g, b = b)
    if (nrow(unique(df)) < n) return(NA)
    km = kmeans(df, centers = n)
    rgb = round(km$centers)
    rgb_hex = apply(rgb, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    paste(rgb_hex, collapse = ", ")
  }

  # ---- TRY trait table (wide) ----
  try_traits = NULL
  if (use_try && !is.null(try_data)) {
    if (missing(trait_species_col) || is.null(trait_species_col) || trait_species_col == "" || !(trait_species_col %in% names(try_data))) {
      stop("trait_species_col must be a valid column name in try_data.")
    }
    trytab = try_data %>%
      dplyr::filter(!is.na(.data[[trait_species_col]]) & .data[[trait_species_col]] != "")
    trait_cols = intersect(c(trait_species_col, "TraitName", "OrigValueStr"), names(trytab))
    trytab = as.data.frame(trytab)[, trait_cols, drop = FALSE]
    suppressWarnings({
      trait_long = fuzzyjoin::stringdist_left_join(
        tibble::tibble(species = species),
        trytab,
        by = c("species" = trait_species_col),
        max_dist = max_dist,
        method = "jw"
      )
    })
    if (nrow(trait_long) > 0) {
      trait_wide = trait_long %>%
        dplyr::filter(!is.na(OrigValueStr)) %>%
        dplyr::mutate(
          OrigValueStr = iconv(OrigValueStr, from = "", to = "UTF-8", sub = NA),
          OrigValueStr = suppressWarnings(as.numeric(OrigValueStr))
        ) %>%
        dplyr::filter(!is.na(TraitName), TraitName != "", !is.na(OrigValueStr)) %>%
        dplyr::group_by(TraitName) %>%
        dplyr::summarise(value = dplyr::first(OrigValueStr), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = TraitName, values_from = value)
      if (nrow(trait_wide) > 0) try_traits = trait_wide
    }
  }

  # ---- Local trait table (any columns, wide) ----
  local_traits = NULL
  if (!is.null(local_trait_df) && is.data.frame(local_trait_df)) {
    if (!(local_species_col %in% names(local_trait_df))) {
      stop("local_species_col must be a valid column name in local_trait_df.")
    }
    localtab = local_trait_df %>%
      dplyr::filter(!is.na(.data[[local_species_col]]) & .data[[local_species_col]] != "")
    trait_cols = setdiff(names(localtab), local_species_col)
    if (length(trait_cols) == 0) stop("No trait columns found in local_trait_df.")
    suppressWarnings({
      trait_match = fuzzyjoin::stringdist_left_join(
        tibble::tibble(species = species),
        localtab,
        by = setNames(local_species_col, "species"),
        max_dist = max_dist,
        ignore_case = TRUE,
        distance_col = "distance"
      )
    })
    if (nrow(trait_match) > 0) {
      trait_match = trait_match %>% dplyr::select(dplyr::all_of(trait_cols))
      local_traits = trait_match[1, , drop = FALSE]
    }
  }

  # ---- Metadata scrape ----
  summary  = if (do_summary)  get_wikipedia_summary(species) else NA
  taxonomy = if (do_taxonomy) get_wikipedia_taxonomy(species) else setNames(rep(NA, 5), c("Kingdom","Phylum","Class","Order","Family"))
  img_url  = if (do_image)    get_wikipedia_image(species) else NA

  # ---- Image handling with remove.bg ----
  img_rgba <- NULL
  outfile  <- NULL

  if (do_image && !is.na(img_url)) {
    outdir <- save_folder %||% if (isTRUE(remove_bg)) tempdir() else NULL
    if (!is.null(outdir) && !dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    base_name <- gsub(" ", "_", species)

    if (isTRUE(remove_bg)) {
      target_path <- if (!is.null(outdir)) file.path(outdir, paste0(base_name, "_nobg.png")) else
        file.path(tempdir(), paste0(base_name, "_nobg.png"))
      res_path <- remove_bg_and_save(image_url = img_url, outfile = target_path)
      if (is.na(res_path)) {
        warning("Background removal failed; falling back to original image.")
        img_rgba <- tryCatch(magick::image_read(img_url), error = function(e) NULL)
        if (!is.null(save_folder) && inherits(img_rgba, "magick-image")) {
          outfile <- file.path(save_folder, paste0(base_name, "_orig.png"))
          magick::image_write(img_rgba, outfile)
        }
      } else {
        outfile <- res_path
        img_rgba <- tryCatch(magick::image_read(res_path), error = function(e) NULL)
      }
    } else {
      img_rgba <- tryCatch(magick::image_read(img_url), error = function(e) NULL)
      if (!is.null(save_folder) && inherits(img_rgba, "magick-image")) {
        outfile <- file.path(save_folder, paste0(base_name, "_orig.png"))
        magick::image_write(img_rgba, outfile)
      }
    }
  }

  if (isTRUE(preview) && inherits(img_rgba, "magick-image")) print(img_rgba)

  palette = if (do_palette && inherits(img_rgba, "magick-image")) get_image_palette_local(img_rgba, n = n_palette) else NA

  # ---- Output row ----
  output = tibble::tibble(
    species = species,
    summary = summary,
    Kingdom = taxonomy["Kingdom"],
    Phylum  = taxonomy["Phylum"],
    Class   = taxonomy["Class"],
    Order   = taxonomy["Order"],
    Family  = taxonomy["Family"],
    img_url = img_url,
    palette = palette,
    image_file = outfile %||% NA_character_
  )
  if (!is.null(try_traits))  output = dplyr::bind_cols(output, try_traits)
  if (!is.null(local_traits)) output = dplyr::bind_cols(output, local_traits)
  output
}
