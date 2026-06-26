#' @keywords internal
"_PACKAGE"

## Imports used without explicit namespace in several functions
#' @importFrom graphics axis contour filled.contour points title
#' @importFrom stats kmeans predict
#' @importFrom utils URLencode tail
NULL

## Tidy-eval / NSE symbols used across the package
utils::globalVariables(c(
  "species","site_id","pred","PCoA1","PCoA2","centrality",
  "Metric","Value","OrigValueStr","TraitName","trait","value",":=",
  # common columns
  "site", "x", "y", "invader", "abundance",

  # assemble_matrices
  "spp_rich",

  # trait space / hull
  "tr1", "tr2", "grp", "id", "rank_val", "in_hull", "d_md",

  # invasion fitness / establishment
  "val", "val_f", "C_mean",

  # modelling
  "make_lambda_long"
))

# # Keep startup hooks minimal: no data I/O, no model fits, no terra calls.
# .onLoad <- function(libname, pkgname) {
#   # Register methods or options only (if any). Avoid side effects.
# }
#
# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("invasimapr loaded. See ?invasimapr for help.")
# }
