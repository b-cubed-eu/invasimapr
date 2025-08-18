#' Compute trait-based interaction strengths and resident abundance matrix
#'
#' @description
#' Computes a trait-based pairwise interaction matrix for all species (residents + invaders)
#' using Gower distance (default) or another dissimilarity metric, with options for:
#' - **Kernelisation** (`distance`, `similarity`, `gaussian`)
#' - **Scaling/standardisation**
#' - **Sparsification** (e.g., k-nearest neighbour graph)
#'
#' Also returns \code{Nstar}, a site × resident abundance matrix from predicted values.
#'
#' @details
#' The trait-based interaction strength matrix (\code{g_all}) is used in community assembly
#' and invasion models to quantify potential biotic effects between species. We recommend
#' \strong{Gower distance} for mixed traits, as it is scale-free and bounded in [0,1].
#'
#' The resident abundance matrix (\code{Nstar}) summarises equilibrium or expected abundances
#' for resident species at each site, based on predictions from the fitted abundance model.
#'
#' **Kernelisation options:**
#' \itemize{
#'   \item `"distance"`: return dissimilarities directly (default for transparent interpretation).
#'   \item `"similarity"`: convert to similarity in [0,1] by \code{1 - scaled_distance}.
#'   \item `"gaussian"`: Gaussian kernel \eqn{K = exp(-D^2 / (2\sigma^2))}, with automatic
#'         \eqn{\sigma} as the median non-zero distance if not supplied.
#' }
#'
#' @param traits data.frame. Species traits, first column = species name, remaining = traits.
#' @param predDF data.frame. Long table with columns `site_id`, `species`, and predicted abundance.
#' @param method character. Distance metric for \code{\link[cluster]{daisy}} (default `"gower"`).
#' @param kernel character. One of `"distance"`, `"similarity"`, `"gaussian"`.
#' @param sigma numeric. Gaussian kernel bandwidth; if NULL, set to median non-zero distance.
#' @param standardise logical. If TRUE, scale distances to [0,1] before kernelisation.
#' @param sparsify_k integer. If non-NULL, retain only k-nearest neighbours per species (symmetrised).
#'
#' @return A list with components:
#' \describe{
#'   \item{g_all}{Processed interaction matrix (distance / similarity / Gaussian kernel).}
#'   \item{raw_distance}{Unscaled pairwise distance matrix (e.g., Gower).}
#'   \item{Nstar}{Site × resident abundance matrix (rows = residents, cols = sites).}
#'   \item{sigma}{Bandwidth used when \code{kernel = "gaussian"}; otherwise \code{NULL}.}
#' }
#'
#' @examples
#' # --- Minimal, reproducible examples ---------------------------------------
#' set.seed(1)
#'
#' # Toy trait table: first column is species name (3 residents + 2 invaders)
#' traits <- data.frame(
#'   species = c("sp1","sp2","sp3","inv1","inv2"),
#'   body_size = c(10, 12, 8, 11, 9),                # numeric
#'   diet      = factor(c("herb","herb","omn","omn","herb")),  # factor
#'   activity  = c(0.2, 0.8, 0.5, 0.6, 0.4)          # numeric
#' )
#'
#' # Long table of predicted resident abundances across 4 sites
#' sites <- paste0("s", 1:4)
#' predDF <- expand.grid(site_id = sites, species = traits$species,
#'                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#' predDF$pred <- rexp(nrow(predDF), rate = 1)
#'
#' # --- 1) Plain distance (default): Gower dissimilarities --------------------
#' res_dist <- compute_interaction_strength(
#'   traits  = traits,
#'   predDF  = predDF,
#'   method  = "gower",
#'   kernel  = "distance"
#' )
#' # Inspect
#' round(res_dist$raw_distance, 3)
#' round(res_dist$g_all, 3)       # same as raw distances for kernel = "distance"
#' res_dist$Nstar[, 1:2]          # first two sites (residents x sites)
#'
#' # --- 2) Similarity in [0,1] (1 - scaled distance) --------------------------
#' res_sim <- compute_interaction_strength(
#'   traits      = traits,
#'   predDF      = predDF,
#'   kernel      = "similarity",
#'   standardise = TRUE            # scale distances to [0,1] before 1 - D
#' )
#' round(res_sim$g_all, 3)        # higher = more similar (stronger interaction)
#'
#' # --- 3) Gaussian kernel on distances ---------------------------------------
#' #      sigma chosen automatically as median non-zero distance
#' res_gauss <- compute_interaction_strength(
#'   traits = traits,
#'   predDF = predDF,
#'   kernel = "gaussian"          # K = exp(-D^2 / (2*sigma^2))
#' )
#' round(res_gauss$g_all, 3)
#' res_gauss$sigma                # bandwidth actually used
#'
#' # --- 4) Sparsify to a k-NN graph (here k = 2), symmetrised -----------------
#' res_sparse <- compute_interaction_strength(
#'   traits     = traits,
#'   predDF     = predDF,
#'   kernel     = "similarity",
#'   standardise = TRUE,
#'   sparsify_k = 2
#' )
#' round(res_sparse$g_all, 3)     # only each species' 2 strongest links retained
#'
#' # Notes:
#' # - Residents are inferred as species whose names do NOT start with "inv".
#' # - For mixed traits, Gower ("gower") is recommended; keep `standardise = TRUE`
#' #   when converting to similarity or when kernelising.
#'
#' @importFrom cluster daisy
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#'
#' @export
compute_interaction_strength <- function(
    traits,
    predDF,
    method = "gower",
    kernel = c("distance", "similarity", "gaussian"),
    sigma = NULL,
    standardise = FALSE,
    sparsify_k = NULL) {
  stopifnot(is.data.frame(traits), ncol(traits) >= 2)
  stopifnot(all(c("site_id", "species", "pred") %in% names(predDF)))

  kernel <- match.arg(kernel)

  # --- Resident species abundances -------------------------------------------
  residents <- setdiff(traits[[1]], grep("^inv", traits[[1]], value = TRUE))
  Nstar <- predDF %>%
    dplyr::filter(species %in% residents) %>%
    dplyr::select(species, site_id, pred) %>%
    tidyr::pivot_wider(names_from = site_id, values_from = pred) %>%
    tibble::column_to_rownames("species") %>%
    as.matrix()

  # --- Pairwise Gower (or other) distances -----------------------------------
  raw_distance <- cluster::daisy(traits[, -1], metric = method) %>% as.matrix()
  rownames(raw_distance) <- colnames(raw_distance) <- traits[[1]]

  # --- Optional scaling ------------------------------------------------------
  scale01 <- function(M) {
    rng <- range(M, na.rm = TRUE)
    if (diff(rng) == 0) {
      return(M * 0)
    }
    (M - rng[1]) / diff(rng)
  }
  D <- if (standardise || (kernel == "similarity" && method != "gower")) scale01(raw_distance) else raw_distance

  # --- Kernelisation ---------------------------------------------------------
  if (kernel == "similarity") {
    g_all <- 1 - D
  } else if (kernel == "gaussian") {
    if (is.null(sigma)) {
      vals <- D[upper.tri(D)]
      sigma <- median(vals[vals > 0], na.rm = TRUE)
    }
    g_all <- exp(-(D^2) / (2 * sigma^2))
    diag(g_all) <- 0
  } else {
    g_all <- D
  }

  # --- Optional sparsification -----------------------------------------------
  if (!is.null(sparsify_k)) {
    sparsify_knn <- function(K, k) {
      n <- nrow(K)
      A <- matrix(0, n, n, dimnames = dimnames(K))
      for (i in 1:n) {
        idx <- order(K[i, ], decreasing = TRUE)[seq_len(min(k, n))]
        A[i, idx] <- K[i, idx]
      }
      pmax(A, t(A))
    }
    g_all <- sparsify_knn(g_all, sparsify_k)
  }

  list(
    g_all        = g_all,
    raw_distance = raw_distance,
    Nstar        = Nstar,
    sigma        = if (kernel == "gaussian") sigma else NULL
  )
}
