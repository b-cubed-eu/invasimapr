#' Compute trait-based competition coefficients from a distance matrix
#'
#' @description
#' Given a species-by-species **trait distance** matrix (e.g., Gower), this function
#' extracts the **invader x resident** block of distances \eqn{d_{ij}} and transforms
#' it to a **competition coefficient** matrix \eqn{a_{ij}} using a Gaussian kernel:
#' \deqn{a_{ij} = \exp\!\left\{ - \frac{d_{ij}^2}{2\,\sigma_t^2} \right\}.}
#'
#' The kernel bandwidth \eqn{\sigma_t} controls how quickly competition decays with
#' trait difference. By default it is estimated from the **resident-resident** distances,
#' so the decay reflects the realised dispersion of the resident pool.
#'
#' @param g_all matrix. Symmetric species x species **trait distance** matrix
#'   with row/column names as species IDs (e.g., from Gower).
#' @param residents character. Vector of resident species IDs (must be in \code{rownames(g_all)}).
#' @param invaders NULL or character. Vector of invader species IDs. If \code{NULL},
#'   defaults to all species in \code{g_all} that are not in \code{residents}.
#' @param sigma_t NULL or numeric. Kernel bandwidth. If \code{NULL}, estimated from the
#'   upper triangle of the resident-resident distances using \code{sigma_method}.
#' @param sigma_method character. How to estimate \eqn{\sigma_t} when missing; one of
#'   \code{c("sd","median","iqr")}. Default \code{"sd"}.
#' @param eps numeric. Tiny positive value used if the estimated \eqn{\sigma_t} is zero;
#'   default \code{1e-8}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{d_ij} \code{(n_inv x n_res)}: invader-resident distance block (from \code{g_all}).
#'   \item \code{a_ij} \code{(n_inv x n_res)}: competition coefficients (Gaussian kernel of distances \eqn{d_{ij}}).
#'   \item \code{sigma_t}: bandwidth used.
#'   \item \code{meta}: list with estimation method and counts of invaders/residents.
#' }
#'
#' @details
#' - Use a **distance** matrix for \code{g_all} (not similarity). For Gower, distances are typically in [0, 1].
#' - Estimating \eqn{\sigma_t} from residents makes the kernel scale **data-driven** and comparable across analyses.
#' - If you prefer a fixed ecological scale, pass \code{sigma_t} explicitly (e.g., a trait threshold of interest).
#'
#' @examples
#' # --- Simulate a toy distance matrix (no external packages) ---
#' set.seed(42)
#' spp <- paste0("sp_", 1:8)
#' traits <- data.frame(
#'   t1 = rnorm(length(spp)),
#'   t2 = runif(length(spp))
#' )
#'
#' # Euclidean distances on simulated traits, scaled to [0, 1]
#' g_all <- as.matrix(stats::dist(traits, method = "euclidean"))
#' g_all <- g_all / max(g_all)                 # scale to [0,1]
#' rownames(g_all) <- colnames(g_all) <- spp
#'
#' # Define residents and invaders
#' residents <- spp[1:5]
#' invaders  <- spp[6:8]
#'
#' # Compute competition kernel (sigma_t estimated from resident-resident distances)
#' comp <- compute_competition_kernel(
#'   g_all     = g_all,
#'   residents = residents,
#'   invaders  = invaders,
#'   sigma_method = "sd"
#' )
#'
#' # Inspect results
#' str(comp$a_ij)           # (invaders x residents) competition coefficients
#' comp$sigma_t             # bandwidth used
#'
#' # Same call without providing invaders (uses complement of residents)
#' comp2 <- compute_competition_kernel(g_all, residents = residents)
#'
#' @importFrom stats sd IQR median
#' @export
compute_competition_kernel <- function(
    g_all,
    residents,
    invaders = NULL,
    sigma_t = NULL,
    sigma_method = c("sd", "median", "iqr"),
    eps = 1e-8) {
  # ---- input checks ----------------------------------------------------------
  if (!is.matrix(g_all)) stop("`g_all` must be a matrix of trait distances.")
  if (is.null(rownames(g_all)) || is.null(colnames(g_all))) {
    stop("`g_all` must have row and column names (species IDs).")
  }

  sigma_method <- match.arg(sigma_method)

  # Ensure residents exist in the matrix
  if (!all(residents %in% rownames(g_all))) {
    stop("Some `residents` are not found in rownames(g_all).")
  }
  if (!all(residents %in% colnames(g_all))) {
    stop("Some `residents` are not found in colnames(g_all).")
  }

  # If invaders not provided, take the complement
  if (is.null(invaders)) {
    invaders <- setdiff(rownames(g_all), residents)
  }
  if (!length(invaders)) stop("No invaders found (empty `invaders`).")
  if (!all(invaders %in% rownames(g_all)) || !all(invaders %in% colnames(g_all))) {
    stop("Some `invaders` are not present in `g_all` rows/cols.")
  }

  # ---- extract blocks --------------------------------------------------------
  # Distances among residents (upper triangle) to estimate sigma_t if needed
  D_rr <- g_all[residents, residents, drop = FALSE]
  ut <- D_rr[upper.tri(D_rr, diag = FALSE)]
  ut <- ut[is.finite(ut)]

  # Distances between each invader (rows) and each resident (cols)
  d_ij <- g_all[invaders, residents, drop = FALSE]

  # ---- bandwidth selection ---------------------------------------------------
  if (is.null(sigma_t)) {
    if (!length(ut)) {
      warning("No resident-resident distances available to estimate sigma_t; using eps.")
      sigma_t <- eps
    } else {
      sigma_t <- switch(sigma_method,
        sd     = stats::sd(ut),
        median = stats::median(ut),
        iqr    = stats::IQR(ut)
      )
      if (!is.finite(sigma_t) || sigma_t <= 0) sigma_t <- eps
    }
  }

  # ---- Gaussian kernel: a_ij = exp(-d_ij^2 / (2*sigma_t^2)) -----------------
  a_ij <- exp(-(d_ij^2) / (2 * sigma_t^2))

  # ---- return ----------------------------------------------------------------
  list(
    d_ij = d_ij,
    a_ij = a_ij,
    sigma_t = sigma_t,
    meta = list(
      sigma_method = if (missing(sigma_t) || is.null(sigma_t)) sigma_method else "user",
      n_invaders   = nrow(d_ij),
      n_residents  = ncol(d_ij)
    )
  )
}
