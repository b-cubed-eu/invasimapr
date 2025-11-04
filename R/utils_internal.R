#' Internal utilities
#'
#' Small helper functions used inside the package. Not part of the public API.
#'
#' @name utils_internal
#' @keywords internal   # <-- hide from index
#' @aliases "%||%" new_invasimapr_fit print.invasimapr_fit .standardise_df .scale_like .row_z
NULL
new_invasimapr_fit = function(x = list()) {
  class(x) = unique(c("invasimapr_fit", class(x)))
  x
}

#' @title Compact print
#' @keywords internal   # <-- hide from index
#' @export
print.invasimapr_fit = function(x, ...) {
  cat("<invasimapr_fit>\n")
  stages = c("inputs","traits","crowding","residents","sensitivities","invaders","fitness","prob","summary")
  present = stages[stages %in% names(x)]
  if (length(present)) cat(" stages:", paste(present, collapse = " → "), "\n")
  if (!is.null(x$meta)) {
    cat(sprintf(" sites: %s | residents: %s | invaders: %s\n",
                x$meta$n_sites %||% NA_integer_,
                x$meta$n_residents %||% NA_integer_,
                x$meta$n_invaders %||% NA_integer_))
  }
  invisible(x)
}
`%||%` = function(a,b) if (!is.null(a)) a else b

# ------------------------------ small helpers ------------------------------

# column-wise z for data.frame (numeric columns only; factors/characters kept)
.standardise_df = function(df, ref_means = NULL, ref_sds = NULL) {
  stopifnot(is.data.frame(df))
  out = df
  num = vapply(df, is.numeric, logical(1))
  if (!any(num)) return(list(data = out, means = numeric(0), sds = numeric(0)))

  X = as.matrix(df[, num, drop = FALSE])
  if (is.null(ref_means) || is.null(ref_sds)) {
    m = colMeans(X, na.rm = TRUE)
    s = apply(X, 2, stats::sd, na.rm = TRUE)
  } else {
    m = ref_means[colnames(X)]
    s = ref_sds  [colnames(X)]
  }
  s[!is.finite(s) | s == 0] = 1
  Z = sweep(sweep(X, 2, m, "-"), 2, s, "/")

  out[, num] = Z
  list(data = out, means = m, sds = s)
}

# scale a df like a reference (means/sds taken from residents)
.scale_like = function(df, ref_means, ref_sds) {
  .standardise_df(df, ref_means, ref_sds)$data
}

# quick row-wise z (matrix)
.row_z = function(M, robust = FALSE) {
  std = standardise_by_site(as.matrix(M), robust = robust)
  std
}

# ------------------------------ sigma_mat_from_vcov helper ------------------------------

# Predictive SD for each (site s, invader i) cell from a glmmTMB fit
sigma_mat_from_vcov <- function(fit, r_is_z, C_is_z, S_is_z, Q_inv, add_resid = TRUE) {
  stopifnot(is.matrix(r_is_z), is.matrix(C_is_z), is.matrix(S_is_z))
  stopifnot(identical(dim(r_is_z), dim(C_is_z)), identical(dim(r_is_z), dim(S_is_z)))
  S <- nrow(r_is_z); I <- ncol(r_is_z)
  sites <- rownames(r_is_z); inv   <- colnames(r_is_z)

  # Align invader trait rows
  stopifnot(is.data.frame(Q_inv), all(c("tr1","tr2") %in% names(Q_inv)))
  Q_inv <- Q_inv[inv, c("tr1","tr2"), drop = FALSE]

  # Fixed effects & vcov
  cf <- tryCatch(glmmTMB::fixef(fit)$cond, error = function(e) NULL)
  if (is.null(cf)) stop("Could not extract fixed effects from `fit`.")
  Vb <- tryCatch(as.matrix(stats::vcov(fit)$cond), error = function(e) NULL)
  if (is.null(Vb)) stop("Could not extract vcov from `fit`.")

  pn <- names(cf); Vb <- Vb[pn, pn, drop = FALSE]
  sig_res <- if (isTRUE(add_resid)) {
    s <- tryCatch(glmmTMB::sigma(fit), error = function(e) NA_real_)
    if (!is.finite(s)) 1 else s
  } else 0

  # Build x-row matching (r_z + C_z + S_z) * (tr1 + tr2)
  build_x_row <- function(param_names, r_sz, C_sz, S_sz, tr1_i, tr2_i){
    x <- setNames(numeric(length(param_names)), param_names)
    add <- function(term, val) if (term %in% names(x)) x[term] <<- x[term] + val
    add("r_z", r_sz);              add("r_z:tr1", r_sz*tr1_i); add("tr1:r_z", r_sz*tr1_i)
    add("r_z:tr2", r_sz*tr2_i);    add("tr2:r_z", r_sz*tr2_i)
    add("C_z", -C_sz);             add("C_z:tr1", -C_sz*tr1_i); add("tr1:C_z", -C_sz*tr1_i)
    add("C_z:tr2", -C_sz*tr2_i);   add("tr2:C_z", -C_sz*tr2_i)
    add("S_z", -S_sz);             add("S_z:tr1", -S_sz*tr1_i); add("tr1:S_z", -S_sz*tr1_i)
    add("S_z:tr2", -S_sz*tr2_i);   add("tr2:S_z", -S_sz*tr2_i)
    x
  }

  Sig <- matrix(NA_real_, S, I, dimnames = list(sites, inv))
  for (j in seq_len(I)) {
    tr1j <- Q_inv$tr1[j]; tr2j <- Q_inv$tr2[j]
    for (i in seq_len(S)) {
      x <- build_x_row(pn, r_is_z[i,j], C_is_z[i,j], S_is_z[i,j], tr1j, tr2j)
      v_mean <- as.numeric(t(x) %*% Vb %*% x)  # variance of the mean predictor
      se_mean <- sqrt(max(v_mean, 0))
      Sig[i,j] <- sqrt(se_mean^2 + sig_res^2)  # predictive SD
    }
  }
  Sig
}

