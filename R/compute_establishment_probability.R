#' Convert invasion fitness to probabilistic establishment (probit/logit/hard)
#'
#' @title Probabilistic establishment from invasion fitness
#'
#' @description
#' Invasion fitness \eqn{\lambda_{is}} integrates **trait space geometry** (distances,
#' overlaps, **convex hull**, **cloud centroid**) with **abiotic suitability**
#' (alignment of invader traits to local environment), **niche crowding** (overlap
#' with resident trait space weighted by composition), and **resident competition**
#' (site saturation). `compute_establishment_probability()` maps \eqn{\lambda_{is}}
#' to **probabilities of establishment** using a unified interface:
#'
#' - **Probit**: \eqn{P = \Phi(\lambda/\sigma)} with \eqn{\sigma} as a scalar,
#'   the residual SD from a fitted auxiliary GLMM, or a cell-wise **predictive SD**.
#' - **Logistic**: \eqn{P = \mathrm{logit}^{-1}(\lambda / \tau)} with scale \eqn{\tau}.
#' - **Hard rule**: \eqn{P = \mathbb{I}\{\lambda>0\}} for a binary map.
#'
#' If `lambda_is` is not supplied, the function will build it from standardized
#' predictors via \eqn{\lambda_{is} = \gamma \, r^{(z)}_{is} - \alpha \, C^{(z)}_{is} - \beta \, S^{(z)}_{is} + \kappa}.
#'
#' @param lambda_is Optional matrix \eqn{S\times I} of invasion fitness (rows = sites,
#'   cols = invaders). If `NULL`, it is computed from components (see below).
#' @param r_is_z,C_is_z,S_is_z Optional matrices \eqn{S\times I} of standardized
#'   abiotic suitability, niche crowding, and saturation. Used only if `lambda_is=NULL`.
#' @param gamma Optional vector (length \eqn{I}) or matrix \eqn{S\times I} for the
#'   abiotic slope \eqn{\gamma}; recycled to \eqn{S\times I} as needed.
#' @param alpha Optional vector (length \eqn{I}) or matrix \eqn{S\times I} of crowding
#'   penalties \eqn{\alpha \ge 0} by convention); recycled if needed.
#' @param beta  Optional vector (length \eqn{I}) of saturation penalties \eqn{\beta}.
#'   If you allow facilitation, pass signed values here.
#' @param kappa Optional scalar calibration offset added to \eqn{\lambda}; default `0`.
#' @param method Character, one of `c("probit","logit","hard")`.
#' @param sigma Numeric **scalar** SD for probit (ignored for other methods). If not
#'   given and `fit` is supplied, `sigma` defaults to `sigma(fit)`.
#' @param tau Numeric **scalar** scale for logistic (denominator inside logit); default `1`.
#' @param fit Optional fitted model (e.g., `glmmTMB`), used to obtain `sigma(fit)` when
#'   `method="probit"` and `sigma` is missing.
#' @param predictive Logical; for `method="probit"`, if `TRUE` the function uses a
#'   **predictive SD** matrix instead of a scalar. Provide it via `sigma_mat`, or set
#'   `use_vcov=TRUE` to compute it from `fit` (see below).
#' @param sigma_mat Optional matrix \eqn{S\times I} of SDs (e.g., predictive SD per cell).
#' @param use_vcov Logical; if `TRUE` and `method="probit"` with `predictive=TRUE`, the
#'   function will attempt to compute **mean-SE** from `vcov(fit)` and add the residual
#'   SD to yield a predictive SD. Requires `fit`, `Q_inv`, and the standardized matrices.
#' @param Q_inv Optional data frame with invader trait-plane scores (`tr1`,`tr2`),
#'   rownames = invader IDs. Required only if `use_vcov=TRUE`.
#' @param site_df Optional site table with columns `site,x,y` for mapping; rownames of
#'   the matrices must match `site_df$site`.
#' @param return_long Logical; if `TRUE`, include a tidy tibble in the output.
#' @param make_plots Logical; if `TRUE`, return ggplot objects (site mean map, invader
#'   ranking, and site×invader heatmap).
#' @param option_label Optional label attached to the tidy output and plot titles.
#'
#' @details
#' **Probit**: \eqn{P_{is}=\Phi(\lambda_{is}/\sigma)}. Use a **scalar** `sigma` (fast),
#' the residual SD from `fit`, or set `predictive=TRUE` to inject a **cell-wise** SD:
#' either supply `sigma_mat` directly or compute it with `use_vcov=TRUE`, which uses
#' the fixed-effect covariance from `vcov(fit)` plus the model residual SD.
#'
#' **Logistic**: \eqn{P_{is}=\mathrm{logit}^{-1}(\lambda_{is}/\tau)}. The scale `tau`
#' controls how sharply probabilities switch around \eqn{\lambda=0}.
#'
#' **Hard rule**: \eqn{P_{is} = \mathbb{I}\{\lambda_{is}>0\}} yields a binary map
#' useful for thresholds and counts (invasibility = \#invaders with \(P=1\) at a site).
#'
#' @return A list with:
#' - `p_is`: matrix \eqn{S\times I} of probabilities (or 0/1 for `hard`);
#' - `lambda_is`: the \eqn{S\times I} fitness matrix used;
#' - `sigma_used`: scalar or matrix actually used for probit (else `NULL`);
#' - `method`, `option_label`;
#' - `prob_long`: tidy tibble (if `return_long=TRUE`);
#' - `plots`: list of ggplots (if `make_plots=TRUE`): `site_mean`, `invader_mean`, `heatmap`.
#'
#' @examples
#' ## Minimal example (toy shapes)
#' set.seed(1)
#' S = 6; I = 4
#' sites   = paste0("s", 1:S)
#' inv     = paste0("i", 1:I)
#' r_is_z  = matrix(rnorm(S*I), S, I, dimnames=list(sites, inv))
#' C_is_z  = matrix(rnorm(S*I), S, I, dimnames=dimnames(r_is_z))
#' S_is_z  = matrix(rep(scale(rnorm(S)), each=I), S, I, dimnames=dimnames(r_is_z))
#' gamma   = setNames(runif(I, 0.5, 1.2), inv)
#' alpha   = setNames(runif(I, 0.2, 1.0), inv)
#' beta    = setNames(runif(I, 0.1, 0.6), inv)
#'
#' # Build lambda internally, then get logistic probabilities
#' out_logit = compute_establishment_probability(
#'   r_is_z=r_is_z, C_is_z=C_is_z, S_is_z=S_is_z,
#'   gamma=gamma, alpha=alpha, beta=beta,
#'   method="logit", tau=1, return_long=FALSE, make_plots=FALSE
#' )
#' str(out_logit$p_is)
#'
#' # Probit with a scalar sigma
#' out_probit = compute_establishment_probability(
#'   r_is_z=r_is_z, C_is_z=C_is_z, S_is_z=S_is_z,
#'   gamma=gamma, alpha=alpha, beta=beta,
#'   method="probit", sigma=1
#' )
#' # View site-mean probability map (requires ggplot2)
#' if (requireNamespace("ggplot2", quietly=TRUE)) print(out_probit$plots$site_mean)
#'
#' # Hard rule (λ>0)
#' out_hard = compute_establishment_probability(
#'   r_is_z=r_is_z, C_is_z=C_is_z, S_is_z=S_is_z,
#'   gamma=gamma, alpha=alpha, beta=beta,
#'   method="hard", return_long=TRUE, make_plots=FALSE
#' )
#' table(out_hard$prob_long$val)  # 0/1 counts
#'
#' @export
compute_establishment_probability = function(
    lambda_is = NULL,
    r_is_z = NULL, C_is_z = NULL, S_is_z = NULL,
    gamma = 1, alpha = NULL, beta = NULL, kappa = 0,
    method = c("probit","logit","hard"),
    sigma = NULL, tau = 1,
    fit = NULL,
    predictive = FALSE, sigma_mat = NULL, use_vcov = FALSE,
    Q_inv = NULL,
    site_df = NULL, return_long = TRUE, make_plots = TRUE,
    option_label = NULL
){
  method = match.arg(method)

  # ---- 0) Build or validate lambda -------------------------------------------
  if (is.null(lambda_is)) {
    # Must have components + coefficients
    for (nm in c("r_is_z","C_is_z","S_is_z")) {
      if (is.null(get(nm))) stop("When lambda_is is NULL, you must provide r_is_z, C_is_z, S_is_z.")
    }
    stopifnot(identical(dim(r_is_z), dim(C_is_z)), identical(dim(r_is_z), dim(S_is_z)))
    S = nrow(r_is_z); I = ncol(r_is_z)
    sites   = rownames(r_is_z); invaders = colnames(r_is_z)
    if (is.null(sites) || is.null(invaders)) stop("Matrices must have rownames (sites) and colnames (invaders).")

    # Recycle gamma/alpha to S×I; beta is per-invader (vector)
    if (length(gamma) == 1L) {
      GI = matrix(gamma, S, I, dimnames=list(sites, invaders))
    } else if (is.vector(gamma)) {
      stopifnot(!is.null(names(gamma)))
      GI = matrix(rep(gamma[invaders], each=S), S, I, dimnames=list(sites, invaders))
    } else {
      stopifnot(is.matrix(gamma), identical(dim(gamma), dim(r_is_z))); GI = gamma
    }
    if (is.null(alpha)) stop("alpha is required when lambda_is is NULL.")
    if (is.vector(alpha)) {
      stopifnot(!is.null(names(alpha)))
      AI = matrix(rep(alpha[invaders], each=S), S, I, dimnames=list(sites, invaders))
    } else {
      stopifnot(is.matrix(alpha), identical(dim(alpha), dim(r_is_z))); AI = alpha
    }
    if (is.null(beta)) stop("beta is required when lambda_is is NULL.")
    stopifnot(is.vector(beta), !is.null(names(beta)))
    BI = beta[invaders]

    lambda_is = GI * r_is_z - AI * C_is_z - sweep(S_is_z, 2, BI, `*`) + kappa
  } else {
    stopifnot(is.matrix(lambda_is), !is.null(rownames(lambda_is)), !is.null(colnames(lambda_is)))
    S = nrow(lambda_is); I = ncol(lambda_is)
    sites   = rownames(lambda_is); invaders = colnames(lambda_is)
  }

  # ---- 1) Choose SD/scale for transformation ---------------------------------
  sigma_used = NULL
  if (method == "probit") {
    if (isTRUE(predictive)) {
      # a) user-supplied predictive SD matrix
      if (!is.null(sigma_mat)) {
        stopifnot(is.matrix(sigma_mat), identical(dim(sigma_mat), dim(lambda_is)))
        sigma_used = sigma_mat
      } else if (isTRUE(use_vcov)) {
        # b) compute mean-SE from vcov(fit) and add residual SD (predictive)
        if (is.null(fit)) stop("use_vcov=TRUE requires a fitted model in `fit`.")
        if (is.null(r_is_z) || is.null(C_is_z) || is.null(S_is_z))
          stop("vcov-based predictive SD needs r_is_z, C_is_z, S_is_z.")
        if (is.null(Q_inv) || !all(c("tr1","tr2") %in% colnames(Q_inv)))
          stop("vcov-based predictive SD needs Q_inv with columns tr1 and tr2 (rownames=invader IDs).")

        # Align trait-plane rows to invader columns
        Q_inv = Q_inv[invaders, c("tr1","tr2"), drop=FALSE]

        # Fixed-effect coefficients and covariance (glmmTMB/lme4 style)
        cf  = tryCatch(stats::coef(fit)$cond, error=function(e) NULL)
        if (is.null(cf)) cf = tryCatch(glmmTMB::fixef(fit)$cond, error=function(e) NULL)
        if (is.null(cf)) stop("Could not extract fixed effects from `fit`.")

        Vb  = tryCatch(as.matrix(stats::vcov(fit)$cond), error=function(e) NULL)
        if (is.null(Vb)) stop("Could not extract vcov from `fit`.")
        sig = tryCatch(glmmTMB::sigma(fit), error=function(e) NA_real_)
        if (!is.finite(sig)) sig = 1

        # Helper to pull coefficients robustly (handles a:b vs b:a)
        get_cf = function(cfvec, a, b=NULL, default=0){
          if (is.null(b)) {
            if (a %in% names(cfvec)) unname(cfvec[a]) else default
          } else {
            nm = c(paste0(a,":",b), paste0(b,":",a))
            nm = nm[nm %in% names(cfvec)]
            if (length(nm)) unname(cfvec[nm[1]]) else default
          }
        }

        # Parameter name universe
        pn = names(cf); Vb = Vb[pn, pn, drop=FALSE]

        # Build x-row for each (s,i): matches the auxiliary linear predictor used for λ
        build_x_row = function(param_names, r_sz, C_sz, S_sz, tr1_i, tr2_i){
          x = setNames(numeric(length(param_names)), param_names)

          add = function(term, val) if (term %in% names(x)) x[term] <= x[term] + val
          # + r terms (positive)
          add("r_z", r_sz)
          add("r_z:tr1", r_sz*tr1_i); add("tr1:r_z", r_sz*tr1_i)
          add("r_z:tr2", r_sz*tr2_i); add("tr2:r_z", r_sz*tr2_i)
          # - C terms (penalty)
          add("C_z", -C_sz)
          add("C_z:tr1", -C_sz*tr1_i); add("tr1:C_z", -C_sz*tr1_i)
          add("C_z:tr2", -C_sz*tr2_i); add("tr2:C_z", -C_sz*tr2_i)
          # - S terms (penalty or signed if you used signed beta when building λ)
          add("S_z", -S_sz)
          add("S_z:tr1", -S_sz*tr1_i); add("tr1:S_z", -S_sz*tr1_i)
          add("S_z:tr2", -S_sz*tr2_i); add("tr2:S_z", -S_sz*tr2_i)
          x
        }

        sigma_used = matrix(NA_real_, S, I, dimnames=dimnames(lambda_is))
        for (j in seq_len(I)) {
          tr1j = Q_inv$tr1[j]; tr2j = Q_inv$tr2[j]
          for (i in seq_len(S)) {
            x  = build_x_row(pn, r_is_z[i,j], C_is_z[i,j], S_is_z[i,j], tr1j, tr2j)
            vv = as.numeric(t(x) %*% Vb %*% x)     # var of mean predictor
            se_mean = sqrt(max(vv, 0))
            sigma_used[i,j] = sqrt(se_mean^2 + sig^2)  # predictive SD
          }
        }
      } else {
        stop("predictive=TRUE requires either `sigma_mat` or `use_vcov=TRUE` + `fit` + predictors + Q_inv.")
      }
    } else {
      # Residual (scalar) probit SD
      if (is.null(sigma) && !is.null(fit)) {
        sigma = tryCatch(glmmTMB::sigma(fit), error=function(e) NA_real_)
      }
      if (!is.finite(sigma) || sigma <= 0) sigma = 1
      sigma_used = sigma
    }
  } else if (method == "logit") {
    if (!is.finite(tau) || tau <= 0) tau = 1
  }

  # ---- 2) Transform λ -> probability -----------------------------------------
  p_is =
    if (method == "probit") {
      if (is.matrix(sigma_used)) {
        stats::pnorm(lambda_is / sigma_used)
      } else {
        stats::pnorm(lambda_is / as.numeric(sigma_used))
      }
    } else if (method == "logit") {
      stats::plogis(lambda_is / tau)
    } else { # hard
      (lambda_is > 0) * 1L
    }

  # ---- 3) Tidy long table -----------------------------------------------------
  prob_long = NULL
  if (isTRUE(return_long)) {
    prob_long = tibble::tibble(
      site    = rep(rownames(p_is), times = ncol(p_is)),
      invader = rep(colnames(p_is), each  = nrow(p_is)),
      val     = as.vector(p_is),
      lambda  = as.vector(lambda_is),
      option  = option_label %||% switch(method,
                                         probit="Probit P(establish)",
                                         logit ="Logistic P(establish)",
                                         hard  ="Establish (λ>0)")
    )
    if (!is.null(site_df) && all(c("site","x","y") %in% names(site_df))) {
      prob_long = dplyr::left_join(prob_long, site_df, by="site")
    }
  }

  if (isTRUE(return_long) && (is.null(prob_long) || NROW(prob_long) == 0L)) {
    prob_long = tibble::tibble(
      site    = rep(rownames(p_is), times = ncol(p_is)),
      invader = rep(colnames(p_is), each  = nrow(p_is)),
      val     = as.vector(p_is),
      lambda  = as.vector(lambda_is),
      option  = option_label %||% method
    )
    if (!is.null(site_df) && all(c("site","x","y") %in% names(site_df))) {
      prob_long = dplyr::left_join(prob_long, site_df, by="site")
    }
  }


  # ---- 4) Plots ---------------------------------------------------------------
  plots = NULL
  if (isTRUE(make_plots) && requireNamespace("ggplot2", quietly=TRUE)) {
    d = if (is.null(prob_long)) {
      tibble::tibble(
        site    = rep(rownames(p_is), times = ncol(p_is)),
        invader = rep(colnames(p_is), each  = nrow(p_is)),
        val     = as.vector(p_is)
      )
    } else prob_long

    # Site mean (invasibility map)
    if (!is.null(site_df) && all(c("x","y") %in% names(d))) {
      site_mean = d |>
        dplyr::group_by(site, x, y) |>
        dplyr::summarise(val = mean(val, na.rm=TRUE), .groups="drop") |>
        ggplot2::ggplot(ggplot2::aes(x, y, fill = val)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(midpoint = if (method=="hard") 0.5 else 0,
                                      name = if (method=="hard") "# establishing" else "Mean P") +
        ggplot2::labs(title = (option_label %||% "Establishment"), x="x", y="y") +
        ggplot2::theme_minimal()
    } else site_mean = NULL

    # Invader mean (invasiveness ranking)
    invader_mean = d |>
      dplyr::group_by(invader) |>
      dplyr::summarise(val = mean(val, na.rm=TRUE), .groups="drop") |>
      ggplot2::ggplot(ggplot2::aes(x = stats::reorder(invader, val), y = val)) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = if (method=="hard") ggplot2::waiver() else scales::percent_format()) +
      ggplot2::labs(x="Invader", y = if (method=="hard") "Share of sites (λ>0)" else "Mean P",
                    title = "Invasiveness (mean across sites)") +
      ggplot2::theme_minimal()

    # Heatmap (site × invader)
    if (method=='hard') {
      d$val_f = factor(d$val, levels = c(0,1), labels = c("0","1"))
    heatmap = d |>
      # dplyr::mutate(val_f = factor(val, levels = c(0,1), labels = c("0","1"))) |>
      ggplot2::ggplot(ggplot2::aes(invader, site, fill = val_f)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(values = c("0" = "darkgrey", "1" = "darkred"), name = "Establish (0/1)") +
      ggplot2::labs(title = "Establishment matrix", x="Invader", y="Site") +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, hjust=1))
  } else heatmap = d |>
      ggplot2::ggplot(ggplot2::aes(invader, site, fill = val)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(name = "P(establish)") +
      ggplot2::labs(title = "Establishment matrix", x="Invader", y="Site") +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle=90, vjust=0.5, hjust=1))

    plots = list(site_mean = site_mean, invader_mean = invader_mean, heatmap = heatmap)
  }

  list(
    p_is        = p_is,
    lambda_is   = lambda_is,
    sigma_used  = sigma_used,
    method      = method,
    option_label= option_label %||% method,
    prob_long   = prob_long,
    plots       = plots
  )
}
