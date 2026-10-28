#' Incremental analysis with strong and extended dominance
#' @param cost,effect Numeric vectors (one per strategy).
#' @param strategies Strategy names.
#' @return Tibble sorted by cost with `status` and `icer` (vs previous
#'   non-dominated strategy).
#' @export
icer_frontier <- function(cost, effect, strategies) {
  df <- tibble::tibble(strategy = strategies, cost = unname(cost), effect = unname(effect))
  df <- df[order(df$cost, -df$effect), ]
  df$status <- ""
  for (i in seq_len(nrow(df))) {
    o <- seq_len(nrow(df))[-i]
    if (any(df$cost[o] <= df$cost[i] & df$effect[o] >= df$effect[i] &
            (df$cost[o] < df$cost[i] | df$effect[o] > df$effect[i]))) df$status[i] <- "Dominated"
  }
  repeat {
    nd <- which(df$status == ""); if (length(nd) < 3) break
    ext <- which(diff(diff(df$cost[nd]) / diff(df$effect[nd])) < 0)
    if (!length(ext)) break
    df$status[nd[ext[1] + 1]] <- "Extendedly dominated"
  }
  nd <- which(df$status == ""); df$icer <- NA_real_
  if (length(nd) > 1) df$icer[nd[-1]] <- diff(df$cost[nd]) / diff(df$effect[nd])
  df$status[nd] <- "On frontier"
  df
}

#' Cost-effectiveness acceptability curves
#' @param cost,effect Matrices: simulations x strategies (named columns).
#' @param wtp_grid Willingness-to-pay values.
#' @return Tibble `wtp`, `arm`, `prob` (probability arm has the highest NMB),
#'   and `frontier` (TRUE for the arm with highest expected NMB).
#' @export
ceac <- function(cost, effect, wtp_grid) {
  arms <- colnames(cost)
  purrr::map_dfr(wtp_grid, function(w) {
    nmb <- w * effect - cost
    best <- apply(nmb, 1, which.max)
    tibble::tibble(wtp = w, arm = arms, prob = tabulate(best, length(arms)) / nrow(cost),
                   frontier = seq_along(arms) == which.max(colMeans(nmb)))
  })
}

#' Expected value of perfect information
#' @inheritParams ceac
#' @param population Population size for population EVPI (`NA` to skip).
#' @return Tibble `wtp`, `evpi` (per person) and `evpi_population`.
#' @export
evpi <- function(cost, effect, wtp_grid, population = NA) {
  purrr::map_dfr(wtp_grid, function(w) {
    nmb <- w * effect - cost
    v <- mean(apply(nmb, 1, max)) - max(colMeans(nmb))
    tibble::tibble(wtp = w, evpi = v, evpi_population = v * population)
  })
}

#' Expected value of partial perfect information (regression method)
#'
#' Non-parametric regression (Strong, Oakley & Brennan 2014). By default
#' (`method = "auto"`) the calculation is delegated to the voi package when it
#' is installed (GAM for up to 4 parameters, Gaussian process above that).
#' Otherwise a built-in GAM with voi's default basis (`te(..., bs = "cr")`) is
#' used: incremental NMB of each strategy vs the first is regressed on the
#' parameter(s); EVPPI = E(max fitted) - max E(fitted).
#' @inheritParams ceac
#' @param params Matrix or data frame of sampled parameter values
#'   (simulations x parameters).
#' @param wtp Willingness to pay.
#' @param groups Optional named list of parameter-name vectors for joint EVPPI.
#' @param method `"auto"` (voi if installed, else built-in), `"voi"` or `"gam"`.
#' @return Tibble `parameter`, `evppi` (per person), `group`, `method`,
#'   sorted descending.
#' @export
evppi <- function(cost, effect, params, wtp, groups = NULL, method = c("auto", "voi", "gam")) {
  method <- match.arg(method)
  if (method == "auto") method <- if (requireNamespace("voi", quietly = TRUE)) "voi" else "gam"
  params <- as.data.frame(params)
  varies <- vapply(params, function(v) length(unique(v)) >= 5, logical(1))
  sets <- c(stats::setNames(as.list(colnames(params)), colnames(params)), groups)
  usable <- vapply(sets, function(cols) all(cols %in% names(params)) && all(varies[cols]), logical(1))
  v <- rep(NA_real_, length(sets))
  if (method == "voi") {
    if (!requireNamespace("voi", quietly = TRUE)) stop("install the voi package or use method = 'gam'")
    out <- voi::evppi(outputs = list(c = as.matrix(cost), e = as.matrix(effect), k = wtp),
                      inputs = params, pars = unname(sets[usable]))
    v[usable] <- pmax(0, out$evppi)
  } else {
    if (!requireNamespace("mgcv", quietly = TRUE)) stop("install the mgcv package for EVPPI")
    nmb <- wtp * effect - cost
    inb <- nmb[, -1, drop = FALSE] - nmb[, 1]
    one <- function(cols) {
      X <- params[, cols, drop = FALSE]; names(X) <- paste0("x", seq_along(cols))
      karg <- if (length(cols) >= 4) ", k = 4" else ""
      f <- stats::as.formula(sprintf("y ~ te(%s, bs = 'cr'%s)", paste(names(X), collapse = ", "), karg))
      fit <- vapply(seq_len(ncol(inb)), function(j) {
        d <- X; d$y <- inb[, j]
        if (stats::sd(d$y) == 0) return(rep(mean(d$y), nrow(d)))
        stats::fitted(mgcv::gam(f, data = d))
      }, numeric(nrow(X)))
      g <- cbind(0, fit)
      max(0, mean(apply(g, 1, max)) - max(colMeans(g)))
    }
    v[usable] <- vapply(sets[usable], one, numeric(1))
  }
  tibble::tibble(parameter = names(sets), evppi = v, group = names(sets) %in% names(groups),
                 method = method) |>
    dplyr::arrange(dplyr::desc(evppi))
}

#' PSA output as matrices for other packages
#'
#' Returns `cost`, `effect` (simulations x strategies) and `parameters`
#' (simulations x parameters) - the inputs expected by
#' `dampack::make_psa_obj()`, `BCEA::bcea()` (as `e`, `c`) and
#' `voi::evppi()` (as `outputs = list(c = , e = , k = )` / `inputs`).
#' @param res A `cea_result`.
#' @export
psa_matrices <- function(res) {
  if (is.null(res$psa)) stop("no PSA in this result")
  a <- res$analysis
  list(cost = t(res$psa[a$cost, , ]), effect = t(res$psa[a$effect, , ]),
       parameters = as.data.frame(res$psa_params), strategies = res$arms)
}

#' Hand PSA results to dampack
#'
#' @param res A `cea_result`.
#' @return A `dampack` PSA object (`dampack::make_psa_obj()`), ready for
#'   `dampack::ceac()`, `dampack::calc_evpi()`, `dampack::owsa()` etc.
#'   dampack converts strategy names with `make.names()`.
#' @export
to_dampack <- function(res) {
  if (!requireNamespace("dampack", quietly = TRUE)) stop("install the dampack package")
  m <- psa_matrices(res)
  suppressWarnings(dampack::make_psa_obj(cost = as.data.frame(m$cost), effectiveness = as.data.frame(m$effect),
                                         parameters = m$parameters, strategies = m$strategies))
}

#' Hand PSA results to BCEA
#'
#' @param res A `cea_result`.
#' @param k Willingness-to-pay grid (default: the analysis grid).
#' @param ... Passed to `BCEA::bcea()`.
#' @return A `bcea` object (reference = the analysis reference arm).
#' @export
to_bcea <- function(res, k = res$settings$wtp_grid, ...) {
  if (!requireNamespace("BCEA", quietly = TRUE)) stop("install the BCEA package")
  m <- psa_matrices(res)
  BCEA::bcea(eff = m$effect, cost = m$cost, ref = match(res$analysis$reference, m$strategies),
             interventions = m$strategies, k = k, ...)
}

#' Hand PSA results to voi
#'
#' @param res A `cea_result`.
#' @param wtp Willingness to pay (default: the analysis WTP).
#' @return `list(outputs, inputs)` for `voi::evppi(outputs, inputs, pars)`,
#'   `voi::evsi()` etc.
#' @export
to_voi <- function(res, wtp = res$settings$wtp) {
  m <- psa_matrices(res)
  list(outputs = list(c = m$cost, e = m$effect, k = wtp), inputs = m$parameters)
}
