#' Define an analysis
#'
#' An analysis is one decision problem: a set of arms (strategies) and a
#' `run` function that returns per-person outcomes for every arm.
#'
#' @param name Short id, e.g. `"base"`.
#' @param title Human-readable title.
#' @param arms Character vector of arm names (first = reference by default).
#' @param run `function(p, data, settings)` returning a numeric matrix with one
#'   column per arm (named as `arms`) and one row per metric (rownames = metric
#'   keys). Must include the `cost` and `effect` rows. Values are per person and
#'   additive (they can be summed across subgroups).
#' @param cost,effect Row names of the cost and effect metrics.
#' @param effect_label Label for the effect, e.g. `"QALY"`.
#' @param reference Reference arm for incremental metrics.
#' @param metrics Tibble describing metrics: `key`, `section`, `label`, `type`
#'   (`"count"` = scaled to `cohort_n`; `"person"` = per person; `"ratio"` =
#'   unscaled) and `digits`. Standard cost-effectiveness rows are appended.
#' @param derive Optional `function(m, settings)` returning extra rows computed
#'   from the metric matrix (e.g. ratios, differences vs reference).
#' @param required Parameter names this analysis needs (for validation).
#' @param cohort_n Counts are reported per this many people.
#' @param unit Text for count units, e.g. `"per 1,000 screened"`.
#' @param owsa Two arms compared in the tornado (`c(arm, comparator)`).
#' @param key_metrics Metric keys shown in summary charts.
#' @param constraints Optional function applied to sampled parameters.
#' @param plots Named list of `function(res, p, data)` returning a ggplot (or
#'   drawing with base graphics) - project-specific figures.
#' @return An object of class `cea_analysis`.
#' @export
cea_analysis <- function(name, arms, run, title = name, cost = "cost", effect = "qaly",
                         effect_label = "QALY", reference = arms[1], metrics = NULL,
                         derive = NULL, required = NULL, cohort_n = 1, unit = "per person",
                         owsa = c(arms[length(arms)], reference), key_metrics = NULL,
                         constraints = NULL, plots = list()) {
  stopifnot(is.function(run), reference %in% arms, all(owsa %in% arms))
  if (is.null(metrics)) metrics <- tibble::tibble(
    key = c(cost, effect), section = "Outcomes", label = c("Cost per person", paste(effect_label, "per person")),
    type = c("person", "person"), digits = c(2, 4))
  structure(list(name = name, title = title, arms = arms, run = run, cost = cost, effect = effect,
                 effect_label = effect_label, reference = reference, metrics = metrics,
                 derive = derive, required = required, cohort_n = cohort_n, unit = unit, owsa = owsa,
                 key_metrics = key_metrics, constraints = constraints, plots = plots),
            class = "cea_analysis")
}

#' Analysis settings
#' @param wtp Willingness to pay per unit of effect.
#' @param wtp_grid Grid for acceptability curves and EVPI.
#' @param n_sim PSA simulations (0 = none).
#' @param seed Random seed.
#' @param owsa_pct Range for one-way SA when `low`/`high` are blank.
#' @param evppi Compute EVPPI (needs the mgcv package).
#' @param evppi_groups Named list of parameter groups for joint EVPPI.
#' @param population Size of the population affected by the decision (for
#'   population EVPI); `NA` to skip.
#' @param ... Any project-specific settings (e.g. discount rates, horizon),
#'   passed to the analysis' `run` function.
#' @export
default_settings <- function(wtp = 100000, wtp_grid = seq(0, 3 * wtp, length.out = 61),
                             n_sim = 1000, seed = 2026, owsa_pct = 0.25, evppi = TRUE,
                             evppi_groups = NULL, population = NA, ...) {
  c(list(wtp = wtp, wtp_grid = wtp_grid, n_sim = n_sim, seed = seed, owsa_pct = owsa_pct,
         evppi = evppi, evppi_groups = evppi_groups, population = population), list(...))
}

std_metrics <- function(a) {
  tibble::tibble(
    key = c("inc_cost", "inc_effect", "icer", "nmb", "inc_nmb"),
    section = "Cost-effectiveness",
    label = c(paste("Incremental cost per person vs", a$reference),
              sprintf("Incremental %ss vs %s", a$effect_label, a$reference),
              sprintf("ICER vs %s ($ per %s)", a$reference, a$effect_label),
              "Net monetary benefit per person ($)",
              paste("Incremental net monetary benefit per person vs", a$reference)),
    type = c("person", "count", "ratio", "person", "person"), digits = c(2, 2, 0, 2, 2))
}

#' Metric table for an analysis (project metrics + standard rows)
#' @param analysis A `cea_analysis`.
#' @export
analysis_metrics <- function(analysis) {
  m <- dplyr::bind_rows(analysis$metrics, std_metrics(analysis))
  m[!duplicated(m$key), ]
}

#' Evaluate an analysis for one parameter set
#' @param analysis A `cea_analysis`.
#' @param p Named list of parameter values.
#' @param data Named list of data tables.
#' @param settings Settings list.
#' @return Metric matrix (metrics x arms), per person.
#' @export
evaluate <- function(analysis, p, data, settings) {
  m <- analysis$run(p, data, settings)
  if (!is.matrix(m)) stop("run() must return a matrix", call. = FALSE)
  if (!identical(colnames(m), analysis$arms)) stop("run() columns must be: ", paste(analysis$arms, collapse = ", "), call. = FALSE)
  miss <- setdiff(c(analysis$cost, analysis$effect), rownames(m))
  if (length(miss)) stop("run() is missing rows: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!is.null(analysis$derive)) m <- rbind(m, analysis$derive(m, settings))
  ref <- analysis$reference
  ic <- m[analysis$cost, ] - m[analysis$cost, ref]
  ie <- m[analysis$effect, ] - m[analysis$effect, ref]
  nmb <- settings$wtp * m[analysis$effect, ] - m[analysis$cost, ]
  m <- rbind(m, inc_cost = ic, inc_effect = ie,
             icer = ifelse(colnames(m) == ref, NA, ic / ie),
             nmb = nmb, inc_nmb = nmb - nmb[ref])
  keys <- analysis_metrics(analysis)$key
  m[c(intersect(keys, rownames(m)), setdiff(rownames(m), keys)), , drop = FALSE]
}

#' Run a complete cost-effectiveness analysis
#'
#' Base case, PSA, incremental analysis, acceptability curves, EVPI, EVPPI and
#' one-way sensitivity analysis.
#' @param analysis A `cea_analysis`.
#' @param params Parameter table.
#' @param data Named list of data tables.
#' @param settings From [default_settings()].
#' @param psa,owsa Run the PSA / one-way SA.
#' @param progress Optional `function(fraction)` for progress reporting.
#' @return Object of class `cea_result`.
#' @export
run_cea <- function(analysis, params, data = list(), settings = default_settings(),
                    psa = TRUE, owsa = TRUE, progress = NULL) {
  prog <- function(x) if (!is.null(progress)) progress(x)
  p0 <- param_list(params)
  base <- evaluate(analysis, p0, data, settings)
  res <- list(analysis = analysis, settings = settings, params = params, p = p0, base = base,
              arms = analysis$arms, time = Sys.time())
  res$incremental <- icer_frontier(base[analysis$cost, ], base[analysis$effect, ], analysis$arms)

  if (psa && settings$n_sim > 0) {
    set.seed(settings$seed)
    varying <- params$parameter[params$dist != "fixed"]
    sims <- vector("list", settings$n_sim); draws <- vector("list", settings$n_sim)
    for (i in seq_len(settings$n_sim)) {
      q <- sample_parameters(params, analysis$constraints)
      sims[[i]] <- tryCatch(evaluate(analysis, q, data, settings), error = function(e) NULL)
      draws[[i]] <- unlist(q[varying])
      if (i %% 25 == 0) prog(0.05 + 0.65 * i / settings$n_sim)
    }
    ok <- !vapply(sims, is.null, logical(1))
    if (sum(!ok)) warning(sum(!ok), " PSA runs failed and were dropped", call. = FALSE)
    res$psa <- simplify2array(sims[ok])
    res$psa_params <- do.call(rbind, draws[ok])
    res$lo <- apply(res$psa, c(1, 2), stats::quantile, 0.025, na.rm = TRUE)
    res$hi <- apply(res$psa, c(1, 2), stats::quantile, 0.975, na.rm = TRUE)
    cm <- t(res$psa[analysis$cost, , ]); em <- t(res$psa[analysis$effect, , ])
    res$ceac <- ceac(cm, em, settings$wtp_grid)
    res$evpi <- evpi(cm, em, settings$wtp_grid, settings$population)
    if (isTRUE(settings$evppi)) {
      prog(0.75)
      res$evppi <- tryCatch(evppi(cm, em, res$psa_params, settings$wtp, settings$evppi_groups),
                            error = function(e) { warning("EVPPI skipped: ", conditionMessage(e), call. = FALSE); NULL })
    }
  }
  if (owsa) { prog(0.85); res$owsa <- run_owsa(analysis, params, data, settings) }
  res$long <- results_long(res)
  res$wide <- results_wide(res)
  prog(1)
  class(res) <- "cea_result"
  res
}

#' One-way sensitivity analysis on incremental net monetary benefit
#' @inheritParams run_cea
#' @param arms Two arms: `c(arm, comparator)`; defaults to `analysis$owsa`.
#' @param top Number of parameters kept for the tornado.
#' @export
run_owsa <- function(analysis, params, data, settings, arms = analysis$owsa, top = 15) {
  p0 <- param_list(params)
  inc <- function(q) { m <- evaluate(analysis, q, data, settings); m["nmb", arms[1]] - m["nmb", arms[2]] }
  base <- inc(p0)
  rg <- owsa_ranges(params, settings$owsa_pct)
  if (!is.null(analysis$required)) rg <- rg[rg$parameter %in% analysis$required, ]
  tab <- purrr::pmap_dfr(rg, function(parameter, low, high) {
    v <- vapply(c(low, high), function(x) {
      q <- p0; q[[parameter]] <- x
      if (!is.null(analysis$constraints)) q <- analysis$constraints(q)
      tryCatch(inc(q), error = function(e) NA_real_)
    }, numeric(1))
    tibble::tibble(parameter = parameter, low_value = low, high_value = high,
                   nmb_at_low = v[1], nmb_at_high = v[2])
  })
  tab$range <- abs(tab$nmb_at_high - tab$nmb_at_low)
  tab <- tab[order(-tab$range), ]
  list(base = base, arms = arms, table = tab,
       top = utils::head(tab[!is.na(tab$range) & tab$range > 0, ], top))
}

#' Run deterministic scenarios
#'
#' @inheritParams run_cea
#' @param scenarios Data frame with `scenario`, `parameter`, `value` (one row
#'   per changed parameter; scenarios may change several parameters).
#' @return Tibble: one row per scenario x arm with cost, effect, incremental
#'   values, ICER and NMB.
#' @export
run_scenarios <- function(analysis, params, data, settings, scenarios) {
  p0 <- param_list(params)
  run_one <- function(changes) {
    q <- p0
    for (i in seq_len(nrow(changes))) {
      if (!changes$parameter[i] %in% names(q)) stop("scenario parameter not found: ", changes$parameter[i], call. = FALSE)
      q[[changes$parameter[i]]] <- changes$value[i]
    }
    m <- evaluate(analysis, q, data, settings)
    g <- function(k) unname(m[k, ])
    tibble::tibble(arm = colnames(m), cost = g(analysis$cost), effect = g(analysis$effect),
                   inc_cost = g("inc_cost"), inc_effect = g("inc_effect"), icer = g("icer"), nmb = g("nmb"))
  }
  base <- run_one(scenarios[0, ]) |> dplyr::mutate(scenario = "Base case", .before = 1)
  others <- purrr::map_dfr(unique(scenarios$scenario), function(s)
    run_one(scenarios[scenarios$scenario == s, ]) |> dplyr::mutate(scenario = s, .before = 1))
  dplyr::bind_rows(base, others)
}
