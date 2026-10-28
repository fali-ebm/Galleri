#' Format numbers for tables
#' @param x Numeric.
#' @param d Decimal places.
#' @export
fmt <- function(x, d = 0) ifelse(is.na(x) | !is.finite(x), "-",
  formatC(x, format = "f", digits = d, big.mark = ","))

metric_scale <- function(metrics, cohort_n) ifelse(metrics$type == "count", cohort_n, 1)

# "Dominated" / "Cost-saving" / "-" (reference) for ICER-type cells
icer_label <- function(res, arm) {
  a <- res$analysis
  if (arm == a$reference) return("-")
  ic <- res$base["inc_cost", arm]; ie <- res$base["inc_effect", arm]
  if (ie <= 0 && ic >= 0) return("Dominated")
  if (ie > 0 && ic < 0) return("Cost-saving (dominant)")
  if (ie < 0 && ic < 0) return(sprintf("%s saved per %s lost", fmt(ic / ie), a$effect_label))
  NA
}

#' Results in long format (one row per metric x arm)
#' @param res A `cea_result`.
#' @export
results_long <- function(res) {
  mt <- analysis_metrics(res$analysis); mt <- mt[mt$key %in% rownames(res$base), ]
  sc <- metric_scale(mt, res$analysis$cohort_n)
  purrr::map_dfr(res$arms, function(a) tibble::tibble(
    section = mt$section, metric = trimws(mt$label), key = mt$key,
    unit = ifelse(mt$type == "count", res$analysis$unit, "per person / ratio"), arm = a,
    value = res$base[mt$key, a] * sc,
    lower_95 = if (is.null(res$lo)) NA_real_ else res$lo[mt$key, a] * sc,
    upper_95 = if (is.null(res$hi)) NA_real_ else res$hi[mt$key, a] * sc))
}

#' Results in wide format ("value (95% UI)" per arm)
#' @param res A `cea_result`.
#' @param ui Include uncertainty intervals.
#' @export
results_wide <- function(res, ui = TRUE) {
  mt <- analysis_metrics(res$analysis); mt <- mt[mt$key %in% rownames(res$base), ]
  sc <- metric_scale(mt, res$analysis$cohort_n)
  out <- tibble::tibble(Section = mt$section, Metric = mt$label,
                        Unit = ifelse(mt$type == "count", res$analysis$unit, "per person / ratio"))
  ratio_keys <- c("icer", mt$key[grepl("^cost_per_", mt$key)])
  for (a in res$arms) out[[a]] <- vapply(seq_len(nrow(mt)), function(i) {
    k <- mt$key[i]; dg <- mt$digits[i]
    if (k == "icer") { lab <- icer_label(res, a); if (!is.na(lab)) return(lab) }
    v <- res$base[k, a]
    if (is.na(v)) return("-")
    s <- fmt(v * sc[i], dg)
    if (ui && !is.null(res$lo) && !k %in% ratio_keys)
      s <- sprintf("%s (%s to %s)", s, fmt(res$lo[k, a] * sc[i], dg), fmt(res$hi[k, a] * sc[i], dg))
    s
  }, character(1))
  out
}

#' Print a compact console summary of a result
#' @param x A `cea_result`.
#' @param ... Unused.
#' @export
print.cea_result <- function(x, ...) {
  cat(sprintf("<cea_result> %s | %d arms | PSA runs: %s | WTP %s\n", x$analysis$title,
              length(x$arms), if (is.null(x$psa)) 0 else dim(x$psa)[3], fmt(x$settings$wtp)))
  inc <- x$incremental
  inc$cost <- fmt(inc$cost, 2); inc$effect <- fmt(inc$effect, 4); inc$icer <- fmt(inc$icer)
  print(as.data.frame(inc), row.names = FALSE)
  if (!is.null(x$evpi)) cat(sprintf("EVPI per person at WTP: %s\n",
    fmt(x$evpi$evpi[which.min(abs(x$evpi$wtp - x$settings$wtp))], 2)))
  invisible(x)
}
