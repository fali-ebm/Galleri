#' Default colours for arms
#' @param arms Arm names.
#' @export
arm_colours <- function(arms) {
  pal <- c("#7F7F7F", "#4C78A8", "#F58518", "#E45756", "#54A24B", "#B279A2", "#72B7B2", "#9D755D")
  stats::setNames(rep(pal, length.out = length(arms)), arms)
}

#' ggplot2 theme used by ceakit
#' @export
theme_cea <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
                   plot.title.position = "plot", strip.text = ggplot2::element_text(face = "bold"))
}

cap1 <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))

#' Cost-effectiveness plane (PSA vs reference)
#' @param res A `cea_result`.
#' @param per Report per this many people.
#' @export
plot_ce_plane <- function(res, per = 1000) {
  if (is.null(res$psa)) return(NULL)
  a <- res$analysis; ref <- a$reference; others <- setdiff(res$arms, ref)
  sims <- purrr::map_dfr(others, function(s) tibble::tibble(arm = s,
    de = per * (res$psa[a$effect, s, ] - res$psa[a$effect, ref, ]),
    dc = per * (res$psa[a$cost, s, ] - res$psa[a$cost, ref, ]) / 1000))
  base <- tibble::tibble(arm = others,
    de = per * (res$base[a$effect, others] - res$base[a$effect, ref]),
    dc = per * (res$base[a$cost, others] - res$base[a$cost, ref]) / 1000)
  cols <- arm_colours(res$arms)
  ggplot2::ggplot(sims, ggplot2::aes(.data$de, .data$dc, colour = .data$arm)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") + ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_abline(slope = res$settings$wtp / 1000, intercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.25, size = 0.8) +
    ggplot2::geom_point(data = base, ggplot2::aes(fill = .data$arm), shape = 23, size = 4, colour = "black") +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(x = sprintf("Incremental %ss per %s (vs %s)", a$effect_label, scales::comma(per), ref),
                  y = sprintf("Incremental cost ($000 per %s)", scales::comma(per)),
                  title = "Cost-effectiveness plane (PSA)",
                  subtitle = sprintf("Dashed line: WTP %s per %s; diamonds: base case",
                                     scales::dollar(res$settings$wtp), a$effect_label)) + theme_cea()
}

#' Cost-effectiveness acceptability curves
#' @param res A `cea_result`.
#' @export
plot_ceac <- function(res) {
  if (is.null(res$ceac)) return(NULL)
  ggplot2::ggplot(res$ceac, ggplot2::aes(.data$wtp / 1000, .data$prob, colour = factor(.data$arm, res$arms))) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_vline(xintercept = res$settings$wtp / 1000, linetype = 2, colour = "grey50") +
    ggplot2::scale_colour_manual(values = arm_colours(res$arms), name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(x = sprintf("Willingness to pay ($000 per %s)", res$analysis$effect_label),
                  y = "Probability cost-effective", title = "Cost-effectiveness acceptability curves") + theme_cea()
}

#' EVPI curve
#' @param res A `cea_result`.
#' @export
plot_evpi <- function(res) {
  if (is.null(res$evpi)) return(NULL)
  ggplot2::ggplot(res$evpi, ggplot2::aes(.data$wtp / 1000, .data$evpi)) +
    ggplot2::geom_line(linewidth = 1.1, colour = "#4C78A8") +
    ggplot2::geom_vline(xintercept = res$settings$wtp / 1000, linetype = 2, colour = "grey50") +
    ggplot2::scale_y_continuous(labels = scales::dollar) +
    ggplot2::labs(x = sprintf("Willingness to pay ($000 per %s)", res$analysis$effect_label),
                  y = "EVPI per person", title = "Expected value of perfect information") + theme_cea()
}

#' EVPPI by parameter
#' @param res A `cea_result`.
#' @param labels Named vector of parameter labels.
#' @param top Number of parameters shown.
#' @export
plot_evppi <- function(res, labels = NULL, top = 15) {
  if (is.null(res$evppi)) return(NULL)
  d <- utils::head(res$evppi[!is.na(res$evppi$evppi), ], top)
  d$name <- if (is.null(labels)) d$parameter else dplyr::coalesce(labels[d$parameter], d$parameter)
  d$name <- factor(d$name, rev(unique(d$name)))
  ggplot2::ggplot(d, ggplot2::aes(.data$evppi, .data$name, fill = .data$group)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "#4C78A8", `TRUE` = "#F58518"),
                               labels = c(`FALSE` = "Single parameter", `TRUE` = "Parameter group"), name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::dollar) +
    ggplot2::labs(x = "EVPPI per person", y = NULL,
                  title = "Expected value of partial perfect information",
                  subtitle = sprintf("WTP %s per %s%s", scales::dollar(res$settings$wtp), res$analysis$effect_label,
                                     if (all(d$evppi < 1e-9)) " - all values ~0: no single parameter changes the decision" else "")) +
    theme_cea()
}

#' Tornado diagram from one-way sensitivity analysis
#' @param res A `cea_result`.
#' @param labels Named vector of parameter labels.
#' @export
plot_tornado <- function(res, labels = NULL) {
  ow <- res$owsa; if (is.null(ow) || !nrow(ow$top)) return(NULL)
  d <- ow$top
  d$name <- if (is.null(labels)) d$parameter else dplyr::coalesce(labels[d$parameter], d$parameter)
  d$name <- factor(d$name, rev(unique(d$name)))
  d <- tidyr::pivot_longer(d, c("nmb_at_low", "nmb_at_high"), names_to = "end", values_to = "nmb")
  d$end <- ifelse(d$end == "nmb_at_low", "Low value", "High value")
  ggplot2::ggplot(d, ggplot2::aes(y = .data$name, xmin = pmin(.data$nmb, ow$base),
                                  xmax = pmax(.data$nmb, ow$base), colour = .data$end)) +
    ggplot2::geom_linerange(linewidth = 7) +
    ggplot2::geom_vline(xintercept = ow$base, linetype = 2) + ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
    ggplot2::scale_colour_manual(values = c("Low value" = "#4C78A8", "High value" = "#F58518"), name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::dollar) +
    ggplot2::labs(y = NULL, x = sprintf("Incremental NMB per person: %s vs %s", ow$arms[1], ow$arms[2]),
                  title = "One-way sensitivity analysis",
                  subtitle = sprintf("Dashed line: base case; WTP %s", scales::dollar(res$settings$wtp))) + theme_cea()
}

#' Bar panels of selected metrics by arm (with 95% UI)
#' @param res A `cea_result`.
#' @param keys Metric keys (defaults to `analysis$key_metrics`).
#' @param title Plot title.
#' @export
plot_metrics <- function(res, keys = res$analysis$key_metrics, title = NULL) {
  if (is.null(keys)) keys <- c(res$analysis$cost, res$analysis$effect, "inc_nmb")
  d <- res$long[res$long$key %in% keys, ]
  d$metric <- cap1(sub("^- ", "", d$metric)); d$metric <- factor(d$metric, unique(d$metric))
  d$arm <- factor(d$arm, res$arms)
  ggplot2::ggplot(d, ggplot2::aes(.data$arm, .data$value, fill = .data$arm)) +
    ggplot2::geom_col(width = 0.7) +
    { if (!all(is.na(d$lower_95))) ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$lower_95, ymax = .data$upper_95), width = 0.2, colour = "grey30") } +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_text(ggplot2::aes(label = scales::comma(.data$value, accuracy = ifelse(abs(.data$value) < 10, 0.1, 1)),
                                    vjust = ifelse(.data$value < 0, 1.4, -0.4)), size = 3.2) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 4, labeller = ggplot2::label_wrap_gen(28)) +
    ggplot2::scale_fill_manual(values = arm_colours(res$arms), name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(c(0.08, 0.15))) +
    ggplot2::labs(x = NULL, y = res$analysis$unit, title = title %||% paste("Key results", res$analysis$unit),
                  subtitle = if (!is.null(res$lo)) "Bars: 95% uncertainty interval from the PSA") +
    theme_cea() + ggplot2::theme(axis.text.x = ggplot2::element_blank())
}

#' Draw a decision tree diagram (base graphics)
#' @param nodes Data frame `id`, `label`, `x`, `y`, `w` and optional `fill`, `tcol`.
#' @param edges Data frame `from`, `to`, `prob` (NA = no label).
#' @param main Title.
#' @param footer Footer text.
#' @param xlim,ylim Plot limits.
#' @param wrap Characters per unit width.
#' @param cex,pcex Text sizes for boxes and probabilities.
#' @export
draw_tree <- function(nodes, edges, main = "", footer = NULL, xlim = range(nodes$x) + c(-2, 2),
                      ylim = range(nodes$y) + c(-0.8, 0.8), wrap = 8.5, cex = 0.82, pcex = 0.7) {
  op <- graphics::par(mar = c(3, 0.5, 2.5, 0.5)); on.exit(graphics::par(op))
  graphics::plot.new(); graphics::plot.window(xlim = xlim, ylim = ylim)
  graphics::title(main, cex.main = 1.05)
  if (is.null(nodes$fill)) nodes$fill <- "#FCCD2A"
  if (is.null(nodes$tcol)) nodes$tcol <- "grey10"
  pos <- stats::setNames(seq_len(nrow(nodes)), nodes$id)
  wrapped <- lapply(seq_len(nrow(nodes)), function(i) strwrap(nodes$label[i], width = round(nodes$w[i] * wrap)))
  hgt <- vapply(wrapped, function(l) max(0.62, 0.3 * length(l) + 0.15), numeric(1))
  ends <- function(i) { a <- pos[edges$from[i]]; b <- pos[edges$to[i]]
    c(nodes$x[a], nodes$y[a] - hgt[a] / 2, nodes$x[b], nodes$y[b] + hgt[b] / 2) }
  for (i in seq_len(nrow(edges))) { e <- ends(i); graphics::arrows(e[1], e[2], e[3], e[4], length = 0.07, col = "grey40") }
  for (i in seq_len(nrow(nodes))) {
    n <- nodes[i, ]
    graphics::rect(n$x - n$w / 2, n$y - hgt[i] / 2, n$x + n$w / 2, n$y + hgt[i] / 2, col = n$fill, border = NA)
    graphics::text(n$x, n$y, paste(wrapped[[i]], collapse = "\n"), cex = cex, col = n$tcol)
  }
  for (i in which(!is.na(edges$prob))) {
    e <- ends(i); lab <- sprintf("p = %.3g", edges$prob[i])
    mx <- e[1] + 0.45 * (e[3] - e[1]); my <- e[2] + 0.45 * (e[4] - e[2])
    sw <- graphics::strwidth(lab, cex = pcex) / 2 + 0.08; sh <- graphics::strheight(lab, cex = pcex) / 2 + 0.05
    graphics::rect(mx - sw, my - sh, mx + sw, my + sh, col = "white", border = NA)
    graphics::text(mx, my, lab, cex = pcex, col = "#B03A00", font = 2)
  }
  if (!is.null(footer)) graphics::mtext(footer, side = 1, line = 1.2, cex = 0.78)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
