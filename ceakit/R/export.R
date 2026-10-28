#' Write all results of an analysis to a folder
#'
#' Writes CSV tables, an Excel workbook, PSA matrices, standard figures and any
#' project figures (`analysis$plots`).
#' @param res A `cea_result`.
#' @param dir Output folder.
#' @param data Data list (passed to project plot functions).
#' @param labels Named vector of parameter labels.
#' @export
write_cea <- function(res, dir, data = list(), labels = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(x) file.path(dir, x)
  if (is.null(labels)) labels <- stats::setNames(dplyr::coalesce(res$params$label, res$params$parameter),
                                                 res$params$parameter)
  tabs <- list(results_table = res$wide, results_long = res$long, incremental = res$incremental,
               ceac = res$ceac, evpi = res$evpi, evppi = res$evppi,
               owsa = if (!is.null(res$owsa)) res$owsa$table, parameters_used = res$params)
  tabs <- Filter(Negate(is.null), tabs)
  for (n in names(tabs)) readr::write_csv(tabs[[n]], out(paste0(n, ".csv")), na = "")
  if (!is.null(res$psa)) {
    m <- psa_matrices(res)
    readr::write_csv(as.data.frame(m$cost), out("psa_cost.csv"))
    readr::write_csv(as.data.frame(m$effect), out("psa_effect.csv"))
    readr::write_csv(m$parameters, out("psa_parameters.csv"))
  }
  openxlsx::write.xlsx(tabs, out("results.xlsx"), overwrite = TRUE)

  save_gg <- function(p, f, w = 10, h = 6.5) if (!is.null(p) && inherits(p, "ggplot"))
    ggplot2::ggsave(out(f), p, width = w, height = h, dpi = 150, bg = "white")
  save_gg(plot_metrics(res), "key_metrics.png", 13, 7)
  save_gg(plot_ce_plane(res), "ce_plane.png", 9, 7)
  save_gg(plot_ceac(res), "ceac.png", 9, 6)
  save_gg(plot_evpi(res), "evpi.png", 9, 6)
  save_gg(plot_evppi(res, labels), "evppi.png", 10, 7)
  save_gg(plot_tornado(res, labels), "tornado.png", 10, 7)
  for (n in names(res$analysis$plots)) {
    f <- out(paste0(n, ".png"))
    grDevices::png(f, width = 2400, height = 1300, res = 170)
    p <- tryCatch(res$analysis$plots[[n]](res, res$p, data), error = function(e) NULL)
    if (inherits(p, "ggplot")) print(p)
    grDevices::dev.off()
  }
  invisible(dir)
}
