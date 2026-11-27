# =============================================================
# Reproduce every result, file and page of the Galleri model, in R only.
#
#   Rscript make_all.R            (from galleri_cea/; about 30-40 minutes)
#   Rscript make_all.R --quick    (skips the 1M-person microsimulation; reuses outputs/microsim)
#
# Steps (each can also be run on its own):
#   1. run.R                          cohort model: checks, base case, PSA, OWSA, scenarios,
#                                     subgroups, accuracy, tests, CHEERS report -> outputs/
#   2. tools/run_microsim.R 1e6 2026  hesim individual-level model -> outputs/microsim/
#   3. tools/calibrate_miscan_eac.R   MISCAN-ESO calibration targets and validation -> outputs/miscan_eac/
#   4. tools/viz_data.R               data for the story pages -> viz/galleri_viz_data.json
#   5. viz/build_*.R                  the three HTML pages in viz/
#   6. tools/export_inputs_xlsx.R     Galleri_model_inputs.xlsx from inputs/ (round trip checked)
#   7. session_info.txt               R and package versions used
# Random numbers: PSA seed in inputs/settings.csv (seed), microsimulation seed 2026.
# =============================================================
quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
step <- function(label, script, args = character(0)) {
  message(sprintf("\n==== %s: Rscript %s %s", label, script, paste(args, collapse = " ")))
  t0 <- Sys.time()
  st <- system2(file.path(R.home("bin"), "Rscript"), c(script, args))
  if (st != 0) stop(sprintf("%s failed (exit %d)", script, st))
  message(sprintf("     done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
step("1 cohort model", "run.R")
if (!quick) step("2 microsimulation", "tools/run_microsim.R", c("1000000", "2026"))
step("3 MISCAN-ESO calibration", "tools/calibrate_miscan_eac.R")
step("4 page data", "tools/viz_data.R")
step("5a burden page", "viz/build_burden.R")
step("5b comparison page", "viz/build_compare.R")
step("5c decision model page", "viz/build_model_page.R")
step("6 inputs workbook", "tools/export_inputs_xlsx.R")
step("6b workbook round trip", "import_inputs_xlsx.R", "--check")
pk <- c("ceakit", "hesim", "data.table", "openxlsx2", "readxl", "readr", "jsonlite", "dplyr", "tidyr", "purrr",
        "tibble", "testthat", "rmarkdown", "targets", "voi", "ggplot2")
ver <- vapply(pk, function(k) if (requireNamespace(k, quietly = TRUE)) as.character(utils::packageVersion(k)) else "not installed", "")
writeLines(c(capture.output(sessionInfo()), "", "Packages used by the pipeline:", sprintf("  %-12s %s", pk, ver)), "session_info.txt")
message("\nAll done. R and package versions: session_info.txt")
