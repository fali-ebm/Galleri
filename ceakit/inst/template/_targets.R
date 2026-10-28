# Optional reproducible pipeline with the targets package:
#   install.packages(c("targets", "tarchetypes")); targets::tar_make()
# Re-runs only what changed when inputs or model code change.
library(targets)
tar_option_set(packages = "ceakit")
list(
  tar_target(input_files, list.files("inputs", full.names = TRUE), format = "file"),
  tar_target(model_files, list.files("R", full.names = TRUE), format = "file"),
  tar_target(results, { input_files; model_files; ceakit::run_project(".", quiet = TRUE) }),
  tar_target(tests, { results; testthat::test_dir("tests/testthat", stop_on_failure = TRUE); TRUE }),
  tar_target(report, { tests; rmarkdown::render("report.Rmd", output_dir = "outputs", quiet = TRUE) },
             format = "file")
)
