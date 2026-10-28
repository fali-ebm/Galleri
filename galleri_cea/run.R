# =============================================================
# Run the whole project: check inputs -> analyses -> outputs -> report
#   Rscript run.R          (or source("run.R") in RStudio)
# =============================================================
library(ceakit)

proj <- load_project(".")
print(check_project(proj))                 # errors stop the run; warnings list DUMMY inputs

results <- run_project(proj, allow_dummy = TRUE)   # set FALSE for final runs
for (r in results) print(r)

# Galleri-specific extra output: results by subgroup (cohort row)
for (r in results) readr::write_csv(
  proj$env$subgroup_table(r$p, proj$data, proj$settings, r$analysis$name, r$analysis$cohort_n),
  file.path("outputs", r$analysis$name, "subgroups.csv"))

# Verification tests (tests/testthat) - must pass before reporting
testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)

# Report (HTML; change output_format to "word_document" for Word)
rmarkdown::render("report.Rmd", output_dir = "outputs", quiet = TRUE)
message("Report: outputs/report.html")
