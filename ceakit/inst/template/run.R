# =============================================================
# Run the whole project: check inputs -> analyses -> outputs -> report
#   Rscript run.R          (or source("run.R") in RStudio)
# =============================================================
library(ceakit)

proj <- load_project(".")
print(check_project(proj))                 # errors stop the run; warnings list DUMMY inputs

results <- run_project(proj, allow_dummy = TRUE)   # set FALSE for final runs
for (r in results) print(r)

# Verification tests (tests/testthat) - must pass before reporting
testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)

# Report (HTML; change output_format to "word_document" for Word)
rmarkdown::render("report.Rmd", output_dir = "outputs", quiet = TRUE)
message("Report: outputs/report.html")
