# =============================================================
# Update the model from the inputs workbook.
#
#   Rscript import_inputs_xlsx.R --check     validate and list every change; writes nothing
#   Rscript import_inputs_xlsx.R             validate, back up, write the CSVs, rebuild baseline
#   Rscript import_inputs_xlsx.R my.xlsx     another workbook (default Galleri_model_inputs.xlsx)
#
# Run from galleri_cea/. Create the workbook with
#   python3 tools/export_inputs_xlsx.py
# then run the model with Rscript run.R. Details: tools/import_inputs.R.
# =============================================================
source(file.path("tools", "import_inputs.R"))
args <- commandArgs(trailingOnly = TRUE)
check <- "--check" %in% args
xlsx <- setdiff(args, "--check"); xlsx <- if (length(xlsx)) xlsx[1] else "Galleri_model_inputs.xlsx"
if (!file.exists(xlsx)) stop("workbook not found: ", xlsx, call. = FALSE)

r <- import_workbook(xlsx, ".", apply = !check)
for (w in r$warnings) message("WARNING: ", w)
if (!r$ok) {
  message("\nNothing was changed. Fix these in the workbook and run again:")
  for (e in r$errors) message("  - ", e)
  quit(status = 1)
}
if (!length(r$changes)) {
  message("The workbook matches the model inputs: nothing to change.")
} else {
  message(if (check) "Changes the import would make:" else "Changes made:")
  for (cc in r$changes) message("  ", cc)
  if (check) message("\nRun without --check to apply.") else message("\nNext: Rscript run.R")
}
