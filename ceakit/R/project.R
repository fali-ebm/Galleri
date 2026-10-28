#' Create a new decision-modelling project from the template
#'
#' Copies the standard project structure (inputs, model, pipeline, tests,
#' report, app) into `path`.
#' @param path Folder to create (must not exist or be empty).
#' @param title Project title used in the README and report.
#' @return The path, invisibly.
#' @export
new_project <- function(path, title = basename(path)) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)))
    stop("folder exists and is not empty: ", path, call. = FALSE)
  tpl <- system.file("template", package = "ceakit")
  if (!nzchar(tpl)) stop("template not found - reinstall ceakit", call. = FALSE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(tpl, recursive = TRUE, all.files = TRUE)
  for (f in files) {
    dest <- file.path(path, sub("^dot-", ".", f))
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(tpl, f), dest, overwrite = TRUE)
    if (grepl("\\.(md|Rmd|R|txt)$", dest)) {
      x <- readLines(dest, warn = FALSE)
      writeLines(gsub("{{title}}", title, x, fixed = TRUE), dest)
    }
  }
  message("Created project in ", normalizePath(path), "\nNext: edit inputs/parameters.csv and R/model.R, then run run.R")
  invisible(path)
}

find_input <- function(dir, stem) {
  f <- list.files(dir, pattern = paste0("^", stem, "\\.(csv|xlsx|xls)$"), full.names = TRUE)
  if (length(f)) f[1] else NA_character_
}

#' Load a project
#'
#' Sources the project's `R/` files, reads `inputs/parameters.*` (applying
#' `inputs/counts.*` if present) and the project's data tables.
#'
#' The project's R files must define `analyses` (a list of [cea_analysis()]
#' objects). They may define `load_data(dir)` (returns a named list of data
#' tables; default: every other table in `inputs/`), `validate_data(data,
#' settings)` (returns `list(errors, warnings)`), `settings` (from
#' [default_settings()]) and `constraints`.
#' @param dir Project folder.
#' @return Object of class `cea_project`.
#' @export
load_project <- function(dir = ".") {
  dir <- normalizePath(dir, mustWork = TRUE)
  env <- new.env(parent = globalenv())
  for (f in sort(list.files(file.path(dir, "R"), pattern = "\\.R$", full.names = TRUE))) sys.source(f, envir = env)
  if (!exists("analyses", envir = env, inherits = FALSE))
    stop("the project's R/ files must define `analyses`", call. = FALSE)
  inputs <- file.path(dir, "inputs")
  pf <- find_input(inputs, "parameters"); if (is.na(pf)) stop("inputs/parameters.csv not found", call. = FALSE)
  params_raw <- read_parameters(pf)
  cf <- find_input(inputs, "counts")
  counts <- if (is.na(cf)) NULL else read_table(cf)
  data <- if (exists("load_data", envir = env, inherits = FALSE)) env$load_data(inputs) else {
    fs <- list.files(inputs, pattern = "\\.(csv|xlsx|xls)$", full.names = TRUE)
    fs <- fs[!tools::file_path_sans_ext(basename(fs)) %in% c("parameters", "counts", "scenarios")]
    stats::setNames(lapply(fs, read_table), tools::file_path_sans_ext(basename(fs)))
  }
  sf <- find_input(inputs, "scenarios")
  structure(list(dir = dir, env = env, analyses = env$analyses, params_raw = params_raw,
                 counts = counts, params = apply_counts(params_raw, counts), data = data,
                 scenarios = if (is.na(sf)) NULL else read_table(sf),
                 settings = if (exists("settings", envir = env, inherits = FALSE)) env$settings else default_settings()),
            class = "cea_project")
}

#' Check a project's inputs and model
#'
#' Validates parameters (against every analysis' `required`), runs the
#' project's `validate_data()` and evaluates each analysis once.
#' @param proj A `cea_project`.
#' @param params,data,settings Override the project's values (used by the app).
#' @return `list(errors, warnings)`.
#' @export
check_project <- function(proj, params = proj$params, data = proj$data, settings = proj$settings) {
  req <- unique(unlist(lapply(proj$analyses, `[[`, "required")))
  v <- validate_parameters(params, req)
  if (exists("validate_data", envir = proj$env, inherits = FALSE)) {
    d <- proj$env$validate_data(data, settings)
    if (is.character(d)) d <- list(errors = d, warnings = character(0))
    v$errors <- c(v$errors, d$errors); v$warnings <- c(v$warnings, d$warnings)
  }
  if (!length(v$errors)) for (a in proj$analyses) {
    err <- tryCatch({ evaluate(a, param_list(params), data, settings); NULL }, error = function(e) conditionMessage(e))
    if (!is.null(err)) v$errors <- c(v$errors, sprintf("analysis '%s' failed: %s", a$name, err))
  }
  v
}

#' Run every analysis in a project and write the outputs
#' @param proj A `cea_project` (or a folder path).
#' @param which Names of analyses to run (default: all).
#' @param out_dir Output folder (relative to the project).
#' @param allow_dummy Run even if DUMMY parameters remain.
#' @param settings Override settings.
#' @param quiet Suppress progress messages.
#' @return Named list of `cea_result`s (also saved as `outputs/results.rds`).
#' @export
run_project <- function(proj = ".", which = NULL, out_dir = "outputs", allow_dummy = TRUE,
                        settings = NULL, quiet = FALSE) {
  if (is.character(proj)) proj <- load_project(proj)
  if (!is.null(settings)) proj$settings <- settings
  v <- check_project(proj)
  if (length(v$errors)) stop("Input problems:\n - ", paste(v$errors, collapse = "\n - "), call. = FALSE)
  if (length(v$warnings)) { if (!allow_dummy) stop(v$warnings, call. = FALSE); if (!quiet) message("Note: ", v$warnings) }
  out <- file.path(proj$dir, out_dir)
  labels <- stats::setNames(dplyr::coalesce(proj$params$label, proj$params$parameter), proj$params$parameter)
  todo <- proj$analyses[names(proj$analyses) %in% (which %||% names(proj$analyses))]
  res <- lapply(todo, function(a) {
    if (!quiet) message("Running ", a$title, " ...")
    r <- run_cea(a, proj$params, proj$data, proj$settings,
                 progress = if (quiet) NULL else function(x) if (x %in% c(0.75, 1)) message(sprintf("  %3.0f%%", 100 * x)))
    write_cea(r, file.path(out, a$name), proj$data, labels)
    if (!is.null(proj$scenarios)) {
      sc <- run_scenarios(a, proj$params, proj$data, proj$settings, proj$scenarios)
      readr::write_csv(sc, file.path(out, a$name, "scenarios.csv")); r$scenarios <- sc
    }
    r
  })
  saveRDS(res, file.path(out, "results.rds"))
  if (!quiet) message("Done. Outputs in ", out)
  invisible(res)
}
