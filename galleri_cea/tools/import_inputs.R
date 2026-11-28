# =============================================================
# Import the inputs workbook (tools/export_inputs_xlsx.R) back into the CSVs
# the model reads. Used by import_inputs_xlsx.R and the test suite.
#
#   import_workbook(xlsx, project, apply = FALSE)
#
# 1. Read every input sheet by its column headers (layout changes elsewhere on
#    a sheet are fine; renamed headers or sheets are an error).
# 2. Validate: types, ranges, keys, completeness of the raw tables.
# 3. Scratch copy of R/, inputs/, data-raw/: write the new CSVs, then run
#    data-raw/inflate_costs.R and data-raw/make_us_inputs.R, load the project,
#    run ceakit::check_project and a base-case calibration smoke test.
# 4. Diff the scratch files against the project (numeric comparison; cells that
#    did not change keep their original text) and report every change.
# 5. apply = TRUE: back up the files that change to
#    data-raw/import_backup/<timestamp>/ and copy the new ones in.
# The Baseline (calibration) sheet is derived and never imported.
# =============================================================

suppressPackageStartupMessages({ library(readxl); library(readr) })

RACE_CODES <- c(all = "All races/ethnicities", nhw = "Non-Hispanic White", nhb = "Non-Hispanic Black",
                nhapi = "Non-Hispanic Asian/Pacific Islander",
                nhaian = "Non-Hispanic American Indian/Alaska Native", hisp = "Hispanic")
NONE_REF <- "none (assumption)"
DISTS <- c("fixed", "beta", "gamma", "normal", "lognormal", "dirichlet")
STATUSES <- c("SOURCED", "ASSUMPTION", "DATA", "DUMMY")

num_fmt <- function(x) ifelse(is.na(x), NA_character_, trimws(formatC(x, digits = 15, format = "fg")))
as_num <- function(x) suppressWarnings(as.numeric(x))
band_start <- function(x) as.character(as_num(sub("^\\s*(\\d+).*$", "\\1", x)))

# ---- Sheet specifications ------------------------------------------------------
# Parameters live on the numbered tree-step tabs ("1 Galleri test" ... "9 Utilities and DALYs"); every tab
# whose name starts with a number is read into inputs/parameters.csv ('In model' is informational).
PARAM_SPEC <- list(sheet_pattern = "^[0-9]+ ", key = "parameter",
  cols = c(Parameter = "parameter", Label = "label", Value = "value", Status = "status",
           `PSA distribution` = "dist", SE = "se", `n (beta/Dirichlet)` = "n", `Dirichlet group` = "group",
           `Low (OWSA)` = "low", `High (OWSA)` = "high", Unit = "unit", Description = "description",
           `Source / derivation` = "source", `Used in` = "analysis", `Refs (see References)` = "refs",
           `Price as published (USD)` = "price_value", `Price year` = "price_year", `In model` = "in_model"),
  num = c("value", "se", "n", "low", "high", "price_value", "price_year"),
  drop = "in_model")
param_rows <- function(d) !is.na(d$parameter) & grepl("^[A-Za-z][A-Za-z0-9_.]*$", d$parameter) &
  !(is.na(d$value) & is.na(d$dist) & is.na(d$status) & is.na(d$price_value))     # skips the notes under each table
# cols: workbook header -> CSV column. keep: rows to keep. num: numeric CSV columns.
# key: CSV columns identifying a row (for diffs). file: target, relative to project.
SPECS <- list(
  parameters = c(PARAM_SPEC, list(file = "inputs/parameters.csv", keep = param_rows)),
  settings = list(sheet = "Settings", file = "inputs/settings.csv", key = "setting",
    cols = c(Key = "setting", Value = "value", Area = "area", Setting = "label", Note = "note", Refs = "refs"),
    num = character(0), keep = function(d) !is.na(d$setting)),
  cohort = list(sheet = "Cohort", file = "inputs/cohort.csv", key = c("population", "risk_group", "sex", "race_eth"),
    cols = c(population = "population", race_eth = "race_eth", sex = "sex", start_age = "start_age",
             risk_group = "risk_group", share = "share", kunzmann_scoped = "kunzmann_scoped", status = "status",
             notes = "notes"),
    num = c("start_age", "share", "kunzmann_scoped"), keep = function(d) !is.na(d$risk_group)),
  scenarios = list(sheet = "Scenarios", file = "inputs/scenarios.csv", key = c("scenario", "parameter"),
    cols = c(Scenario = "scenario", Parameter = "parameter", `Scenario value` = "value", Refs = "refs"),
    num = "value", keep = function(d) !is.na(d$scenario) & !is.na(d$parameter)),
  counts = list(sheet = "Workup counts", file = "inputs/counts.csv", key = c("branch", "outcome"),
    cols = c(Branch = "branch", Outcome = "outcome", `Sets parameter` = "parameter", `n (patients)` = "n",
             `What to count` = "notes"),
    num = "n", keep = function(d) !is.na(d$branch) & !is.na(d$outcome)),
  references = list(sheet = "References", file = "inputs/references.csv", key = "ref",
    cols = c(`#` = "ref", Reference = "citation", `DOI / URL` = "url", PMID = "pmid"),
    num = c("ref", "pmid"), keep = function(d) !is.na(as_num(d$ref))),
  seer_inc = list(sheet = "SEER incidence", file = "data-raw/raw/seer_csr_2013_2017.csv",
    key = c("site", "race", "sex", "age_band"),
    cols = c(Site = "site", Race = "race", Sex = "sex", `Age band (start)` = "age_band", `Rate per 100,000` = "rate"),
    num = c("age_band", "rate"), keep = function(d) !is.na(d$sex) & !is.na(d$age_band),
    fix = function(d) { d$age_band <- band_start(d$age_band); d }),
  seer_aa = list(sheet = "SEER race ratios", file = "data-raw/raw/seer_csr_age_adjusted_by_race.csv",
    key = c("site", "race"),
    cols = c(Site = "site", `Race/ethnicity` = "race", `Male, age-adjusted per 100,000` = "male",
             `Female, age-adjusted per 100,000` = "female"),
    num = c("male", "female"), keep = function(d) !is.na(d$race) & !is.na(d$male)),
  seer_stage = list(sheet = "SEER stage", file = "data-raw/raw/seer_stage_by_race_sex.csv", key = c("site", "race", "sex"),
    cols = c(Site = "site", Race = "race", Sex = "sex", `Localized %` = "localized", `Regional %` = "regional",
             `Distant %` = "distant", `Unstaged %` = "unstaged"),
    num = c("localized", "regional", "distant", "unstaged"), keep = function(d) !is.na(d$sex)),
  eac = list(sheet = "EAC share", file = "data-raw/raw/seer_eac_share_by_race.csv", key = "race",
    cols = c(race_eth = "race", `Adenocarcinoma share of esophageal cancer` = "eac_share"),
    num = "eac_share", keep = function(d) d$race %in% names(RACE_CODES)),
  seer_mort = list(sheet = "SEER mortality", file = "data-raw/raw/seer_us_mortality_2013_2017.csv",
    key = c("site", "race", "sex", "age_band"),
    cols = c(Site = "site", Race = "race", Sex = "sex", `Age band (start)` = "age_band",
             `US death rate per 100,000` = "rate"),
    num = c("age_band", "rate"), keep = function(d) !is.na(d$sex) & !is.na(d$age_band),
    fix = function(d) { d$age_band <- band_start(d$age_band); d }),
  seer_mort_aa = list(sheet = "SEER mortality ratios", file = "data-raw/raw/seer_us_mortality_age_adjusted_by_race.csv",
    key = c("site", "race"),
    cols = c(Site = "site", `Race/ethnicity` = "race", `Male, age-adjusted death rate per 100,000` = "male",
             `Female, age-adjusted death rate per 100,000` = "female"),
    num = c("male", "female"), keep = function(d) !is.na(d$race) & !is.na(d$male)),
  cpi = list(sheet = "CPI medical care", file = "data-raw/raw/cpi_medical_care.csv", key = "year",
    cols = c(Year = "year", `Annual average (index)` = "index", `First half` = "h1", `Second half` = "h2"),
    num = c("year", "index", "h1", "h2"), keep = function(d) !is.na(as_num(d$year))),
  eq5d = list(sheet = "EQ-5D norms", file = "data-raw/raw/eq5d_us_norms.csv", key = "age_from",
    cols = c(`Age from` = "age_from", `Age to` = "age_to", `EQ-5D-5L utility` = "utility"),
    num = c("age_from", "age_to", "utility"), keep = function(d) !is.na(as_num(d$age_from))),
  yll = list(sheet = "WHO standard YLL", file = "data-raw/raw/who_ghe_standard_yll.csv", key = "age",
    cols = c(`Age (band start; 85 = 85+)` = "age", `Standard YLL per death` = "yll"),
    num = c("age", "yll"), keep = function(d) !is.na(as_num(d$age)))
)

# ---- Reading -------------------------------------------------------------------
read_one <- function(xlsx, sheet, sp, errs) {
  d <- read_excel(xlsx, sheet, col_types = "text", .name_repair = "minimal")
  miss <- setdiff(names(sp$cols), names(d))
  if (length(miss)) {
    errs$add(sprintf("%s: column header(s) not found: %s", sheet, paste(sQuote(miss, FALSE), collapse = ", ")))
    return(NULL)
  }
  d <- as.data.frame(d[names(sp$cols)], stringsAsFactors = FALSE); names(d) <- unname(sp$cols); d
}

read_spec <- function(xlsx, sp, errs) {
  if (!is.null(sp$sheet_pattern)) {
    sh <- grep(sp$sheet_pattern, excel_sheets(xlsx), value = TRUE)
    if (!length(sh)) { errs$add("no numbered tree-step tabs (e.g. '1 Galleri test') found"); return(NULL) }
    parts <- lapply(sh, read_one, xlsx = xlsx, sp = sp, errs = errs)
    if (any(vapply(parts, is.null, TRUE))) return(NULL)
    d <- do.call(rbind, parts); sp$sheet <- "Tree-step tabs"
  } else {
    if (!sp$sheet %in% excel_sheets(xlsx)) { errs$add(sprintf("sheet '%s' is missing", sp$sheet)); return(NULL) }
    d <- read_one(xlsx, sp$sheet, sp, errs); if (is.null(d)) return(NULL)
  }
  d[] <- lapply(d, function(x) { x <- trimws(x); x[x == ""] <- NA; x })
  d <- d[sp$keep(d), , drop = FALSE]
  if (!is.null(sp$drop)) d <- d[setdiff(names(d), sp$drop)]
  if (!is.null(sp$fix)) d <- sp$fix(d)
  for (k in sp$num) {
    bad <- !is.na(d[[k]]) & is.na(as_num(d[[k]]))
    if (any(bad)) errs$add(sprintf("%s: non-numeric %s in row(s) for %s", sp$sheet, k,
                                   paste(utils::head(apply(d[bad, sp$key, drop = FALSE], 1, paste, collapse = "/"), 5), collapse = ", ")))
    d[[k]] <- num_fmt(as_num(d[[k]]))
  }
  rownames(d) <- NULL
  d
}

read_nchs <- function(xlsx, errs) {
  if (!"NCHS life tables" %in% excel_sheets(xlsx)) { errs$add("sheet 'NCHS life tables' is missing"); return(NULL) }
  d <- read_excel(xlsx, "NCHS life tables", col_types = "text", .name_repair = "minimal")
  d <- d[!is.na(as_num(d$Age)), , drop = FALSE]
  labs <- setdiff(names(d), "Age"); out <- list()
  for (h in labs) {
    sx <- sub("^.* (male|female)$", "\\1", h); lab <- sub(" (male|female)$", "", h)
    code <- names(RACE_CODES)[match(lab, RACE_CODES)]
    if (is.na(code) || !sx %in% c("male", "female")) { errs$add(sprintf("NCHS life tables: unrecognised column '%s'", h)); next }
    q <- as_num(d[[h]])
    if (anyNA(q)) errs$add(sprintf("NCHS life tables: non-numeric qx in '%s'", h))
    out[[h]] <- data.frame(race_eth = code, sex = sx, age = num_fmt(as_num(d$Age)), qx = num_fmt(q), stringsAsFactors = FALSE)
  }
  do.call(rbind, unname(out))
}

# ---- Validation ----------------------------------------------------------------
collector <- function() { v <- character(0); list(add = function(x) v <<- c(v, x), get = function() v) }

validate_inputs <- function(X, errs, warns) {
  n <- function(x) as_num(x)
  P <- X$parameters
  if (!is.null(P)) {
    if (any(dup <- duplicated(P$parameter))) errs$add(paste("Parameters: duplicate", paste(unique(P$parameter[dup]), collapse = ", ")))
    for (i in seq_len(nrow(P))) {
      k <- P$parameter[i]; v <- n(P$value[i]); lo <- n(P$low[i]); hi <- n(P$high[i]); d <- tolower(ifelse(is.na(P$dist[i]), "", P$dist[i]))
      if (is.na(v) && is.na(P$price_value[i])) errs$add(sprintf("Parameters %s: Value is empty", k))
      if (!d %in% DISTS) errs$add(sprintf("Parameters %s: PSA distribution '%s' must be one of %s", k, P$dist[i], paste(DISTS, collapse = "/")))
      if (!toupper(ifelse(is.na(P$status[i]), "", P$status[i])) %in% STATUSES) errs$add(sprintf("Parameters %s: status must be one of %s", k, paste(STATUSES, collapse = "/")))
      if (identical(P$unit[i], "proportion"))
        for (z in c(v, lo, hi)) if (!is.na(z) && (z < 0 || z > 1)) { errs$add(sprintf("Parameters %s: proportions must be 0-1", k)); break }
      if (!is.na(v) && !is.na(lo) && lo > v) errs$add(sprintf("Parameters %s: Low (%s) is above Value (%s)", k, lo, v))
      if (!is.na(v) && !is.na(hi) && hi < v) errs$add(sprintf("Parameters %s: High (%s) is below Value (%s)", k, hi, v))
      if (d %in% c("gamma", "normal", "lognormal") && !isTRUE(n(P$se[i]) > 0)) errs$add(sprintf("Parameters %s: %s needs SE > 0", k, d))
      if (d == "beta" && !isTRUE(n(P$n[i]) > 0) && !isTRUE(n(P$se[i]) > 0)) errs$add(sprintf("Parameters %s: beta needs n or SE > 0", k))
      if (d == "dirichlet" && (is.na(P$group[i]) || !isTRUE(n(P$n[i]) > 0))) errs$add(sprintf("Parameters %s: dirichlet needs a group and n > 0", k))
      if (xor(is.na(P$price_value[i]), is.na(P$price_year[i]))) errs$add(sprintf("Parameters %s: fill both Price as published and Price year, or neither", k))
    }
    P$analysis <- c(`Both analyses` = "both", `Comparison (Kunzmann arm)` = "comparison", `Net burden analysis` = "burden")[P$analysis] %|na|% tolower(P$analysis)
    P$refs[!is.na(P$refs) & P$refs == NONE_REF] <- NA
    X$parameters <- P
  }
  S <- X$settings
  if (!is.null(S)) {
    g <- function(k) n(S$value[S$setting == k])
    rng <- list(wtp = c(0, Inf), disc_cost = c(0, 0.2), disc_qaly = c(0, 0.2), disc_daly = c(0, 0.2), max_age = c(51, 100),
                screen_stop_age = c(50, 99), screen_interval = c(1, 10), kunz_share_male = c(0.001, 0.999),
                kunz_share_female = c(0.001, 0.999), kunz_capture_male = c(0.001, 0.999), kunz_capture_female = c(0.001, 0.999),
                kunz_egd_rounds = c(1, 50), n_sim = c(1, 1e6), seed = c(-Inf, Inf), owsa_pct = c(0.001, 0.99), cost_year = c(2000, 2100))
    for (k in names(rng)) {
      v <- g(k)
      if (!length(v)) next
      if (is.na(v) || v < rng[[k]][1] || v > rng[[k]][2]) errs$add(sprintf("Settings %s: '%s' must be a number in %s-%s", k, S$value[S$setting == k], rng[[k]][1], rng[[k]][2]))
    }
    for (k in c("n_sim", "screen_interval", "max_age", "screen_stop_age", "cost_year"))
      if (length(g(k)) && !is.na(g(k)) && g(k) != round(g(k))) errs$add(sprintf("Settings %s must be a whole number", k))
    lst <- function(k) trimws(strsplit(S$value[S$setting == k] %||% "", ",")[[1]])
    if (any(is.na(n(lst("subgroup_ages"))) | n(lst("subgroup_ages")) < 50 | n(lst("subgroup_ages")) > 99))
      errs$add("Settings subgroup_ages: comma-separated ages 50-99")
    if (!all(tolower(lst("subgroup_sexes")) %in% c("male", "female"))) errs$add("Settings subgroup_sexes: male and/or female")
    if (!all(tolower(lst("subgroup_races")) %in% names(RACE_CODES)))
      errs$add(paste("Settings subgroup_races: codes from", paste(names(RACE_CODES), collapse = ", ")))
    for (sx in c("male", "female")) {
      sh <- g(paste0("kunz_share_", sx)); ca <- g(paste0("kunz_capture_", sx))
      if (length(sh) && length(ca) && !anyNA(c(sh, ca)) && ca < sh)
        warns$add(sprintf("Settings: kunz_capture_%s < kunz_share_%s means the flagged group is at LOWER risk", sx, sx))
    }
  }
  C <- X$cohort
  if (!is.null(C) && !is.null(S)) {
    for (i in seq_len(nrow(C))) {                       # Kunzmann shares follow Settings
      k <- n(S$value[S$setting == paste0("kunz_share_", C$sex[i])])
      if (C$risk_group[i] %in% c("kunzmann_high", "kunzmann_low") && length(k) && !is.na(k))
        C$share[i] <- num_fmt(if (C$risk_group[i] == "kunzmann_high") k else 1 - k)
    }
    if (!all(C$race_eth %in% names(RACE_CODES))) errs$add("Cohort: race_eth must be a race/ethnicity code")
    if (!all(C$sex %in% c("male", "female"))) errs$add("Cohort: sex must be male or female")
    a <- n(C$start_age); if (any(is.na(a) | a < 50 | a > 99)) errs$add("Cohort: start_age must be 50-99")
    for (pp in unique(C$population)) if (abs(sum(n(C$share[C$population == pp])) - 1) > 1e-6)
      errs$add(sprintf("Cohort: shares for '%s' must sum to 1", pp))
    X$cohort <- C
  }
  if (!is.null(X$scenarios)) X$scenarios$refs[!is.na(X$scenarios$refs) & X$scenarios$refs == NONE_REF] <- NA
  if (!is.null(X$scenarios) && !is.null(P)) {
    bad <- setdiff(X$scenarios$parameter, P$parameter)
    if (length(bad)) errs$add(paste("Scenarios: unknown parameter(s)", paste(bad, collapse = ", ")))
    if (anyNA(n(X$scenarios$value))) errs$add("Scenarios: every scenario row needs a numeric Scenario value")
  }
  if (!is.null(X$counts)) {
    v <- n(X$counts$n)
    if (any(!is.na(v) & (v < 0 | v != round(v)))) errs$add("Workup counts: n must be a whole number >= 0")
  }
  R <- X$references
  if (!is.null(R)) {
    if (any(duplicated(R$ref))) errs$add("References: duplicate reference numbers")
    if (!is.null(P)) {
      used <- unique(trimws(unlist(strsplit(stats::na.omit(c(P$refs, X$scenarios$refs, S$refs)), "[,;]"))))
      used <- used[used != "" & !grepl("-", used)]
      miss <- setdiff(used, R$ref)
      if (length(miss)) errs$add(paste("Refs columns cite reference numbers not on the References sheet:", paste(miss, collapse = ", ")))
    }
  }
  pos <- function(sheet, d, cols) for (k in cols) if (any(is.na(n(d[[k]])) | n(d[[k]]) <= 0))
    errs$add(sprintf("%s: %s must be positive numbers", sheet, k))
  need <- function(sheet, d, cols, want) {
    have <- do.call(paste, c(d[cols], sep = "/")); miss <- setdiff(want, have)
    if (length(miss)) errs$add(sprintf("%s: missing rows %s", sheet, paste(utils::head(miss, 6), collapse = ", ")))
    if (any(duplicated(have))) errs$add(sprintf("%s: duplicate rows", sheet))
  }
  g3 <- function(sites, races, sexes) as.vector(outer(outer(sites, races, paste, sep = "/"), sexes, paste, sep = "/"))
  if (!is.null(X$seer_inc)) { pos("SEER incidence", X$seer_inc, "rate")
    for (k in g3(c("esophagus", "stomach", "all_sites"), c("all", "white", "black"), c("male", "female")))
      if (sum(do.call(paste, c(X$seer_inc[c("site", "race", "sex")], sep = "/")) == k) < 2)
        errs$add(sprintf("SEER incidence: %s needs at least two age bands", k)) }
  aa_races <- c("all", "white", "nhw", "black", "api", "aian", "hispanic")
  if (!is.null(X$seer_aa)) { pos("SEER race ratios", X$seer_aa, c("male", "female"))
    need("SEER race ratios", X$seer_aa, c("site", "race"), as.vector(outer(c("esophagus", "stomach", "all_sites"), aa_races, paste, sep = "/"))) }
  if (!is.null(X$seer_stage)) {
    need("SEER stage", X$seer_stage, c("site", "race", "sex"), g3(c("esophagus", "stomach"), c("all", "white", "black"), c("male", "female")))
    st <- X$seer_stage; tot <- n(st$localized) + n(st$regional) + n(st$distant) + n(st$unstaged)
    if (any(is.na(tot) | n(st$localized) < 0 | n(st$regional) < 0 | n(st$distant) < 0)) errs$add("SEER stage: percentages must be >= 0")
    else if (any(abs(tot - 100) > 5)) warns$add("SEER stage: some rows do not sum to about 100%")
  }
  if (!is.null(X$eac)) { need("EAC share", X$eac, "race", names(RACE_CODES))
    e <- n(X$eac$eac_share); if (any(is.na(e) | e < 0 | e > 1)) errs$add("EAC share: must be 0-1") }
  if (!is.null(X$seer_mort)) { pos("SEER mortality", X$seer_mort, "rate")
    for (k in g3(c("esophagus", "stomach"), c("all", "white", "black"), c("male", "female")))
      if (sum(do.call(paste, c(X$seer_mort[c("site", "race", "sex")], sep = "/")) == k) < 2)
        errs$add(sprintf("SEER mortality: %s needs at least two age bands", k)) }
  if (!is.null(X$seer_mort_aa)) { pos("SEER mortality ratios", X$seer_mort_aa, c("male", "female"))
    need("SEER mortality ratios", X$seer_mort_aa, c("site", "race"), as.vector(outer(c("esophagus", "stomach"), aa_races, paste, sep = "/"))) }
  if (!is.null(X$cpi)) {
    pos("CPI medical care", X$cpi, "index")
    if (!is.null(P) && !is.null(S)) {
      yrs <- unique(stats::na.omit(c(n(P$price_year), n(S$value[S$setting == "cost_year"]))))
      miss <- setdiff(pmin(yrs, max(n(X$cpi$year))), n(X$cpi$year))
      if (length(miss)) errs$add(paste("CPI medical care: no index for year(s)", paste(miss, collapse = ", ")))
    }
  }
  if (!is.null(X$nchs)) {
    need("NCHS life tables", X$nchs, c("race_eth", "sex", "age"),
         as.vector(outer(outer(names(RACE_CODES), c("male", "female"), paste, sep = "/"), 50:99, paste, sep = "/")))
    q <- n(X$nchs$qx); if (any(is.na(q) | q <= 0 | q >= 1)) errs$add("NCHS life tables: qx must be between 0 and 1")
  }
  if (!is.null(X$eq5d)) {
    u <- n(X$eq5d$utility); if (any(is.na(u) | u <= 0 | u > 1)) errs$add("EQ-5D norms: utility must be in (0, 1]")
    cover <- unlist(Map(seq, n(X$eq5d$age_from), n(X$eq5d$age_to)))
    if (!setequal(cover, 50:99) || any(duplicated(cover))) errs$add("EQ-5D norms: bands must cover each age 50-99 exactly once")
  }
  if (!is.null(X$yll)) { pos("WHO standard YLL", X$yll, "yll")
    if (!all(seq(50, 85, 5) %in% n(X$yll$age))) errs$add("WHO standard YLL: need ages 50, 55, ..., 85") }
  X
}

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
`%|na|%` <- function(a, b) { a <- unname(a); a[is.na(a)] <- b[is.na(a)]; a }

# ---- Writing, diffing ------------------------------------------------------------
read_csv_chr <- function(f) if (file.exists(f)) as.data.frame(read_csv(f, col_types = cols(.default = "c"), na = ""),
                                                              stringsAsFactors = FALSE) else NULL

# Align new to the old file's columns; keep the old text of cells that are numerically unchanged.
merge_old <- function(new, old, key) {
  if (is.null(old)) return(new)
  extra <- setdiff(names(old), names(new))
  for (k in extra) new[[k]] <- old[[k]][match(do.call(paste, new[key]), do.call(paste, old[key]))]
  new <- new[c(names(old), setdiff(names(new), names(old)))]
  m <- match(do.call(paste, new[key]), do.call(paste, old[key]))
  for (k in intersect(names(new), names(old))) {
    o <- old[[k]][m]; v <- new[[k]]
    same <- !is.na(m) & ((is.na(o) & is.na(v)) | (!is.na(o) & !is.na(v) & (o == v |
            (!is.na(as_num(o)) & !is.na(as_num(v)) & abs(as_num(o) - as_num(v)) <= 1e-9 * pmax(1, abs(as_num(o)))))))
    v[same] <- o[same]; new[[k]] <- v
  }
  new
}

diff_tables <- function(old, new, key, label, max_lines = 40) {
  out <- character(0)
  if (is.null(old)) return(sprintf("%s: new file (%d rows)", label, nrow(new)))
  ko <- do.call(paste, c(old[key], sep = "/")); kn <- do.call(paste, c(new[key], sep = "/"))
  for (k in setdiff(kn, ko)) out <- c(out, sprintf("%s: + row %s", label, k))
  for (k in setdiff(ko, kn)) out <- c(out, sprintf("%s: - row %s", label, k))
  for (k in intersect(kn, ko)) {
    a <- old[ko == k, , drop = FALSE][1, ]; b <- new[kn == k, , drop = FALSE][1, ]
    for (cc in setdiff(intersect(names(a), names(b)), key)) {
      x <- a[[cc]]; y <- b[[cc]]
      if (identical(x, y) || (is.na(x) && is.na(y))) next
      if (!is.na(as_num(x)) && !is.na(as_num(y)) && abs(as_num(x) - as_num(y)) <= 1e-9 * max(1, abs(as_num(x)))) next
      out <- c(out, sprintf("%s: %s %s: %s -> %s", label, k, cc, ifelse(is.na(x), "(blank)", x), ifelse(is.na(y), "(blank)", y)))
    }
  }
  if (length(out) > max_lines) out <- c(utils::head(out, max_lines), sprintf("%s: ... and %d more changes", label, length(out) - max_lines))
  out
}

# ---- Main --------------------------------------------------------------------------
import_workbook <- function(xlsx, project = ".", apply = FALSE, quiet = FALSE) {
  say <- function(...) if (!quiet) message(...)
  project <- normalizePath(project, mustWork = TRUE); xlsx <- normalizePath(xlsx, mustWork = TRUE)
  errs <- collector(); warns <- collector()
  X <- lapply(SPECS, read_spec, xlsx = xlsx, errs = errs)
  X$nchs <- read_nchs(xlsx, errs)
  X <- validate_inputs(X, errs, warns)
  if (length(errs$get())) return(invisible(list(ok = FALSE, errors = errs$get(), warnings = warns$get(), changes = character(0))))

  # Baseline sheet is derived: flag edits
  B_old <- read_csv_chr(file.path(project, "inputs/baseline_risks.csv"))
  if ("Baseline (calibration)" %in% excel_sheets(xlsx) && !is.null(B_old)) {
    bw <- read_excel(xlsx, "Baseline (calibration)", col_types = "text", .name_repair = "minimal")
    bw <- bw[!is.na(as_num(bw$Age)), ]
    if (nrow(bw) == nrow(B_old)) {
      chk <- c(5, 8, 9)                                      # mort_other, inc_eso, inc_stom
      ref <- cbind(as_num(B_old$mort_other), as_num(B_old$inc_eso) * 1e5, as_num(B_old$inc_stom) * 1e5)
      got <- sapply(chk, function(j) as_num(bw[[j]]))
      nd <- sum(abs(got - ref) > 1e-6 * pmax(1, abs(ref)), na.rm = TRUE)
      if (nd) warns$add(sprintf(paste("Baseline (calibration): %d cell(s) differ from the model's baseline. That sheet is derived and never",
                                      "imported (edit the raw-data sheets instead); re-export the workbook to refresh it."), nd))
    }
  }

  # Scratch copy and pipeline
  scratch <- tempfile("galleri_import_"); dir.create(scratch)
  for (d in c("R", "inputs", "data-raw")) file.copy(file.path(project, d), scratch, recursive = TRUE)
  unlink(file.path(scratch, "data-raw", "import_backup"), recursive = TRUE)
  for (nm in names(SPECS)) if (!is.null(X[[nm]])) {
    sp <- SPECS[[nm]]; f <- file.path(scratch, sp$file)
    write_csv(merge_old(X[[nm]], read_csv_chr(file.path(project, sp$file)), sp$key), f, na = "")
  }
  write_csv(merge_old(X$nchs, read_csv_chr(file.path(project, "data-raw/raw/nchs_lifetables_2024_qx.csv")),
                      c("race_eth", "sex", "age")), file.path(scratch, "data-raw/raw/nchs_lifetables_2024_qx.csv"), na = "")
  rscript <- file.path(R.home("bin"), "Rscript")
  old_wd <- setwd(scratch); on.exit(setwd(old_wd), add = TRUE)
  for (sc in c("data-raw/inflate_costs.R", "data-raw/make_us_inputs.R")) {
    o <- suppressWarnings(system2(rscript, sc, stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(o, "status")) && attr(o, "status") != 0) errs$add(paste0(sc, " failed: ", paste(utils::tail(o, 5), collapse = " | ")))
  }
  setwd(old_wd)
  if (!length(errs$get())) {
    res <- tryCatch({
      pj <- ceakit::load_project(scratch); ck <- ceakit::check_project(pj)
      for (e in ck$errors) errs$add(paste("model check:", e)); for (w in ck$warnings) warns$add(paste("model check:", w))
      if (!length(ck$errors)) {
        p <- ceakit::param_list(pj$params)
        ct <- pj$env$calibration_table(p, pj$data, pj$settings, pj$env$ARMS[1])
        rel <- max(abs(ct$model_clinical - ct$seer) / pmax(ct$seer, 1e-12))
        if (rel > 1e-6 || any(ct$clipped > 0)) warns$add(sprintf("calibration: base case does not reproduce SEER exactly (max rel. error %.2g; clipped years %d)", rel, sum(ct$clipped > 0)))
        out <- pj$env$run_population(p, pj$data, pj$settings, pj$env$ARMS)$total
        if (any(!is.finite(out))) errs$add("model run: non-finite results in the base case")
      }
    }, error = function(e) errs$add(paste("model load/run failed:", conditionMessage(e))))
  }
  if (length(errs$get())) return(invisible(list(ok = FALSE, errors = errs$get(), warnings = warns$get(), changes = character(0))))

  # Diff scratch vs project
  files <- c(vapply(SPECS, `[[`, "", "file"), nchs = "data-raw/raw/nchs_lifetables_2024_qx.csv",
             baseline = "inputs/baseline_risks.csv")
  keys <- c(lapply(SPECS, `[[`, "key"), list(nchs = c("race_eth", "sex", "age"), baseline = c("race_eth", "sex", "age_from")))
  changes <- character(0); changed_files <- character(0)
  for (nm in names(files)) {
    o <- read_csv_chr(file.path(project, files[[nm]])); n_ <- read_csv_chr(file.path(scratch, files[[nm]]))
    d <- diff_tables(o, n_, keys[[nm]], basename(files[[nm]]), max_lines = if (nm == "baseline") 10 else 60)
    if (length(d)) { changes <- c(changes, d); changed_files <- c(changed_files, files[[nm]]) }
  }
  if (apply && length(changed_files)) {
    bk <- file.path(project, "data-raw", "import_backup", format(Sys.time(), "%Y%m%d-%H%M%S"))
    for (f in changed_files) {
      dir.create(dirname(file.path(bk, f)), recursive = TRUE, showWarnings = FALSE)
      if (file.exists(file.path(project, f))) file.copy(file.path(project, f), file.path(bk, f))
      file.copy(file.path(scratch, f), file.path(project, f), overwrite = TRUE)
    }
    say("Backed up the previous files to ", bk)
  }
  unlink(scratch, recursive = TRUE)
  invisible(list(ok = TRUE, errors = character(0), warnings = warns$get(), changes = changes, files = changed_files,
                 applied = apply && length(changed_files) > 0))
}
