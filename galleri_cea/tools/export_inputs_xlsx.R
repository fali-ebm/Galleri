# Export every model input to one Excel workbook (the mirror of import_inputs_xlsx.R).
#
#   Rscript tools/export_inputs_xlsx.R [output.xlsx]     (run from galleri_cea/)
#
# Reads inputs/*.csv and data-raw/raw/*.csv only; default output Galleri_model_inputs.xlsx.
# Uses openxlsx2 (CRAN). Formula cells (Kunzmann relative risks, upper GI totals, late
# shares, Cohort links to Settings) are computed by Excel/LibreOffice when the file is
# opened; the importer never relies on them. Round trip check:
#   Rscript tools/export_inputs_xlsx.R && Rscript import_inputs_xlsx.R --check   (expects no changes)
suppressPackageStartupMessages(library(openxlsx2))
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args)) args[1] else "Galleri_model_inputs.xlsx"
rd <- function(f) utils::read.csv(f, colClasses = "character", check.names = FALSE, na.strings = character(0), encoding = "UTF-8")
NONE_REF <- "none (assumption)"

# ---- styles ------------------------------------------------------------------
FONT <- "Arial"
HEX <- function(h) wb_color(hex = paste0("FF", h))
STATUS_FILL <- c(SOURCED = "D9EAD3", ASSUMPTION = "FCE5CD", DATA = "CFE2F3", DUMMY = "F4CCCC")
font <- function(wb, sh, dims, color = "000000", size = 10, bold = FALSE, italic = FALSE, underline = FALSE)
  wb_add_font(wb, sh, dims, name = FONT, size = as.character(size), color = HEX(color),
              bold = if (bold) "1" else "", italic = if (italic) "1" else "", underline = if (underline) "single" else "")
fill <- function(wb, sh, dims, hex) wb_add_fill(wb, sh, dims, color = HEX(hex))
align <- function(wb, sh, dims, wrap = FALSE, vertical = "top", horizontal = NULL)
  wb_add_cell_style(wb, sh, dims, wrap_text = if (wrap) "1" else "0", vertical = vertical, horizontal = horizontal)
numfmt <- function(wb, sh, dims, f) wb_add_numfmt(wb, sh, dims, numfmt = f)
cell <- function(r, c) wb_dims(rows = r, cols = c)
put <- function(wb, sh, r, c, v) if (length(v) && !is.na(v) && !identical(v, "")) wb_add_data(wb, sh, v, dims = cell(r, c), col_names = FALSE) else wb
header <- function(wb, sh, ncol) {
  d <- wb_dims(rows = 1, cols = seq_len(ncol))
  wb <- font(wb, sh, d, "FFFFFF", bold = TRUE); wb <- fill(wb, sh, d, "1F3A5F")
  wb <- align(wb, sh, d, wrap = TRUE, vertical = "center"); wb_set_row_heights(wb, sh, 1, 30)
}

# A number written as a number when the text is one; integers stay integers; "" -> NA.
num <- function(x) {
  if (is.na(x) || x == "") return(NA)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) return(x)
  if (v == round(v) && !grepl("[.eE]", x)) return(as.integer(v))
  v
}
row_df <- function(vals) { d <- as.data.frame(lapply(vals, function(v) if (is.null(v)) NA else v), stringsAsFactors = FALSE); names(d) <- paste0("V", seq_along(vals)); d }

# Generic table sheet: header row, data rows (a data.frame with the final column names),
# base font, blue for value columns, wrap, number formats, widths, freeze, Excel table.
sheet <- function(wb, title, df, widths, fmts = list(), value_cols = integer(0), wrap_cols = integer(0),
                  freeze = c(2, 1), table = TRUE) {
  wb <- wb_add_worksheet(wb, title)
  nr <- nrow(df); nc <- ncol(df)
  if (table && nr) wb <- wb_add_data_table(wb, title, df, table_style = "TableStyleLight1", with_filter = TRUE,
                                          table_name = gsub("[ -]", "_", gsub("[()]", "", title)), na = NULL)
  else wb <- wb_add_data(wb, title, df, na = NULL)
  wb <- header(wb, title, nc)
  if (nr) {
    body <- wb_dims(rows = 2:(nr + 1), cols = seq_len(nc))
    wb <- font(wb, title, body); wb <- align(wb, title, body)
    for (c in value_cols) wb <- font(wb, title, wb_dims(rows = 2:(nr + 1), cols = c), "0000FF")
    for (c in wrap_cols) wb <- align(wb, title, wb_dims(rows = 2:(nr + 1), cols = c), wrap = TRUE)
    for (c in names(fmts)) wb <- numfmt(wb, title, wb_dims(rows = 2:(nr + 1), cols = as.integer(c)), fmts[[c]])
  }
  wb <- wb_set_col_widths(wb, title, cols = seq_along(widths), widths = widths)
  wb_freeze_pane(wb, title, first_active_row = freeze[1], first_active_col = freeze[2])
}
note <- function(wb, sh, r, txt, size = 9, italic = TRUE, bold = FALSE, color = "000000") {
  wb <- put(wb, sh, r, 1, txt); font(wb, sh, cell(r, 1), color, size = size, italic = italic, bold = bold)
}

wb <- wb_workbook()

# ---- README ------------------------------------------------------------------
lines <- list(
  c("Galleri upper GI screening model: all inputs", "title"),
  c("Source: github.com/fali-ebm/Galleri, folder galleri_cea. Exported from the files the model reads (tools/export_inputs_xlsx.R).", ""),
  c("", ""),
  c("How to update the model from this workbook", "h"),
  c("1. Edit blue cells (and add or delete rows on the numbered tree-step tabs, Scenarios, References and the raw-data sheets). Black cells are formulas or labels; the Baseline (calibration) sheet is rebuilt from the raw-data sheets, so edits there are ignored.", ""),
  c("2. Save the workbook as galleri_cea/Galleri_model_inputs.xlsx (keep the sheet names and column headers).", ""),
  c("3. From galleri_cea/, preview: Rscript import_inputs_xlsx.R --check   (validates and lists every change, writes nothing).", ""),
  c("4. Apply: Rscript import_inputs_xlsx.R   (validates on a scratch copy, backs up the old CSVs, writes the new ones, rebuilds the SEER/NCHS baseline).", ""),
  c("5. Run the model: Rscript run.R", ""),
  c("", ""),
  c("Model in one paragraph", "h"),
  c("Annual-cycle cohort model, lifetime horizon (to age 100), US payer perspective, 2025 USD, 3% discounting of costs and QALYs. Arms: usual care; Kunzmann score once -> one-time EGD if score >= 8; Galleri every year: +ve -> EGD (EGD +ve = true positive; EGD -ve = false positive -> rest of the workup; any EGD can perforate), -ve -> true negative (stays well) or false negative (symptomatic cancer later or found at a later round), then usual care. Screening runs from the start age to 85. Natural history is calibrated to SEER for every subgroup: onset of preclinical upper GI cancer is solved so that, without screening, symptomatic diagnoses reproduce SEER incidence by age, sex and race/ethnicity, and their stage reproduces SEER stage at diagnosis. Base case: 50-year-old non-Hispanic White men.", ""),
  c("Revised tree (Nov 2026, in the engine): a 2-year hold on repeat EGDs after a negative EGD; a negative EGD can hide a visible lesion that was missed or a cancer not yet visible; Bayesian update of each person's false-positive propensity; EUS sampling as a coin toss; tumor board -> laparoscopy or watchful waiting, with laparoscopy less likely after each negative result. Its inputs are on the tree-step tabs.", ""),
  c("Microsimulation (Nov 2026): the same inputs also drive an individual-level hesim model (R/microsim.R; run Rscript tools/run_microsim.R 1000000 2026). It follows each simulated person in continuous time and cross-checks the cohort model, which runs the PSA, subgroups and value of information.", ""),
  c("MISCAN-ESO scenario (Nov 2026): esophageal adenocarcinoma gets its own natural history calibrated to the CISNET MISCAN-ESO model (sex-specific preclinical sojourn, early -> late progression before symptoms); inputs on tab 6, calibration script tools/calibrate_miscan_eac.R.", ""),
  c("", ""),
  c("Sheets", "h"),
  c("1 Galleri test ... 10 Utilities and DALYs - every scalar input, one tab per step of the decision tree: value, one-way range, PSA distribution, status, whether the current model uses it (In model), and source.", ""),
  c("References - numbered sources; Refs columns elsewhere point to these numbers.", ""),
  c("Settings - analysis settings: willingness to pay, discounting, horizon, screening schedule, Kunzmann split, subgroup grid, PSA size.", ""),
  c("Cohort - base-case cohort rows (start age, sex, race/ethnicity). The Kunzmann shares are linked to Settings.", ""),
  c("Scenarios - deterministic scenario analyses (parameter overrides).", ""),
  c("Workup counts - template for local Galleri +ve workup counts; filled counts override the workup probabilities.", ""),
  c("Baseline (calibration) - life tables and SEER calibration targets by race/ethnicity x sex x single year of age 50-99 (600 rows).", ""),
  c("SEER incidence, SEER race ratios, SEER stage, EAC share, SEER mortality, SEER mortality ratios, NCHS life tables, EQ-5D norms, WHO standard YLL - the raw published tables the baseline is built from. CPI medical care - inflation index for costs.", ""),
  c("", ""),
  c("Status legend", "h"),
  c("SOURCED", "SOURCED"), c("ASSUMPTION", "ASSUMPTION"), c("DATA", "DATA"), c("DUMMY", "DUMMY"),
  c("", ""),
  c("Conventions", "h"),
  c("Blue text = editable input values; black = formulas, labels or derived values. Probabilities are stored as fractions and shown as %. Costs in USD (Medicare rates as payer proxy; years 2023-2026, not yet inflated to a common year). Rates in the baseline sheet are per person-year; per-100,000 columns are provided for reading.", ""),
  c("Yellow cells on the Workup counts sheet are the ones to fill in (n column).", ""),
  c("Full derivations and references: galleri_cea/inputs/SOURCES.md.", ""))
desc <- c(SOURCED = "Taken from, or derived arithmetically from, a published source.",
          ASSUMPTION = "No direct source; value and reasoning documented so it can be challenged.",
          DATA = "Set from local data (counts.csv).", DUMMY = "Placeholder with no source (none remain).")
wb <- wb_add_worksheet(wb, "README")
for (i in seq_along(lines)) {
  txt <- lines[[i]][1]; kind <- lines[[i]][2]
  if (nzchar(txt)) wb <- put(wb, "README", i, 1, txt)
  wb <- font(wb, "README", cell(i, 1), size = if (kind == "title") 14 else if (kind == "h") 11 else 10,
             bold = kind %in% c("title", "h", names(STATUS_FILL)))
  wb <- align(wb, "README", cell(i, 1), wrap = TRUE)
  if (kind %in% names(STATUS_FILL)) {
    wb <- fill(wb, "README", cell(i, 1), STATUS_FILL[[kind]])
    wb <- put(wb, "README", i, 2, desc[[kind]]); wb <- font(wb, "README", cell(i, 2))
  }
}
wb <- wb_set_col_widths(wb, "README", cols = 1:2, widths = c(110, 70))

# ---- Parameters: one tab per step of the decision tree -----------------------
# Tabs whose names start with a number are read by the importer into inputs/parameters.csv.
STEPS <- list(
  list("1 Galleri test", "Annual Galleri test. Chance of a positive result with an esophagus/stomach (E/S) signal for each true health state, how sensitivity falls after round 1, and whether false positives repeat in the same person.",
       c("sens_g_ugi", "sens_g_other", "fpr_g_es", "sens_g_rep_ratio", "p_fn_persist", "fp_kappa", "hr_ca_after_fp")),
  list("2 EGD", "Galleri +ve: EGD unless inside the hold after a negative EGD. What the EGD finds, what it misses (visible lesion missed vs cancer not yet visible), and perforation.",
       c("sens_egd", "f_invisible", "t_visible", "egd_hold_years", "p_perf_egd", "p_death_perf")),
  list("3 ENT and PET CT", "EGD -ve or skipped: routed to ENT exam or PET CT. What each finds.",
       c("p_ent_route", "p_ent_pos", "p_pet_pb", "p_pet_incid", "p_pet_ugi", "p_incid_cancer", "p_miss", "p_egd_pos")),
  list("4 Pancreatobiliary", "PET pancreatobiliary +ve: EUS (FNA/FNB a coin toss), tumor board after a negative or unsampled EUS, laparoscopy or watchful waiting. After a prior negative EUS or laparoscopy, laparoscopy becomes less likely.",
       c("p_pb_malig", "p_eus_samp", "eus_tissue_sens", "p_tb_lap", "r_lap_after_eus", "r_lap_after_lap", "p_reeus_after_eus", "p_reeus_after_lap")),
  list("5 Kunzmann arm", "Kunzmann score once at the start age; a one-time EGD for a score of 8 or more. The flagged share and sensitivity are on the Settings tab; this switch runs the sex-specific split as a scenario.",
       c("kunz_sexspec")),
  list("6 Natural history", "Preclinical sojourn, stage at screen detection and at symptoms, and cancer mortality by stage. The MISCAN-ESO rows give esophageal adenocarcinoma its own natural history in the miscan_eac scenario.",
       c("sojourn_years", "sojourn_eac", "site_sojourn", "p_late_screen", "or_late_clin", "m_early", "m_late",
         "miscan_eac", "sojourn_eac_male", "sojourn_eac_female", "eac_dwell_ratio", "m_early_eac", "m_late_eac")),
  list("7 Procedure harms", "Serious adverse events of each procedure and radiation from PET CT. EGD perforation is on the EGD tab.",
       c("sae_ent", "sae_eus", "sae_eus_samp", "p_lap_cx", "p_death_sae", "msv_pet", "risk_per_sv", "loss_rad_ca")),
  list("8 Time to each step", "Days from one workup step to the next (time spent waiting for the all-clear).",
       c("t_egd", "t_ent", "t_pet", "t_eus", "t_tb", "t_lap")),
  list("9 Costs", "Costs in USD, inflated to the Settings cost year. Rows with a published price: edit Price as published and Price year.",
       c("c_galleri", "galleri_cost_share", "c_kunz", "c_egd", "c_ent", "c_pet", "c_eus", "c_eus_samp", "c_tb", "c_lap", "c_surv",
         "c_incid_fu", "c_symp_dx", "c_tx_first_early", "c_tx_first_late", "c_fu_early", "c_fu_late", "c_eol", "c_sae", "c_perf")),
  list("10 Utilities and DALYs", "Quality-of-life multipliers, disutilities of procedures and false positives, and disability weights.",
       c("u_mult_early", "u_mult_late", "d_proc", "d_sae", "dfp_endo", "dfp_pet", "dfp_ext", "fp_tto", "dw_anx", "p_worry",
         "dur_worry", "dw_early", "dw_late")))
PR <- rd("inputs/parameters.csv"); P <- lapply(split(PR, PR$parameter), as.list)[PR$parameter]   # file order
for (k in setdiff(c("refs", "price_value", "price_year"), names(PR))) PR[[k]] <- ""
allk <- unlist(lapply(STEPS, `[[`, 3))
# STEPS keys present in parameters.csv; parameters missing from STEPS are an error
STEPS <- lapply(STEPS, function(s) { s[[3]] <- intersect(s[[3]], names(P)); s })
allk <- unlist(lapply(STEPS, `[[`, 3))
if (!setequal(allk, names(P)) || anyDuplicated(allk))
  stop("STEPS and parameters.csv differ: ", paste(union(setdiff(allk, names(P)), setdiff(names(P), allk)), collapse = ", "))
IN_MODEL_NOTE <- c(p_egd_pos = "Yes (validation target only)")
USED_IN <- c(both = "Both analyses", comparison = "Comparison (Kunzmann arm)")
hdr <- c("Parameter", "Label", "Value", "Low (OWSA)", "High (OWSA)", "PSA distribution", "SE", "n (beta/Dirichlet)",
         "Dirichlet group", "Status", "In model", "Unit", "Description", "Source / derivation", "Refs (see References)",
         "Used in", "Price as published (USD)", "Price year")
g0 <- function(p, k) { v <- p[[k]]; if (is.null(v) || is.na(v)) "" else v }
for (s in STEPS) {
  title <- s[[1]]; keys <- s[[3]]
  wb <- wb_add_worksheet(wb, title)
  wb <- wb_add_data(wb, title, as.data.frame(t(hdr)), col_names = FALSE)
  wb <- header(wb, title, length(hdr))
  for (j in seq_along(keys)) {
    k <- keys[j]; p <- P[[k]]; r <- j + 1
    vals <- list(k, g0(p, "label"), num(g0(p, "value")), num(g0(p, "low")), num(g0(p, "high")), g0(p, "dist"), num(g0(p, "se")),
                 num(g0(p, "n")), if (nzchar(g0(p, "group"))) g0(p, "group") else NA, g0(p, "status"),
                 if (k %in% names(IN_MODEL_NOTE)) IN_MODEL_NOTE[[k]] else "Yes", g0(p, "unit"), g0(p, "description"), g0(p, "source"),
                 if (nzchar(g0(p, "refs"))) g0(p, "refs") else NONE_REF,
                 if (g0(p, "analysis") %in% names(USED_IN)) USED_IN[[g0(p, "analysis")]] else g0(p, "analysis"),
                 num(g0(p, "price_value")), num(g0(p, "price_year")))
    for (c in seq_along(vals)) wb <- put(wb, title, r, c, vals[[c]])
    rowd <- wb_dims(rows = r, cols = 1:18)
    wb <- font(wb, title, rowd); wb <- align(wb, title, rowd)
    for (c in c(2, 13, 14)) wb <- align(wb, title, cell(r, c), wrap = TRUE)
    wb <- font(wb, title, wb_dims(rows = r, cols = 3:8), "0000FF")
    u <- g0(p, "unit"); v <- suppressWarnings(as.numeric(g0(p, "value"))); if (is.na(v)) v <- 0
    f <- if (u == "proportion" && v >= 0.01) "0.0%" else if (u %in% c("proportion", "per procedure", "per year")) "0.000%" else
         if (startsWith(u, "USD")) "$#,##0" else if (u == "count") "#,##0" else "0.00###"
    if (u == "per year" && startsWith(k, "m_")) f <- "0.0%"
    wb <- numfmt(wb, title, wb_dims(rows = r, cols = 3:5), f)
    if (nzchar(g0(p, "price_value"))) {
      wb <- font(wb, title, wb_dims(rows = r, cols = 17:18), "0000FF"); wb <- numfmt(wb, title, cell(r, 17), "$#,##0")
      wb <- font(wb, title, cell(r, 3))            # value is computed from the published price on import
    }
    st <- g0(p, "status"); wb <- fill(wb, title, cell(r, 10), if (st %in% names(STATUS_FILL)) STATUS_FILL[[st]] else "FFFFFF")
  }
  wb <- wb_set_col_widths(wb, title, cols = 1:18, widths = c(20, 34, 12, 11, 11, 13, 9, 10, 10, 13, 15, 12, 52, 60, 14, 14, 14, 10))
  wb <- wb_freeze_pane(wb, title, first_active_row = 2, first_active_col = 2)
  wb <- wb_add_filter(wb, title, rows = 1, cols = 1:18)
  n0 <- length(keys) + 3
  notes <- c(s[[2]],
    "Edit the blue cells: Value is the base case; Low/High are one-way sensitivity bounds (blank = +/-25%); PSA distribution with SE or n sets the probabilistic range (beta and Dirichlet use n as the effective sample size).",
    "Status: SOURCED = from a published source; ASSUMPTION = reasoned, no direct source; DUMMY = placeholder with no source yet (edit it).")
  for (j in seq_along(notes)) {
    r <- n0 + j - 1
    wb <- put(wb, title, r, 1, notes[j]); wb <- font(wb, title, cell(r, 1), bold = j == 1)
    wb <- wb_merge_cells(wb, title, wb_dims(rows = r, cols = 1:14)); wb <- align(wb, title, cell(r, 1), wrap = TRUE)
    wb <- wb_set_row_heights(wb, title, r, 30)
  }
}

# ---- Settings ------------------------------------------------------------------
ST <- rd("inputs/settings.csv")
FMT <- c(wtp = "$#,##0", disc_cost = "0.0%", disc_qaly = "0.0%", disc_daly = "0.0%", kunz_share_male = "0.0%",
         kunz_share_female = "0.0%", kunz_capture_male = "0.0%", kunz_capture_female = "0.0%", owsa_pct = "0%",
         n_sim = "#,##0", cost_year = "0", kunz_egd_rounds = "0")
sh <- "Settings"; wb <- wb_add_worksheet(wb, sh)
wb <- wb_add_data(wb, sh, as.data.frame(t(c("Area", "Setting", "Key", "Value", "Note", "Refs"))), col_names = FALSE)
wb <- header(wb, sh, 6)
row_of <- integer(0)
for (j in seq_len(nrow(ST))) {
  r <- j + 1; x <- ST[j, ]; row_of[x$setting] <- r
  vals <- list(x$area, x$label, x$setting, num(x$value), x$note, if (nzchar(x$refs)) x$refs else NA)
  for (c in 1:6) wb <- put(wb, sh, r, c, vals[[c]])
  wb <- font(wb, sh, wb_dims(rows = r, cols = 1:6)); wb <- align(wb, sh, wb_dims(rows = r, cols = 1:6))
  for (c in c(2, 5)) wb <- align(wb, sh, cell(r, c), wrap = TRUE)
  wb <- font(wb, sh, cell(r, 4), "0000FF"); wb <- align(wb, sh, cell(r, 4), wrap = TRUE, horizontal = "left")
  if (x$setting %in% names(FMT)) wb <- numfmt(wb, sh, cell(r, 4), FMT[[x$setting]])
}
D <- function(k) paste0("D", row_of[[k]])
lst <- function(k) sprintf('(LEN(%s)-LEN(SUBSTITUTE(%s,",",""))+1)', D(k), D(k))
derived <- list(
  list("Kunzmann score", "Implied relative risk, flagged, men", sprintf("=%s/%s", D("kunz_capture_male"), D("kunz_share_male")), "within-sex; applied to EAC share of esophageal cancer", "8, 9", "0.00"),
  list("Kunzmann score", "Implied relative risk, not flagged, men", sprintf("=(1-%s)/(1-%s)", D("kunz_capture_male"), D("kunz_share_male")), "within-sex", "8, 9", "0.00"),
  list("Kunzmann score", "Implied relative risk, flagged, women", sprintf("=%s/%s", D("kunz_capture_female"), D("kunz_share_female")), "within-sex", "8, 9", "0.00"),
  list("Kunzmann score", "Implied relative risk, not flagged, women", sprintf("=(1-%s)/(1-%s)", D("kunz_capture_female"), D("kunz_share_female")), "within-sex", "8, 9", "0.00"),
  list("Subgroup grid", "Number of subgroups", sprintf("=%s*%s*%s", lst("subgroup_ages"), lst("subgroup_sexes"), lst("subgroup_races")), "deterministic runs, each re-calibrated", NA, "0"),
  list("Model structure", "Cycle length", 1L, "year (fixed in model code)", NA, "0"),
  list("Model structure", "Half-cycle correction", "Yes", "mid-cycle alive fraction (fixed in model code)", NA, NA),
  list("Economic", "Perspective", "US payer", "Medicare rates as proxy (fixed)", NA, NA),
  list("Screening", "Galleri arm", "Galleri every round; +ve -> full workup tree", "fixed in model code", "1, 2", NA),
  list("Screening", "Kunzmann arm", "Score once at start; EGD for score >= 8 (one-time unless kunz_egd_rounds > 1)", "kunz_egd_rounds above", NA, NA),
  list("Base case", "Start age, sex, race/ethnicity", "set on the Cohort sheet", "", NA, NA))
r <- nrow(ST) + 3
wb <- put(wb, sh, r, 1, "Derived and fixed (not imported)"); wb <- font(wb, sh, cell(r, 1), bold = TRUE)
wb <- fill(wb, sh, wb_dims(rows = r, cols = 1:6), "E3E8EF")
for (d in derived) {
  r <- r + 1
  wb <- put(wb, sh, r, 1, d[[1]]); wb <- put(wb, sh, r, 2, d[[2]])
  if (is.character(d[[3]]) && startsWith(d[[3]], "=")) wb <- wb_add_formula(wb, sh, x = sub("^=", "", d[[3]]), dims = cell(r, 4))
  else wb <- put(wb, sh, r, 4, d[[3]])
  wb <- put(wb, sh, r, 5, d[[4]]); wb <- put(wb, sh, r, 6, d[[5]])
  wb <- font(wb, sh, wb_dims(rows = r, cols = 1:6)); wb <- align(wb, sh, wb_dims(rows = r, cols = 1:6))
  for (c in c(2, 5)) wb <- align(wb, sh, cell(r, c), wrap = TRUE)
  wb <- align(wb, sh, cell(r, 4), wrap = TRUE, horizontal = "left")
  if (!is.na(d[[6]])) wb <- numfmt(wb, sh, cell(r, 4), d[[6]])
}
wb <- wb_set_col_widths(wb, sh, cols = 1:6, widths = c(18, 44, 16, 40, 48, 10))
wb <- wb_freeze_pane(wb, sh, first_active_row = 2, first_active_col = 1)
wb <- note(wb, sh, r + 2, "Only rows with a Key are imported (blue Value cells). Lists are comma-separated.")

# ---- Cohort --------------------------------------------------------------------
C <- rd("inputs/cohort.csv")
Cn <- as.data.frame(lapply(C, function(col) { v <- lapply(col, num); if (all(vapply(v, function(z) is.na(z) || is.numeric(z), TRUE))) unlist(v) else col }),
                    check.names = FALSE, stringsAsFactors = FALSE)
wb <- sheet(wb, "Cohort", Cn, c(12, 10, 8, 10, 16, 10, 16, 13, 80), fmts = list(`6` = "0.0%"), value_cols = 2:4, wrap_cols = 9)
for (i in seq_len(nrow(C))) {
  kr <- row_of[paste0("kunz_share_", C$sex[i])]
  link <- if (!is.na(kr)) switch(C$risk_group[i], kunzmann_high = sprintf("=Settings!D%d", kr), kunzmann_low = sprintf("=1-Settings!D%d", kr), NULL)
  if (!is.null(link)) wb <- wb_add_formula(wb, "Cohort", x = sub("^=", "", link), dims = cell(i + 1, 6))
  else wb <- font(wb, "Cohort", cell(i + 1, 6), "0000FF")
}
wb <- note(wb, "Cohort", nrow(C) + 3, "Shares of the Kunzmann groups are linked to Settings (kunz_share_male / kunz_share_female) and set from it on import. Ref [8].")

# ---- Scenarios -----------------------------------------------------------------
SC <- rd("inputs/scenarios.csv")
if (!"refs" %in% names(SC)) SC$refs <- ""
lab <- vapply(names(P), function(k) P[[k]]$label, "")
basev <- vapply(SC$parameter, function(k) as.numeric(P[[k]]$value), 0)
scn <- data.frame(Scenario = SC$scenario, Parameter = SC$parameter, Label = unname(lab[SC$parameter]),
                  `Base value` = unname(basev), `Scenario value` = as.numeric(SC$value),
                  Refs = ifelse(nzchar(SC$refs), SC$refs, NONE_REF), check.names = FALSE, stringsAsFactors = FALSE)
wb <- sheet(wb, "Scenarios", scn, c(58, 16, 32, 13, 14, 16), value_cols = c(1, 2, 5, 6))

# ---- Workup counts --------------------------------------------------------------
CT <- rd("inputs/counts.csv")
ctd <- data.frame(Branch = CT$branch, Outcome = CT$outcome, `Sets parameter` = ifelse(nzchar(CT$parameter), CT$parameter, NA),
                  `n (patients)` = suppressWarnings(as.numeric(CT$n)), `What to count` = CT$notes, check.names = FALSE, stringsAsFactors = FALSE)
wb <- sheet(wb, "Workup counts", ctd, c(12, 20, 16, 13, 70), value_cols = 4, wrap_cols = 5, table = FALSE)
wb <- fill(wb, "Workup counts", wb_dims(rows = 2:(nrow(ctd) + 1), cols = 4), "FFFF00")
k <- nrow(ctd) + 3
for (j in 0:2) wb <- note(wb, "Workup counts", k + j, c(
  "Fill the yellow n column with numbers of patients. Fill every row of a branch, or leave the whole branch blank (a partly filled branch stops the run).",
  "Example: egd positive = 42, egd negative = 58 gives P(EGD +ve) = (42 + 0.5) / (100 + 1) = 0.42 with a beta(42.5, 58.5) PSA distribution.",
  "Then run Rscript import_inputs_xlsx.R (writes inputs/counts.csv).")[j + 1], size = 10, italic = FALSE)

# ---- References (placed after Workup counts) --------------------------------------
REFS <- rd("inputs/references.csv")
used <- list()
for (k in names(P)) for (x in trimws(strsplit(g0(P[[k]], "refs"), ",")[[1]])) if (nzchar(x)) used[[x]] <- c(used[[x]], k)
extra <- list(`8` = "Kunzmann split (Settings, Cohort)", `9` = "Kunzmann RR on EAC share (Settings)",
  `17` = "utility_well (EQ-5D norms)", `19` = "daly_ref_le (WHO standard YLL)", `21` = "scenario: p_incid_cancer",
  `26` = c("inc_eso, eac_share, late_eso (Baseline); SEER incidence/stage sheets", "mort_ugi (SEER mortality sheets)"),
  `27` = c("inc_stom, late_stom (Baseline); SEER incidence/stage sheets", "mort_ugi (SEER mortality sheets)"),
  `28` = "inc_other (Baseline); SEER incidence/race-ratio sheets", `30` = "mort_other (Baseline); NCHS life tables sheet",
  `38` = "false-positive pathway disutilities (dfp_*)", `39` = "scenario: fp_tto = 0", `40` = "scenario: fp_tto = 0",
  `41` = "p_worry", `42` = "dur_worry", `36` = "cost inflation (CPI medical care sheet; Settings cost_year)",
  `37` = "repeat-round positivity (context)")
for (x in names(extra)) used[[x]] <- c(used[[x]], extra[[x]])
sh <- "References"; wb <- wb_add_worksheet(wb, sh)
wb <- wb_add_data(wb, sh, as.data.frame(t(c("#", "Reference", "DOI / URL", "PMID", "Used for (derived, not imported)"))), col_names = FALSE)
wb <- header(wb, sh, 5)
for (j in seq_len(nrow(REFS))) {
  r <- j + 1; x <- REFS[j, ]
  vals <- list(as.integer(x$ref), x$citation, if (nzchar(x$url)) x$url else NA, if (nzchar(x$pmid)) as.integer(x$pmid) else NA,
               paste(used[[x$ref]], collapse = ", "))
  for (c in 1:5) wb <- put(wb, sh, r, c, vals[[c]])
  wb <- font(wb, sh, wb_dims(rows = r, cols = 1:4), "0000FF"); wb <- font(wb, sh, cell(r, 5))
  wb <- align(wb, sh, wb_dims(rows = r, cols = 1:5)); for (c in c(2, 3, 5)) wb <- align(wb, sh, cell(r, c), wrap = TRUE)
  if (nzchar(x$url)) { wb <- wb_add_hyperlink(wb, sh, cell(r, 3), target = x$url); wb <- font(wb, sh, cell(r, 3), "0563C1", underline = TRUE) }
  if (nzchar(x$pmid)) { wb <- wb_add_hyperlink(wb, sh, cell(r, 4), target = sprintf("https://pubmed.ncbi.nlm.nih.gov/%s/", x$pmid))
                        wb <- font(wb, sh, cell(r, 4), "0563C1", underline = TRUE) }
}
wb <- wb_set_col_widths(wb, sh, cols = 1:5, widths = c(5, 80, 45, 11, 45))
wb <- wb_freeze_pane(wb, sh, first_active_row = 2, first_active_col = 2)
r0 <- nrow(REFS) + 3
wb <- put(wb, sh, r0, 2, "Inputs marked 'none (assumption)' have no published source; their reasoning is in the Source / derivation column of the numbered tree-step tabs. Citations and DOIs verified in PubMed (Nov 2026) except where noted.")
wb <- font(wb, sh, cell(r0, 2), size = 9, italic = TRUE); wb <- align(wb, sh, cell(r0, 2), wrap = TRUE)

# ---- Baseline (calibration): derived, not imported ---------------------------------
B <- rd("inputs/baseline_risks.csv")
RL <- c(all = "All races/ethnicities", nhw = "Non-Hispanic White", nhb = "Non-Hispanic Black",
        nhapi = "Non-Hispanic Asian/Pacific Islander", nhaian = "Non-Hispanic American Indian/Alaska Native", hisp = "Hispanic")
nb <- function(k) as.numeric(B[[k]])
bd <- data.frame(B$race_eth, unname(RL[B$race_eth]), B$sex, as.integer(B$age_from), nb("mort_other"), nb("mort_all"), nb("mort_ugi"),
                 nb("inc_eso") * 1e5, nb("inc_stom") * 1e5, NA_real_, nb("inc_other") * 1e5, nb("eac_share"), nb("late_eso"),
                 nb("late_stom"), nb("utility_well"), nb("daly_ref_le"), stringsAsFactors = FALSE)
names(bd) <- c("race_eth", "Race/ethnicity", "Sex", "Age", "Other-cause death prob. (cause-deleted) [30, 26, 27]",
  "All-cause qx (NCHS 2024) [30]", "Upper GI cancer death prob. (SEER US mortality) [26, 27]",
  "SEER esophagus incidence per 100,000 [26]", "SEER stomach incidence per 100,000 [27]",
  "Upper GI target per 100,000 (formula)", "Other cancers incidence per 100,000 [26-28]",
  "Adenocarcinoma share of esophageal [26]", "Late-stage share at dx, esophagus [26]",
  "Late-stage share at dx, stomach [27]", "Utility, general population [17] (EQ-5D norms sheet)",
  "Standard YLL per death [19] (WHO standard YLL sheet)")
sh <- "Baseline (calibration)"
wb <- sheet(wb, sh, bd, c(9, 30, 8, 6, 14, 12, 14, 14, 14, 14, 14, 13, 13, 13, 12, 12),
            fmts = list(`5` = "0.000000", `6` = "0.000000", `7` = "0.000000", `8` = "0.00", `9` = "0.00", `10` = "0.00", `11` = "0.0",
                        `12` = "0.0%", `13` = "0.0%", `14` = "0.0%", `15` = "0.000", `16` = "0.00"), freeze = c(2, 5))
nrB <- nrow(bd)
wb <- wb_add_formula(wb, sh, x = sprintf("H%d+I%d", 2:(nrB + 1), 2:(nrB + 1)), dims = wb_dims(rows = 2:(nrB + 1), cols = 10))
wb <- wb_set_row_heights(wb, sh, 1, 66)
wb <- note(wb, sh, nrB + 3, "Numbers in [ ] refer to the References sheet. Race/ethnicity groups without SEER age-specific rates are scaled by age-adjusted ratios [26-28]; stage for groups other than White/Black uses all races.")
wb <- note(wb, sh, nrB + 4, "DERIVED SHEET - not imported. Rebuilt from the raw-data sheets by data-raw/make_us_inputs.R on import; edit those sheets instead.",
           size = 10, italic = FALSE, bold = TRUE, color = "C00000")

# ---- Raw source tables ---------------------------------------------------------------
raw <- function(f) rd(file.path("data-raw/raw", f))
band <- function(a) ifelse(a == "85", "85+", paste0(a, "-", as.integer(a) + 4))
src <- list()
R1 <- raw("seer_csr_2013_2017.csv")
wb <- sheet(wb, "SEER incidence", data.frame(Site = R1$site, Race = R1$race, Sex = R1$sex, `Age band (start)` = band(R1$age_band),
            `Rate per 100,000` = as.numeric(R1$rate), check.names = FALSE), c(12, 10, 10, 16, 16), fmts = list(`5` = "#,##0.0"), value_cols = 5)
src[["SEER incidence"]] <- c(nrow(R1), "[26, 27, 28] SEER Cancer Statistics Review 1975-2017, Tables 2.7 (all sites), 8.7 (esophagus), 24.7 (stomach): age-specific SEER 21 incidence 2013-2017.")
R2 <- raw("seer_csr_age_adjusted_by_race.csv")
wb <- sheet(wb, "SEER race ratios", data.frame(Site = R2$site, `Race/ethnicity` = R2$race, `Male, age-adjusted per 100,000` = as.numeric(R2$male),
            `Female, age-adjusted per 100,000` = as.numeric(R2$female), check.names = FALSE), c(12, 14, 18, 18),
            fmts = list(`3` = "#,##0.0", `4` = "#,##0.0"), value_cols = 3:4)
src[["SEER race ratios"]] <- c(nrow(R2), "[26, 27, 28] SEER CSR 1975-2017 Tables 2.15, 8.15, 24.15: age-adjusted SEER 21 incidence 2013-2017 by race/ethnicity and sex.")
R3 <- raw("seer_stage_by_race_sex.csv")
wb <- sheet(wb, "SEER stage", data.frame(Site = R3$site, Race = R3$race, Sex = R3$sex, `Localized %` = as.integer(R3$localized),
            `Regional %` = as.integer(R3$regional), `Distant %` = as.integer(R3$distant), `Unstaged %` = as.integer(R3$unstaged),
            `Late share among staged (formula)` = NA_real_, check.names = FALSE), c(12, 10, 10, 12, 12, 12, 12, 16), value_cols = 4:7)
i3 <- 2:(nrow(R3) + 1)
wb <- wb_add_formula(wb, "SEER stage", x = sprintf("(E%d+F%d)/(D%d+E%d+F%d)", i3, i3, i3, i3, i3), dims = wb_dims(rows = i3, cols = 8))
wb <- numfmt(wb, "SEER stage", wb_dims(rows = i3, cols = 8), "0.0%")
src[["SEER stage"]] <- c(nrow(R3), "[26, 27] SEER CSR 1975-2017 Tables 8.8, 24.8: stage distribution 2010-2016 (SEER Summary Stage 2000).")
R4 <- raw("seer_eac_share_by_race.csv")
wb <- sheet(wb, "EAC share", data.frame(race_eth = R4$race, `Race/ethnicity` = unname(RL[R4$race]),
            `Adenocarcinoma share of esophageal cancer` = as.numeric(R4$eac_share), check.names = FALSE), c(10, 40, 22),
            fmts = list(`3` = "0.0%"), value_cols = 3)
src[["EAC share"]] <- c(nrow(R4), "[26] SEER CSR 1975-2017 Table 8.22: histology by race, 2013-2017 (White used for NH White, Black for NH Black).")
R5 <- raw("nchs_lifetables_2024_qx.csv")
grp <- unique(R5[, c("race_eth", "sex")])
nd <- data.frame(Age = 50:99)
for (j in seq_len(nrow(grp))) {
  s <- R5[R5$race_eth == grp$race_eth[j] & R5$sex == grp$sex[j], ]
  nd[[paste(RL[[grp$race_eth[j]]], grp$sex[j])]] <- as.numeric(s$qx[match(50:99, as.integer(s$age))])
}
fm <- setNames(as.list(rep("0.000000", nrow(grp))), as.character(seq_len(nrow(grp)) + 1))
wb <- sheet(wb, "NCHS life tables", nd, c(6, rep(16, nrow(grp))), fmts = fm, value_cols = seq_len(nrow(grp)) + 1, freeze = c(2, 2), table = FALSE)
wb <- wb_set_row_heights(wb, "NCHS life tables", 1, 54)
src[["NCHS life tables"]] <- c(nrow(nd), "[30] Arias E, Xu JQ, Tejada-Vera B. United States Life Tables, 2024. NVSR vol 75 no 5 (Tables 2-3, 5-6, 8-9, 11-12, 14-15, 17-18).")
R8 <- raw("seer_us_mortality_2013_2017.csv")
wb <- sheet(wb, "SEER mortality", data.frame(Site = R8$site, Race = R8$race, Sex = R8$sex, `Age band (start)` = band(R8$age_band),
            `US death rate per 100,000` = as.numeric(R8$rate), check.names = FALSE), c(12, 10, 10, 16, 18), fmts = list(`5` = "#,##0.0"), value_cols = 5)
src[["SEER mortality"]] <- c(nrow(R8), "[26, 27] SEER CSR 1975-2017 Tables 8.7, 24.7: U.S. age-specific death rates 2013-2017. Used to remove upper GI cancer deaths from the life tables (cause-deleted other-cause mortality).")
R9 <- raw("seer_us_mortality_age_adjusted_by_race.csv")
wb <- sheet(wb, "SEER mortality ratios", data.frame(Site = R9$site, `Race/ethnicity` = R9$race,
            `Male, age-adjusted death rate per 100,000` = as.numeric(R9$male), `Female, age-adjusted death rate per 100,000` = as.numeric(R9$female),
            check.names = FALSE), c(12, 14, 20, 20), fmts = list(`3` = "#,##0.0", `4` = "#,##0.0"), value_cols = 3:4)
src[["SEER mortality ratios"]] <- c(nrow(R9), "[26, 27] SEER CSR 1975-2017 Tables 8.15, 24.15: U.S. age-adjusted death rates 2013-2017 by race/ethnicity; scale the all-races/White curves like incidence.")
R10 <- raw("cpi_medical_care.csv")
wb <- sheet(wb, "CPI medical care", data.frame(Year = as.integer(R10$year), `Annual average (index)` = as.numeric(R10$index),
            `First half` = as.numeric(R10$h1), `Second half` = as.numeric(R10$h2), check.names = FALSE), c(8, 20, 12, 12),
            fmts = list(`2` = "0.000", `3` = "0.000", `4` = "0.000"), value_cols = 2)
src[["CPI medical care"]] <- c(nrow(R10), "[36] BLS CPI-U Medical Care, U.S. city average (CUUS0000SAM, semiannual; annual = mean of halves), via FRED. Annual average is the column used.")
R6 <- raw("eq5d_us_norms.csv")
wb <- sheet(wb, "EQ-5D norms", data.frame(`Age from` = as.integer(R6$age_from), `Age to` = as.integer(R6$age_to),
            `EQ-5D-5L utility` = as.numeric(R6$utility), check.names = FALSE), c(10, 10, 16), fmts = list(`3` = "0.000"), value_cols = 1:3)
src[["EQ-5D norms"]] <- c(nrow(R6), "[17] Jiang R, Janssen MFB, Pickard AS. Qual Life Res 2021;30:803-16 (US EQ-5D-5L norms by age band). Bands must cover ages 50-99 once each.")
R7 <- raw("who_ghe_standard_yll.csv")
wb <- sheet(wb, "WHO standard YLL", data.frame(`Age (band start; 85 = 85+)` = as.integer(R7$age), `Standard YLL per death` = as.numeric(R7$yll),
            check.names = FALSE), c(22, 20), fmts = list(`2` = "0.00"), value_cols = 1:2)
src[["WHO standard YLL"]] <- c(nrow(R7), "[19] WHO GHE 2019 DALY methods, Table 2.1 standard life table: remaining life expectancy at the start of each 5-year band.")
for (nm in names(src)) wb <- note(wb, nm, as.integer(src[[nm]][1]) + 3, paste("Source:", src[[nm]][2]))

# README pointer to the References sheet
wb <- put(wb, "README", length(lines) + 2, 1, "References - numbered citations with DOI/PubMed links; the Refs columns on other sheets point to these numbers.")
wb <- font(wb, "README", cell(length(lines) + 2, 1))

# sheet order: README, tree-step tabs, Settings, Cohort, Scenarios, Workup counts, References, Baseline, raw tables
ord <- c("README", vapply(STEPS, `[[`, "", 1), "Settings", "Cohort", "Scenarios", "Workup counts", "References",
         "Baseline (calibration)", names(src))
wb <- wb_set_order(wb, match(ord, wb$get_sheet_names()))
wb_save(wb, OUT, overwrite = TRUE)
cat("saved", OUT, "\n")
