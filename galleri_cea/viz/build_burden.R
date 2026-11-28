# Build viz/galleri_burden.html from viz/galleri_viz_data.json.
#   Rscript viz/build_burden.R        (run from galleri_cea/, after Rscript tools/viz_data.R)
# The page template (viz/burden_template.html) holds the D3 code; this script computes
# the numbers it shows and inlines them as one JSON object in place of __DATA__.
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))   # UTF-8 text (dashes, arrows) whatever the system locale
suppressPackageStartupMessages(library(jsonlite))
here <- "viz"
D <- fromJSON(file.path(here, "galleri_viz_data.json"), simplifyVector = FALSE)
G <- Filter(function(r) r$age <= D$meta$stop_age, D$by_age$Galleri)
R <- list(); for (r in D$galleri) if (r$arm == "Galleri") R[[r$key]] <- r
P <- list(); for (r in D$params) P[[r$parameter]] <- r$value
col <- function(rows, f) vapply(rows, f, numeric(1))

# Whole people per age column for the dot field: largest-remainder rounding keeps the total.
largest_remainder <- function(v, total) {
  fl <- floor(v); rem <- total - sum(fl)
  o <- order(v - fl, decreasing = TRUE)
  if (rem > 0) fl[o[seq_len(rem)]] <- fl[o[seq_len(rem)]] + 1
  fl
}
groups <- list(
  found      = col(G, function(r) r$found_ugi + r$found_other - r$overdx_ugi - r$overdx_other),
  overdx     = col(G, function(r) r$overdx_ugi + r$overdx_other),
  missed     = col(G, function(r) r$tp + r$fp_other_ca - r$found_ugi - r$found_other),
  falsealarm = col(G, function(r) r$free_inv))
tot <- lapply(groups, function(v) round(sum(v)))
alloc <- Map(largest_remainder, groups, tot)
ages <- lapply(seq_along(G), function(i) {
  r <- G[[i]]
  # the whole-people counts (alloc) replace the fractional overdx, as the page expects
  modifyList(list(age = r$age, tests = r$tests, alive = r$alive, tp = r$tp, fn = r$fn, fn_symp = r$fn_symp,
                  found_ugi = r$found_ugi, found_other = r$found_other, overdx = r$overdx_ugi + r$overdx_other,
                  pos = r$tp + r$fp), lapply(alloc, `[`, i))
})
val <- function(k) R[[k]]$value
n_fp <- val("free_investigated")
ledger <- list(
  n_fp = n_fp, ent = val("ent_free") / n_fp, pet = val("pet_ct_free") / n_fp, eus = val("eus_free") / n_fp,
  samp = val("eus_sampling_free") / n_fp, incid = val("incidental_followups_free") / n_fp,
  days = val("mean_days_fp"), cost = val("cost_fp") * 1e5 / n_fp, msv = val("dose_sv_free") * 1000 / n_fp,
  qaly_days = val("qaly_fp") / n_fp * 365.25, perf_1_in = 1 / P$p_perf_egd,
  cost_total = val("cost_fp") * 1e5, qaly_total = val("qaly_fp"))
fmt <- function(x, d) formatC(x, format = "f", digits = d, big.mark = ",")
ci <- function(k, d = 0) { r <- R[[k]]; list(v = fmt(r$value, d), lo = fmt(r$lower_95, d), hi = fmt(r$upper_95, d)) }
table <- list(
  list("Galleri tests over 35 years", ci("galleri_tests")),
  list("People with a positive result, worked up", ci("investigated")),
  list("Cancers found (upper GI and other sites)", ci("screen_detected")),
  list("  of which overdiagnosed", ci("overdiagnosed")),
  list("Cancer-free people worked up (false alarms)", ci("free_investigated")),
  list("PET CT scans in cancer-free people", ci("pet_ct_free")),
  list("EGD perforations, all people worked up", ci("perforations", 2)),
  list("Radiation-induced cancers, all people worked up", ci("radiation_cancers", 2)),
  list("False-negative test results", ci("screen_fn")),
  list("  presenting with symptoms within a year", ci("fn_symptomatic_1y")),
  list("QALYs lost to the workup and its harms, everyone worked up", ci("qaly_harms", 1)),
  list("QALYs lost by cancer-free people worked up", ci("qaly_fp", 1)),
  list("DALYs incurred by cancer-free people worked up", ci("daly_fp", 1)))
pair <- function(lbl, k) list(lbl, val(k))
fpb <- list(
  counts = list(pair("Negative EGDs", "egd_free"), pair("Negative PET CT scans", "pet_ct_free"),
    pair("Negative ENT exams", "ent_free"), pair("Negative endoscopic ultrasounds", "eus_free"),
    pair("  of which with a biopsy", "eus_sampling_free"), pair("Negative diagnostic laparoscopies", "laparoscopy_free"),
    pair("Watchful-waiting imaging", "surveillance_free"), pair("Follow-ups of a harmless PET finding", "incidental_followups_free"),
    pair("Days waiting for the all-clear", "days_free"), pair("Still worried 6 months later", "fp_worry_cases"),
    pair("EGD perforations", "perforations_free"), pair("Other serious complications", "fp_sae_other"),
    pair("Deaths from complications", "proc_deaths_free"), pair("Radiation-induced cancers, later in life", "radiation_cancers_free")),
  comp = list(
    list("Workup: procedures and uncertainty", val("fp_qaly_pathway"), val("fp_daly_anxiety")),
    list("Persisting worry", val("fp_qaly_worry"), val("fp_daly_worry")),
    list("Radiation-induced cancers", val("fp_qaly_radiation"), val("fp_daly_radiation")),
    list("Deaths from complications", val("fp_qaly_deaths"), val("fp_daly_deaths")),
    list("Non-fatal complications", val("fp_qaly_complications"), val("fp_daly_complications"))),
  qaly = list(val("qaly_fp"), R$qaly_fp$lower_95, R$qaly_fp$upper_95),
  daly = list(val("daly_fp"), R$daly_fp$lower_95, R$daly_fp$upper_95),
  icer = val("icer"), qaly_net = val("qaly_net"))
SC <- list(); for (r in D$scenarios) if (r$arm == "Galleri") SC[[r$scenario]] <- r$icer
fpb$icer_eq5d <- SC[["False positives: per-procedure loss only (EQ-5D finds none)"]]
keep <- c("sex", "race_eth", "race_label", "age", "tests", "tp", "fp", "fp_other_cancer", "fp_cancer_free", "tn", "fn",
          "prevalence_per_100k_tested", "ppv_upper_gi", "ppv_any_cancer", "npv", "cancer_free_per_cancer_found")
acc <- lapply(D$accuracy, function(r) r[keep])
data <- list(fpb = fpb, acc = acc, ages = ages, tot = tot, links = D$workup_links, ledger = ledger, table = table,
  tests_total = val("galleri_tests"), investigated = val("investigated"), found = val("screen_detected"),
  overdx = val("overdiagnosed"), fn = val("screen_fn"), fn_symp = val("fn_symptomatic_1y"),
  price = P$c_galleri, cost_year = D$meta$cost_year,
  net = utils::read.csv("outputs/harmonized/burden.csv", stringsAsFactors = FALSE), us_pop = P$us_pop_50plus)

js <- as.character(toJSON(data, auto_unbox = TRUE, digits = NA, null = "null", na = "null"))
tpl <- paste(readLines(file.path(here, "burden_template.html"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
html <- sub("__DATA__", js, tpl, fixed = TRUE)
writeLines(html, file.path(here, "galleri_burden.html"), useBytes = TRUE)
cat("wrote viz/galleri_burden.html", nchar(html), "found", tot$found, "overdx", tot$overdx,
    "missed", tot$missed, "falsealarm", tot$falsealarm, "\n")
