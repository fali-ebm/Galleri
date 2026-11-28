# =============================================================
# The two harmonized models (run after run.R and tools/run_microsim.R):
#
#   A. NET BURDEN MODEL - what annual Galleri screening does to a population, in natural units,
#      undiscounted: benefits (extra stage I-II diagnoses, late-stage diagnoses and cancer deaths
#      averted, life-years) against harms (people without cancer worked up, invasive procedures,
#      complications, overdiagnosis, radiation, patient time), net QALYs and DALYs, and the
#      endoscopy capacity a national programme would need.
#   B. COMPARATIVE MODEL - Galleri (and the Kunzmann -> EGD strategy) against usual care
#      ("do nothing"): discounted costs and QALYs, ICER, net monetary benefit, the efficiency
#      frontier, and the primary endpoint, cost per extra stage I-II diagnosis (stage I-II =
#      endoscopically resectable or stage II). Headline = microsimulation (1M people) when
#      available; PSA, VOI and subgroups come from the cohort model (run.R).
#
# Both use one natural history (SEER-calibrated, MISCAN-ESO for EAC, two-stage preclinical
# squamous and gastric cancer), the same Galleri accuracy (PMA test, stage-specific) and the
# same workup tree; they differ only in perspective and discounting.
#
#   Rscript tools/run_harmonized.R     -> outputs/harmonized/{burden.csv, comparative.csv,
#                                         frontier.csv, harmonized.json}
# =============================================================
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
proj <- ceakit::load_project("."); E <- proj$env
p <- ceakit::param_list(proj$params); S <- proj$settings; D <- proj$data
out <- "outputs/harmonized"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
PER <- 1e5

# Populations: the base case (cohort.csv: 50-year-old non-Hispanic White men) and the US
# population at 50 (all races; men and women equally weighted), used for national figures.
us_cohort <- function(s) {
  m <- E$subgroup_cohort(50, "male", "all", s); f <- E$subgroup_cohort(50, "female", "all", s)
  m$share <- m$share / 2; f$share <- f$share / 2; rbind(m, f)
}
run_pop <- function(arms, s, cohort = NULL) {
  cohort <- cohort %||% E$cohort_rows(D)
  A <- E$run_population(p, D, s, arms, cohort)$total
  rbind(A, E$derive_galleri(A, s))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- A. Net burden (undiscounted) -------------------------------------------------
S0 <- S; S0$disc_cost <- 0; S0$disc_qaly <- 0; S0$disc_daly <- 0
patient_hours <- function(m) m["galleri_tests", ] * p$pt_galleri + m["egd", ] * p$pt_egd + m["eus", ] * p$pt_eus +
  m["pet_ct", ] * p$pt_pet + m["ent", ] * p$pt_ent + m["laparoscopy", ] * p$pt_lap +
  (m["incidental_followups", ] + m["surveillance_imaging", ] + m["tumor_board", ]) * p$pt_other
burden <- function(m, pop) {
  g <- "Galleri"; u <- "Usual care"; d <- function(k) m[k, g] - m[k, u]
  hrs <- patient_hours(m); hrs_free <- m["galleri_tests", g] * 0 +
    m["egd_free", g] * p$pt_egd + m["eus_free", g] * p$pt_eus + m["pet_ct_free", g] * p$pt_pet + m["ent_free", g] * p$pt_ent +
    m["laparoscopy_free", g] * p$pt_lap + (m["incidental_followups_free", g] + m["surveillance_free", g]) * p$pt_other
  # national, steady state: a stationary 50+ population screened since 50 has, each year, the
  # cohort's lifetime events per life-year lived x the US 50+ population (us_pop_50plus)
  yr <- function(k) p$us_pop_50plus * m[k, g] / m["ly_total", g]
  rows <- list(
    c("Benefit", "Extra stage I-II upper GI diagnoses (endoscopically resectable or stage II)", d("stage12_ugi")),
    c("Benefit", "  - esophageal (primary focus)", d("stage12_eso")),
    c("Benefit", "    - adenocarcinoma", d("stage12_eac")), c("Benefit", "    - squamous cell", d("stage12_escc")),
    c("Benefit", "  - gastric (concomitant)", d("stage12_gastric")),
    c("Benefit", "Stage III-IV upper GI diagnoses averted", -d("stage34_ugi")),
    c("Benefit", "  - esophageal", -d("stage34_eso")),
    c("Benefit", "Upper GI cancer deaths averted", -d("ca_deaths")),
    c("Benefit", "Life-years gained (undiscounted)", d("ly_total")),
    c("Benefit", "Other-site cancers found by the workup", m["found_other", g]),
    c("Harm", "False-positive workups (people without cancer; repeats counted)", m["free_investigated", g]),
    c("Harm", "  - repeat false positives (same person)", m["repeat_false_positives", g]),
    c("Harm", "Invasive procedures in people without cancer", m["invasive_free", g]),
    c("Harm", "EGDs (all)", m["egd", g]), c("Harm", "EUS (all)", m["eus", g]), c("Harm", "PET CT (all)", m["pet_ct", g]),
    c("Harm", "Diagnostic laparoscopies", m["laparoscopy", g]),
    c("Harm", "Serious complications", m["sae_total", g]), c("Harm", "Deaths from procedure complications", m["proc_deaths", g]),
    c("Harm", "Overdiagnosed cancers", m["overdiagnosed", g]),
    c("Harm", "Radiation-induced cancers", m["radiation_cancers", g]),
    c("Harm", "Patient time, all screening and workup (hours)", hrs[[g]] - hrs[[u]]),
    c("Harm", "  - of which workup of people without cancer (hours)", hrs_free),
    c("Net", "Net QALYs (undiscounted)", d("qaly_total")),
    c("Net", "Net DALYs averted (undiscounted)", -d("daly_total")),
    c("Net", "QALYs lost by people without cancer worked up", m["qaly_fp", g]),
    c("Ratio", "People screened per extra stage I-II diagnosis", 1 / d("stage12_ugi")),
    c("Ratio", "People screened per cancer death averted", 1 / -d("ca_deaths")),
    c("Ratio", "False-positive workups per extra stage I-II diagnosis", m["free_investigated", g] / d("stage12_ugi")),
    c("Ratio", "False-positive workups per death averted", m["free_investigated", g] / -d("ca_deaths")),
    c("Ratio", "EGDs per extra stage I-II diagnosis", m["egd", g] / d("stage12_ugi")),
    c("Ratio", "Patient hours per extra stage I-II diagnosis", (hrs[[g]] - hrs[[u]]) / d("stage12_ugi")),
    c("Capacity (US, per year)", "Galleri tests", yr("galleri_tests")),
    c("Capacity (US, per year)", "EGDs", yr("egd")),
    c("Capacity (US, per year)", "  - share of current US upper endoscopy volume (7.46M)", yr("egd") / p$us_upper_endo),
    c("Capacity (US, per year)", "EUS", yr("eus")),
    c("Capacity (US, per year)", "  - share of current US upper EUS volume (290,655)", yr("eus") / p$us_eus),
    c("Capacity (US, per year)", "PET CT", yr("pet_ct")),
    c("Capacity (US, per year)", "Extra stage I-II upper GI diagnoses", p$us_pop_50plus * d("stage12_ugi") / m["ly_total", g]),
    c("Capacity (US, per year)", "Patient hours", p$us_pop_50plus * (hrs[[g]] - hrs[[u]]) / m["ly_total", g]))
  x <- rbindlist(lapply(rows, function(r) data.table(section = r[1], measure = r[2], per_person = as.numeric(r[3]))))
  x[, population := pop]
  x[, per_100k := ifelse(section %in% c("Ratio", "Capacity (US, per year)"), NA_real_, per_person * PER)]
  x[, value := ifelse(section %in% c("Ratio", "Capacity (US, per year)"), per_person, per_100k)]
  x[, .(population, section, measure, value)]
}
mb <- run_pop(c("Usual care", "Galleri"), S0)
mu <- run_pop(c("Usual care", "Galleri"), S0, us_cohort(S0))
B <- rbind(burden(mb, "Base case: 50-year-old non-Hispanic White men"), burden(mu, "US 50-year-olds, all races, both sexes"))
fwrite(B, file.path(out, "burden.csv"))

# ---- B. Comparative model vs usual care (discounted) --------------------------------
arms <- E$ARMS
mc <- run_pop(arms, S)
cmp <- rbindlist(lapply(arms, function(a) data.table(
  arm = a, cost = mc["cost_total", a], qaly = mc["qaly_total", a], ly = mc["ly_total", a],
  stage12_ugi = PER * mc["stage12_ugi", a], stage12_eso = PER * mc["stage12_eso", a], stage34_ugi = PER * mc["stage34_ugi", a],
  stage12_gained_ugi = PER * mc["stage12_gained_ugi", a], stage12_gained_eso = PER * mc["stage12_gained_eso", a],
  stage12_share_eso = mc["stage12_share_eso", a],
  ca_deaths_averted = PER * mc["ca_deaths_averted", a], free_investigated = PER * mc["free_investigated", a],
  inc_cost = mc["cost_total", a] - mc["cost_total", "Usual care"], inc_qaly = mc["qaly_total", a] - mc["qaly_total", "Usual care"],
  cost_per_stage12_gained = mc["cost_per_stage12_gained", a], cost_per_stage12_gained_eso = mc["cost_per_stage12_gained_eso", a])))
cmp[, icer_vs_usual_care := ifelse(arm == "Usual care", NA, inc_cost / inc_qaly)]
cmp[, nmb_100k := qaly * S$wtp - cost]
cmp[, nmb_vs_usual_care := nmb_100k - nmb_100k[arm == "Usual care"]]
# efficiency frontier: drop strongly dominated, then extendedly dominated strategies
fr <- cmp[order(cost)][, .(arm, cost, qaly)]; fr[, status := "frontier"]
repeat {
  f <- fr[status == "frontier"]; ch <- FALSE
  for (i in seq_len(nrow(f))[-1]) if (f$qaly[i] <= max(f$qaly[seq_len(i - 1)])) { fr[arm == f$arm[i], status := "dominated"]; ch <- TRUE; break }
  if (ch) next
  f <- fr[status == "frontier"]
  if (nrow(f) > 2) { ic <- diff(f$cost) / diff(f$qaly)
    j <- which(diff(ic) < 0)
    if (length(j)) { fr[arm == f$arm[j[1] + 1], status := "extendedly dominated"]; next } }
  break
}
f <- fr[status == "frontier"]; fr[status == "frontier", icer_on_frontier := c(NA, diff(f$cost) / diff(f$qaly))]
fwrite(cmp, file.path(out, "comparative.csv")); fwrite(fr, file.path(out, "frontier.csv"))

# headline from the microsimulation (same inputs) when it has been run
ms <- if (file.exists("outputs/microsim/summary.csv")) fread("outputs/microsim/summary.csv") else NULL
ceac <- if (file.exists("outputs/comparison/ceac.csv")) fread("outputs/comparison/ceac.csv") else NULL
J <- list(meta = list(wtp = S$wtp, disc = c(cost = S$disc_cost, qaly = S$disc_qaly), us_pop_50plus = p$us_pop_50plus,
                      galleri_stage_sens = lapply(E$galleri_stage_sens(p, D), function(v) if (length(v) > 1) as.list(v) else v)),
          burden = B, comparative = cmp, frontier = fr, microsim = ms, ceac = ceac)
write_json(J, file.path(out, "harmonized.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("\n== Net burden (undiscounted), base case, per 100,000 ==\n"); print(B[population == B$population[1]], nrows = 80)
cat("\n== Comparative (discounted 3%), per person ==\n"); print(cmp[, .(arm, cost, qaly, inc_cost, inc_qaly, icer_vs_usual_care, stage12_gained_eso, cost_per_stage12_gained)])
print(fr)
if (!is.null(ms)) { cat("\n== Microsimulation headline ==\n"); print(ms) }
