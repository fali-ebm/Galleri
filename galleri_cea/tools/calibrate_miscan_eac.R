# Calibrate the esophageal adenocarcinoma (EAC) natural history to CISNET MISCAN-ESO and
# check it against published targets.
#
#   Rscript tools/calibrate_miscan_eac.R          (run from galleri_cea/)
#
# What it does, step by step (every number written to outputs/miscan_eac/):
#   1. Targets. MISCAN-ESO (Erasmus/UW) mean preclinical EAC sojourn: men 4.9 y, women 4.5 y
#      (Kong 2014, CEBP 23:997, doi:10.1158/1055-9965.EPI-13-1233; profile calibration target
#      4-5 y). SEER esophagus stage at diagnosis by race and sex (late = regional + distant
#      among staged). SEER age-specific esophageal incidence x adenocarcinoma share.
#   2. Rates. miscan_rates() turns (sojourn S, late share L, dwell ratio r) into the three
#      exponential rates of the two-stage preclinical process: PE -> symptoms early (hE),
#      PE -> PL (lam), PL -> symptoms late (hL). Closed form, so the fit is exact.
#   3. Onset. calibrate() solves EAC onset year by year so that symptomatic EAC per living
#      person equals the SEER target at every age (and the rest of upper GI cancer likewise).
#   4. Checks. Mean sojourn recovered from the yearly transition matrices; symptomatic EAC
#      incidence and late share by age vs SEER; lifetime EAC cases and deaths under usual care
#      vs the CISNET comparative-model averages for White men and women (Rubenstein 2022,
#      Gastroenterology 163:163, doi:10.1053/j.gastro.2022.03.037): 654 / 536 per 100,000 men,
#      125 / 103 per 100,000 women, birth cohort 1950 followed from age 40.
#   MISCAN-ESO's Barrett's esophagus layer (prevalence 1.4% at 60-64, progression 0.18%/y,
#   dwell times 6.7 / 1.0 / 1.1 y) is not part of this model: Galleri and a one-time EGD for
#   cancer are modelled, not Barrett's surveillance. Those targets are listed as out of scope.
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
suppressPackageStartupMessages({ library(ceakit) })
proj <- load_project("."); E <- proj$env; D <- proj$data; S <- proj$settings
p <- param_list(proj$params); p$miscan_eac <- 1
out <- "outputs/miscan_eac"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- 2. rates by sex (base case race: non-Hispanic White) --------------------------
rate_row <- function(sex, race = "nhw", start = 50) {
  row <- E$subgroup_cohort(start, sex, race, S)[2, ]                  # kunzmann_low row: rr does not change L
  curve <- E$baseline_for(D, row, S); m <- E$eac_rates(p, curve, row, S)
  data.frame(sex = sex, race_eth = race, sojourn_target_y = m$S, late_share_target = m$L_seen,
             late_share_without_competing_death = m$L, dwell_ratio = m$r,
             rate_PE_to_symptoms_early = m$hE, rate_PE_to_PL = m$lam, rate_PL_to_symptoms_late = m$hL,
             mean_years_in_PE = 1 / m$a, mean_years_in_PL = 1 / m$hL,
             share_found_late_if_screened_once = (m$lam / m$hL) / (1 + m$lam / m$hL))   # steady-state PL / (PE + PL)
}
rates <- rbind(rate_row("male"), rate_row("female"))
write.csv(rates, file.path(out, "rates.csv"), row.names = FALSE)

# Mean sojourn recovered from the yearly matrices (no other-cause death): expected years
# preclinical, counting each year a person starts preclinical (annual-cycle analogue of S).
sojourn_from_matrix <- function(m) {
  y <- E$miscan_year(m); x <- c(1, 0); tot <- 0
  for (k in 1:400) { tot <- tot + sum(x); x <- c(x[1] * y[["pe_pe"]], x[1] * y[["pe_pl"]] + x[2] * y[["pl_pl"]]) }
  tot - 0.5                                                           # presentation mid-year on average
}

# ---- 3-4. usual care, EAC only, for a cohort ------------------------------------------
eac_flows <- function(C) {
  wk <- E$workup_summary(p); C <- E$kunz_cohort(C, E$kunz_settings(p, S)); res <- NULL; err <- 0
  for (i in seq_len(nrow(C))) {
    row <- C[i, ]; prep <- E$prepare_row(row, D, p, S); curve <- prep$curve
    s <- E$simulate_arm("Usual care", row, curve, prep$cal, prep$p, S, wk, prep$V)
    tr <- s$trace; n <- nrow(curve)
    f <- t(vapply(seq_len(n), function(t) {
      P <- E$tm_at(curve, prep$cal, prep$p, t); x <- tr[t, ]
      # includes this year's onsets that present before the next yearly test (pre_fe / pre_fl, within_year_onset)
      clin_e <- x[["pre_ee"]] * (P["pre_ee", "early_e"] + P["pre_ee", "late_e"]) + x[["pre_el"]] * P["pre_el", "late_e"] +
                x[["pre_fe"]] * P["pre_fe", "early_e"] + x[["pre_fl"]] * P["pre_fl", "late_e"]
      late_e <- x[["pre_ee"]] * P["pre_ee", "late_e"] + x[["pre_el"]] * P["pre_el", "late_e"] + x[["pre_fl"]] * P["pre_fl", "late_e"]
      death_e <- x[["early_e"]] * P["early_e", "dead_cancer"] + x[["late_e"]] * P["late_e", "dead_cancer"]
      alive <- sum(x[!E$STATES %in% E$DEAD])
      c(clin = clin_e, late = late_e, death = death_e, alive = alive)
    }, numeric(4)))
    d <- data.frame(age = curve$age_from, f, seer = E$target_eac(curve, row, S), seer_late = curve$late_eso, share = row$share)
    err <- max(err, abs(d$clin / d$alive - d$seer) / d$seer)            # per cohort row: model vs SEER, every age
    res <- if (is.null(res)) d * 0 + d * row$share else res + d * row$share   # share-weighted (age, seer_late identical)
    res$age <- curve$age_from; res$seer_late <- curve$late_eso
  }
  attr(res, "max_rel_err") <- err
  res
}
by_sex <- lapply(c(male = "male", female = "female"), function(sx) {
  C <- E$subgroup_cohort(50, sx, "nhw", S); d <- eac_flows(C)
  d$model_incidence <- d$clin / d$alive; d$seer_incidence <- d$seer / 1   # seer already share-weighted per person-year
  d$model_late_share <- d$late / d$clin; d$sex <- sx; d
})
ba <- do.call(rbind, by_sex)
write.csv(ba[, c("sex", "age", "seer_incidence", "model_incidence", "seer_late", "model_late_share", "clin", "death", "alive")],
          file.path(out, "by_age.csv"), row.names = FALSE)

# ---- target table ----------------------------------------------------------------------
fit_row <- function(target, source, kind, value, model, note = "") data.frame(target = target, source = source, kind = kind,
  target_value = value, model_value = model, difference_pct = ifelse(is.na(value) | is.na(model), NA, 100 * (model - value) / value), note = note)
mtab <- list()
for (sx in c("male", "female")) {
  r <- rates[rates$sex == sx, ]; m <- E$miscan_rates(r$sojourn_target_y, r$late_share_without_competing_death, r$dwell_ratio)
  d <- by_sex[[sx]]; lab <- if (sx == "male") "men" else "women"
  mtab[[length(mtab) + 1]] <- fit_row(paste("Mean preclinical EAC sojourn,", lab, "(years)"),
    "MISCAN-ESO (Kong 2014): men 4.9, women 4.5; profile target 4-5", "calibration", r$sojourn_target_y,
    1 / m$a + m$L / m$hL, sprintf("annual-cycle equivalent %.2f y", sojourn_from_matrix(m)))
  mtab[[length(mtab) + 1]] <- fit_row(paste("Late share at symptomatic EAC diagnosis, NH White", lab),
    "SEER esophagus stage 2010-2016 (regional + distant among staged)", "calibration", r$late_share_target,
    sum(d$late) / sum(d$clin), "lifetime, incidence-weighted; constant by age in the model")
  mtab[[length(mtab) + 1]] <- fit_row(paste("Symptomatic EAC incidence by age 50-99, NH White", lab),
    "SEER 21 2013-2017 esophagus incidence x adenocarcinoma share", "calibration", NA,
    100 * attr(d, "max_rel_err"), "model_value = largest relative error, %, over ages and cohort rows")
  cis <- if (sx == "male") c(654, 536) else c(125, 103)
  mtab[[length(mtab) + 1]] <- fit_row(paste("Lifetime EAC cases per 100,000 White", lab, "(no screening)"),
    "CISNET 3-model average, cohort born 1950 from age 40 (Rubenstein 2022)", "external validation", cis[1], 1e5 * sum(d$clin),
    "model: NH White from age 50, SEER 2013-2017 rates; CISNET also counts ages 40-49")
  mtab[[length(mtab) + 1]] <- fit_row(paste("Lifetime EAC deaths per 100,000 White", lab, "(no screening)"),
    "CISNET 3-model average, cohort born 1950 from age 40 (Rubenstein 2022)", "external validation", cis[2], 1e5 * sum(d$death),
    "model deaths use SEER esophagus stage-specific survival (m_early_eac, m_late_eac)")
  mtab[[length(mtab) + 1]] <- fit_row(paste("EAC deaths per EAC case, White", lab),
    "CISNET 3-model average (Rubenstein 2022)", "external validation", cis[2] / cis[1], sum(d$death) / sum(d$clin))
}
mtab[[length(mtab) + 1]] <- fit_row("Barrett's esophagus prevalence at 60-64; progression to EAC 0.18%/y; BE dwell times",
  "MISCAN-ESO profile", "out of scope", NA, NA, "no Barrett's states: the model screens for cancer, not Barrett's surveillance")
tab <- do.call(rbind, mtab)
write.csv(tab, file.path(out, "targets.csv"), row.names = FALSE)

# ---- figure ---------------------------------------------------------------------------
png(file.path(out, "miscan_eac_calibration.png"), width = 2000, height = 1200, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
dm <- by_sex$male
plot(dm$age, 1e5 * dm$seer_incidence, type = "l", lwd = 6, col = "grey75", xlab = "Age", ylab = "Symptomatic EAC per 100,000 per year",
     main = "EAC incidence, NH White men")
lines(dm$age, 1e5 * dm$model_incidence, lwd = 2, col = "#1F3A5F")
legend("topleft", c("SEER target", "Model (MISCAN-ESO scenario)"), lwd = c(6, 2), col = c("grey75", "#1F3A5F"), bty = "n")
x <- tab[tab$kind == "external validation" & !grepl("per EAC case", tab$target), ]
bp <- barplot(rbind(x$target_value, x$model_value), beside = TRUE, col = c("grey75", "#1F3A5F"), border = NA,
              names.arg = c("Cases, men", "Deaths, men", "Cases, women", "Deaths, women"), ylab = "Lifetime per 100,000",
              main = "Lifetime EAC, no screening")
text(bp, rbind(x$target_value, x$model_value), round(rbind(x$target_value, x$model_value)), pos = 3, cex = 0.8, xpd = TRUE)
legend("topright", c("CISNET (Rubenstein 2022, from 40)", "Model (from 50)"), fill = c("grey75", "#1F3A5F"), border = NA, bty = "n")
invisible(dev.off())

print(rates, digits = 3); print(tab[, c("target", "kind", "target_value", "model_value", "difference_pct")], digits = 4)
cat("wrote", out, "\n")
