# Builds the DUMMY baseline_risks.csv and cohort.csv templates.
# Replace the outputs with real data; this script only documents how the
# placeholders were made. Run from the app folder: Rscript data-raw/make_dummy_inputs.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

# Kunzmann et al. 2018 (UK Biobank): score >= 8 flags 29.5% of people and
# contains 77.4% of esophageal adenocarcinoma cases. Used here to split the
# dummy cancer risk between the two Kunzmann risk groups.
k_share <- 0.295; k_capture <- 0.774
group_rr <- tibble(
  risk_group = c("all", "kunzmann_high", "kunzmann_low"),
  rr_ugi     = c(1, k_capture / k_share, (1 - k_capture) / (1 - k_share)))

bands <- tibble(age_from = seq(40, 95, by = 5)) |>
  mutate(age_to = ifelse(age_from == 95, 100, age_from + 4), age_mid = age_from + 2)

baseline <- expand_grid(group_rr, sex = c("male", "female"), bands) |>
  mutate(
    sex_rr_mort = ifelse(sex == "male", 1, 0.7),
    sex_rr_ugi  = ifelse(sex == "male", 1, 0.5),
    mort_other  = pmin(0.6, 0.0012 * exp(0.09 * (age_mid - 40)) * sex_rr_mort),
    inc_ugi     = 0.0002 * exp(0.06 * (age_mid - 60)) * sex_rr_ugi * rr_ugi,
    prev_ugi    = 0.0015 * exp(0.05 * (age_mid - 60)) * sex_rr_ugi * rr_ugi,
    prev_other  = 0.0005 * exp(0.05 * (age_mid - 60)),
    utility_well = pmax(0.6, 0.90 - 0.004 * (age_mid - 40)),
    daly_ref_le = pmax(2, 88.9 - 0.95 * age_mid),
    status = "DUMMY", source = "") |>
  mutate(across(c(mort_other, inc_ugi, prev_ugi, prev_other), ~ signif(.x, 4)),
         across(c(utility_well, daly_ref_le), ~ round(.x, 3))) |>
  select(risk_group, sex, age_from, age_to, mort_other, inc_ugi, prev_ugi,
         prev_other, utility_well, daly_ref_le, status, source)
baseline <- mutate(baseline, across(c(mort_other, inc_ugi, prev_ugi, prev_other), ~ format(.x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)))
write_csv(baseline, "inputs/baseline_risks.csv", na = "")

cohort <- bind_rows(
  tibble(population = "screening", risk_group = rep(c("kunzmann_high", "kunzmann_low"), each = 2),
         sex = rep(c("male", "female"), 2), start_age = 60,
         share = round(c(k_share, k_share, 1 - k_share, 1 - k_share) / 2, 5),
         kunzmann_scoped = rep(c(1, 0), each = 2)),
  tibble(population = "galleri_positive", risk_group = "all", sex = c("male", "female"),
         start_age = 65, share = c(0.6, 0.4), kunzmann_scoped = 0)) |>
  mutate(status = "DUMMY", notes = "")
write_csv(cohort, "inputs/cohort.csv", na = "")
