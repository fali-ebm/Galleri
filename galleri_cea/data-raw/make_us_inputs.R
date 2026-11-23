# Builds inputs/baseline_risks.csv and inputs/cohort.csv from published US data.
# Run from the project folder:  Rscript data-raw/make_us_inputs.R
# Raw inputs (data-raw/raw/) and their sources are listed in inputs/SOURCES.md.
#
# baseline_risks.csv: one row per race/ethnicity x sex x single year of age 50-99
#   mort_other  all-cause death probability qx (NCHS US life tables 2024)
#   inc_eso / inc_stom / inc_other  SEER clinical incidence per person-year
#               (2013-2017). These are CALIBRATION TARGETS: the model fits
#               cancer onset so that symptomatic presentation reproduces them.
#   eac_share   adenocarcinoma share of esophageal cancer (Kunzmann acts on it)
#   late_eso / late_stom  share of staged symptomatic diagnoses that are
#               regional or distant (SEER 2010-2016, CSR Tables 8.8 / 24.8):
#               stage-at-diagnosis CALIBRATION TARGETS by site, race and sex.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })
raw <- function(f) read_csv(file.path("data-raw/raw", f), show_col_types = FALSE)

RACES <- c(all = "All races/ethnicities", nhw = "Non-Hispanic White", nhb = "Non-Hispanic Black",
           nhapi = "Non-Hispanic Asian/Pacific Islander", nhaian = "Non-Hispanic American Indian/Alaska Native",
           hisp = "Hispanic")
ages <- 50:99

# ---- SEER age-specific rates (5-year bands) -> single years ------------------
# Log-linear interpolation between band midpoints; flat beyond the first and
# last midpoints (85+ is treated as one band centred at 87.5).
seer <- raw("seer_csr_2013_2017.csv")
interp <- function(band, rate) exp(approx(band + 2.5, log(rate), xout = ages + 0.5, rule = 2)$y) / 1e5
single <- seer |> group_by(site, race, sex) |> reframe(age = ages, rate = interp(age_band, rate))

# ---- Race/ethnicity scaling (SEER CSR tables 8.15, 24.15, 2.15) ---------------
# Age-specific rates exist only for all races, White and Black. Other groups
# use a proportional assumption: the all-race (or White) age curve x the ratio
# of age-adjusted rates.
aa <- raw("seer_csr_age_adjusted_by_race.csv") |> pivot_longer(c(male, female), names_to = "sex", values_to = "aa")
ratio <- function(site_, sex_, num, den) {
  a <- aa[aa$site == site_ & aa$sex == sex_, ]
  a$aa[a$race == num] / a$aa[a$race == den]
}
curve_for <- function(site_, sex_, g) {
  base <- function(r) single$rate[single$site == site_ & single$race == r & single$sex == sex_]
  switch(g,
    all    = base("all"),
    nhw    = base("white") * ratio(site_, sex_, "nhw", "white"),
    nhb    = base("black"),                                   # Black (all ethnicities) as proxy
    nhapi  = base("all") * ratio(site_, sex_, "api", "all"),
    nhaian = base("all") * ratio(site_, sex_, "aian", "all"),
    hisp   = base("all") * ratio(site_, sex_, "hispanic", "all"))
}

lt <- raw("nchs_lifetables_2024_qx.csv")

# SEER stage at diagnosis by site, race and sex (all races / White / Black only):
# NHW uses White, NHB uses Black, other groups use all races. ASSUMPTION: no
# variation by age within a race-sex group (the site mix still varies by age).
stage <- raw("seer_stage_by_race_sex.csv") |>
  mutate(late = (regional + distant) / (localized + regional + distant))
stage_race <- c(all = "all", nhw = "white", nhb = "black", nhapi = "all", nhaian = "all", hisp = "all")
late_for <- function(site_, g, s) stage$late[stage$site == site_ & stage$race == stage_race[[g]] & stage$sex == s]
eac <- raw("seer_eac_share_by_race.csv")

# EQ-5D-5L US norms by age (Jiang 2021); WHO GHE 2019 standard YLL per death
u_norm <- function(a) case_when(a < 55 ~ 0.816, a < 65 ~ 0.815, a < 75 ~ 0.824, TRUE ~ 0.811)
yll <- c(`50` = 38.10, `55` = 33.33, `60` = 28.66, `65` = 24.12, `70` = 19.76, `75` = 15.65, `80` = 11.96, `85` = 7.05)

grid <- expand_grid(race_eth = names(RACES), sex = c("male", "female"))
baseline <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
  g <- grid$race_eth[i]; s <- grid$sex[i]
  tibble(race_eth = g, sex = s, age_from = ages,
         inc_eso = curve_for("esophagus", s, g), inc_stom = curve_for("stomach", s, g),
         inc_all = curve_for("all_sites", s, g),
         late_eso = late_for("esophagus", g, s), late_stom = late_for("stomach", g, s))
})) |>
  left_join(rename(lt, age_from = age, mort_other = qx), by = c("race_eth", "sex", "age_from")) |>
  left_join(rename(eac, race_eth = race), by = "race_eth") |>
  mutate(age_to = age_from,
         inc_other = pmax(0, inc_all - inc_eso - inc_stom),
         utility_well = u_norm(age_from),
         daly_ref_le = unname(yll[as.character(pmin(85, 5 * (age_from %/% 5)))]),
         status = "SOURCED",
         source = "SEER CSR 1975-2017 Tables 2.7/8.7/24.7 (+2.15/8.15/24.15 race ratios), 8.8/24.8 (stage); NCHS US Life Tables 2024") |>
  mutate(across(c(inc_eso, inc_stom, inc_other), ~ signif(.x, 5)), across(c(late_eso, late_stom), ~ round(.x, 4))) |>
  select(race_eth, sex, age_from, age_to, mort_other, inc_eso, inc_stom, inc_other, eac_share, late_eso, late_stom,
         utility_well, daly_ref_le, status, source)
stopifnot(!anyNA(baseline$mort_other), !anyNA(baseline$inc_eso), nrow(baseline) == 6 * 2 * length(ages))
write_csv(baseline, "inputs/baseline_risks.csv", na = "")

# ---- Base-case cohort: 50-year-old non-Hispanic White men --------------------
# Split by Kunzmann score (29.5% score >= 8; Kunzmann 2018). The subgroup grid
# (age x sex x race/ethnicity) is built by the model from settings$subgroups.
k_share <- 0.295
cohort <- tibble(population = "screening", race_eth = "nhw", sex = "male", start_age = 50,
                 risk_group = c("kunzmann_high", "kunzmann_low"), share = c(k_share, 1 - k_share),
                 kunzmann_scoped = c(1, 0), status = "ASSUMPTION",
                 notes = "Base case: 50-year-old non-Hispanic White men; Kunzmann 29.5% flagged (same in every subgroup)")
write_csv(cohort, "inputs/cohort.csv", na = "")
message("wrote inputs/baseline_risks.csv (", nrow(baseline), " rows) and inputs/cohort.csv")
