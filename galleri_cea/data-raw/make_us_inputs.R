# Builds the US baseline_risks.csv and cohort.csv from published sources.
# Run from the project folder:  Rscript data-raw/make_us_inputs.R
# Every number below is cited in inputs/SOURCES.md.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

ages <- 50:100

# ---- All-cause death probability qx, SSA Period Life Table 2022 (2025 TR) ----
qx_m <- c(0.005666, 0.006069, 0.006539, 0.007073, 0.007675, 0.008348, 0.009051, 0.009822, 0.010669, 0.011548,
          0.012458, 0.013403, 0.014450, 0.015571, 0.016737, 0.017897, 0.019017, 0.020213, 0.021569, 0.023088,
          0.024828, 0.026705, 0.028761, 0.031116, 0.033861, 0.037088, 0.041126, 0.045241, 0.049793, 0.054768,
          0.060660, 0.067027, 0.073999, 0.081737, 0.090458, 0.100525, 0.111793, 0.124494, 0.138398, 0.153207,
          0.169704, 0.187963, 0.208395, 0.230808, 0.253914, 0.277402, 0.300882, 0.324326, 0.347332, 0.369430,
          0.391927)
qx_f <- c(0.003407, 0.003642, 0.003917, 0.004238, 0.004619, 0.005040, 0.005493, 0.005987, 0.006509, 0.007067,
          0.007658, 0.008305, 0.008991, 0.009681, 0.010343, 0.011018, 0.011743, 0.012532, 0.013512, 0.014684,
          0.016025, 0.017468, 0.019195, 0.021195, 0.023452, 0.025980, 0.029153, 0.032394, 0.035888, 0.039676,
          0.044156, 0.049087, 0.054635, 0.061066, 0.068431, 0.076841, 0.086205, 0.096851, 0.109019, 0.121867,
          0.135805, 0.151108, 0.168020, 0.186340, 0.206432, 0.228086, 0.250406, 0.273699, 0.296984, 0.319502,
          0.342716)

# ---- SEER 21 age-specific incidence per 100,000, 2013-2017 (CSR 1975-2017) ----
band_from <- c(50, 55, 60, 65, 70, 75, 80, 85)            # 85 = 85+
seer <- tibble(
  band_from = rep(band_from, 2), sex = rep(c("male", "female"), each = 8),
  eso   = c(7.5, 14.2, 22.2, 30.3, 37.8, 41.3, 45.6, 40.7,   1.7, 3.4, 4.7, 6.4, 8.6, 10.3, 12.4, 12.8),  # Table 8.7
  stom  = c(10.2, 16.3, 23.4, 32.3, 45.4, 58.7, 68.6, 70.7,  6.5, 8.7, 11.5, 15.6, 21.0, 27.8, 34.7, 36.1), # Table 24.7
  all   = c(508.9, 860.1, 1312.0, 1888.1, 2296.8, 2623.9, 2803.7, 2748.9,
            618.8, 791.6, 1030.8, 1345.2, 1614.5, 1815.0, 1896.5, 1726.0))                              # Table 2.7

# ---- Other sourced constants ----
mst_ugi   <- 2.37    # preclinical mean sojourn time, gastric cancer (Bae 2014)
mst_other <- 2.37    # ASSUMPTION: same sojourn for other cancers the workup can find
eac_share <- 0.631   # adenocarcinoma share of esophageal cancer, SEER 2004-15 (Then 2020)
k_share   <- 0.295   # Kunzmann score >= 8: share of people flagged (Kunzmann 2018)
k_capture <- 0.774   # ... and share of EAC cases they contain
rr_eac    <- c(kunzmann_high = k_capture / k_share, kunzmann_low = (1 - k_capture) / (1 - k_share), all = 1)
# EQ-5D-5L US population norms by age (Jiang 2021), not sex-specific
u_norm <- function(a) case_when(a < 55 ~ 0.816, a < 65 ~ 0.815, a < 75 ~ 0.824, TRUE ~ 0.811)
# WHO GHE 2019 standard YLL per death by age group (85+ is one group)
yll <- c(`50` = 38.10, `55` = 33.33, `60` = 28.66, `65` = 24.12, `70` = 19.76, `75` = 15.65, `80` = 11.96, `85` = 7.05)

baseline <- expand_grid(risk_group = c("all", "kunzmann_high", "kunzmann_low"), sex = c("male", "female"),
                        age_from = ages) |>
  mutate(band_from = pmin(85, 5 * (age_from %/% 5))) |>
  left_join(seer, by = c("band_from", "sex")) |>
  mutate(
    age_to = age_from,
    mort_other = ifelse(sex == "male", qx_m[age_from - 49], qx_f[age_from - 49]),
    # Kunzmann relative risk acts on the EAC part of esophageal cancer only
    rr = rr_eac[risk_group],
    inc_eso = eso * (eac_share * rr + (1 - eac_share)),
    inc_ugi = (inc_eso + stom) / 1e5,
    prev_ugi = inc_ugi * mst_ugi,
    prev_other = (all - eso - stom) / 1e5 * mst_other,
    utility_well = u_norm(age_from),
    daly_ref_le = unname(yll[as.character(band_from)]),
    status = "SOURCED",
    source = "qx: SSA 2022 period life table; incidence: SEER CSR 1975-2017 Tables 2.7/8.7/24.7; see SOURCES.md") |>
  mutate(across(c(mort_other, inc_ugi, prev_ugi, prev_other), ~ signif(.x, 4))) |>
  select(risk_group, sex, age_from, age_to, mort_other, inc_ugi, prev_ugi, prev_other,
         utility_well, daly_ref_le, status, source)
baseline <- mutate(baseline, across(c(mort_other, inc_ugi, prev_ugi, prev_other),
                                    ~ format(.x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)))
write_csv(baseline, "inputs/baseline_risks.csv", na = "")

# ---- Cohort: US adults 50-79, Census Vintage 2024 (NC-EST2024-SYASEXN) ----
pop <- tibble(band_from = rep(seq(50, 75, 5), 2), sex = rep(c("male", "female"), each = 6),
              n = c(10192640, 10051985, 10363023, 9297585, 7402084, 5435868,
                    10293667, 10307943, 10938774, 10191886, 8553407, 6557458)) |>
  mutate(share = n / sum(n), start_age = band_from + 2)

screening <- expand_grid(pop, tibble(risk_group = c("kunzmann_high", "kunzmann_low"),
                                     k = c(k_share, 1 - k_share), kunzmann_scoped = c(1, 0))) |>
  mutate(population = "screening", share = share * k, status = "ASSUMPTION",
         notes = "Census 2024 age-sex shares x Kunzmann 29.5% flagged (same in every age-sex group: assumption)")
gpos <- pop |> mutate(population = "galleri_positive", risk_group = "all", kunzmann_scoped = 0, status = "ASSUMPTION",
                      notes = "Galleri +ve cohort assumed to mirror the 50-79 age-sex mix; replace with your own cohort")
cohort <- bind_rows(screening, gpos) |>
  mutate(share = round(share, 6)) |>
  group_by(population) |> mutate(share = round(share + (row_number() == 1) * (1 - sum(share)), 6)) |> ungroup() |>
  select(population, risk_group, sex, start_age, share, kunzmann_scoped, status, notes)
write_csv(cohort, "inputs/cohort.csv", na = "")
message("wrote inputs/baseline_risks.csv (", nrow(baseline), " rows) and inputs/cohort.csv (", nrow(cohort), " rows)")
