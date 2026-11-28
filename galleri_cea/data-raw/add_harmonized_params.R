# One-off, reproducible edit of inputs/parameters.csv and inputs/references.csv for the
# harmonized model: MISCAN-ESO base case, within-year onset timing, Bayesian EGD hold, and the
# patient-time and capacity inputs of the net burden analysis. Safe to rerun (idempotent).
Sys.setlocale("LC_CTYPE", "C.UTF-8")
P <- read.csv("inputs/parameters.csv", colClasses = "character", check.names = FALSE)
setv <- function(k, v) P$value[P$parameter == k] <<- v
setv("miscan_eac", "1")
P$label <- sub(" \\(MISCAN-ESO scenario\\)", " (MISCAN-ESO)", P$label)
row <- function(parameter, label, value, status, unit, description, source, refs, analysis = "both",
                dist = "fixed", low = "", high = "")
  data.frame(parameter, label, value, status, dist, se = "", n = "", group = "", low, high, unit, description,
             source, analysis, refs, price_value = "", price_year = "", check.names = FALSE)
new <- rbind(
  row("within_year_onset", "Onset spread within the year (switch)", "1", "ASSUMPTION", "switch",
      "1 = cancer onsets happen at a uniform time within each yearly cycle, so some present (or progress) before the next annual test, as in continuous time; 0 = all onsets at the cycle end (the original cohort timing)",
      "Matches the microsimulation's continuous-time natural history; removes most of the cohort-microsim gap in screen-found cancers", ""),
  row("egd_rule", "EGD decision rule after a negative EGD (switch)", "1", "ASSUMPTION", "switch",
      "1 = Bayesian: after a negative EGD, a new positive Galleri gets an EGD when P(upper GI cancer | history) >= egd threshold; 0 = fixed hold of egd_hold_years rounds",
      "Posterior from usual-care prevalence, regrowth after a negative EGD (preclinical sojourn), Galleri sensitivity and the person's false-positive history (beta-binomial)", "7, 43"),
  row("egd_threshold", "EGD probability threshold", "-1", "ASSUMPTION", "probability",
      "Probability of upper GI cancer above which the EGD is done; -1 = set by the model at each age from QALYs: L / (L + G), G = QALYs gained by finding the cancer now, L = QALYs an EGD costs a person without cancer",
      "Decision-analytic threshold (Pauker-Kassirer); 0.3 and 0.5 tested as scenarios", ""),
  row("egd_memory_years", "Years a negative EGD is remembered (Bayesian rule)", "6", "ASSUMPTION", "years",
      "Rounds after a negative EGD over which the Bayesian rule is applied; after that the EGD is always done",
      "About one EAC sojourn plus a year (Kong 2014: 4.5-4.9 y); the posterior is near its unconditional value by then", "43"),
  row("pt_galleri", "Patient time: Galleri blood draw", "1.5", "ASSUMPTION", "hours",
      "Hours of patient time per Galleri test (travel, waiting, draw)", "Jonas 2008: travel 1.1 h per visit; draw and wait assumed 0.4 h", "57", "burden", "fixed", "1", "3"),
  row("pt_egd", "Patient time: EGD", "8", "ASSUMPTION", "hours",
      "Hours per sedated EGD incl. fasting, travel, time at the centre and recovery (a day off, usually with a companion)",
      "Jonas 2008 (colonoscopy, sedated): travel 1.1 h, at the centre 2.8 h, recovery 0.8 h; no bowel prep for EGD; rounded up to a working day", "57", "burden", "fixed", "5", "12"),
  row("pt_eus", "Patient time: EUS (with or without sampling)", "8", "ASSUMPTION", "hours", "Hours per sedated EUS", "As pt_egd (Jonas 2008)", "57", "burden", "fixed", "5", "12"),
  row("pt_pet", "Patient time: PET CT", "4", "ASSUMPTION", "hours", "Hours per PET CT incl. uptake time and travel", "Assumed: 2 h at the scanner plus travel (Jonas 2008)", "57", "burden", "fixed", "3", "6"),
  row("pt_ent", "Patient time: ENT evaluation", "2.5", "ASSUMPTION", "hours", "Hours per ENT visit with nasopharyngolaryngoscopy", "Assumed: one clinic visit plus travel (Jonas 2008)", "57", "burden", "fixed", "1.5", "4"),
  row("pt_lap", "Patient time: diagnostic laparoscopy", "24", "ASSUMPTION", "hours", "Hours per diagnostic laparoscopy incl. day-case stay and recovery", "Assumed: day surgery and 2 days of recovery", "", "burden", "fixed", "16", "48"),
  row("pt_other", "Patient time: other visit (imaging, tumor board result, follow-up)", "3", "ASSUMPTION", "hours", "Hours per other workup contact (incidental-finding imaging, watchful-waiting scan, sampling)", "Assumed (Jonas 2008 travel and waiting)", "57", "burden", "fixed", "2", "5"),
  row("us_pop_50plus", "US population aged 50+", "123300000", "SOURCED", "people", "Population used to scale per-100,000 results to the US", "US Census Bureau, Vintage 2024 population estimates", "60", "burden"),
  row("us_upper_endo", "US upper endoscopies per year", "7459419", "SOURCED", "procedures", "Annual upper GI endoscopies (EGD) in the US: the capacity reference", "Peery 2022 (Gastroenterology update 2021), 2019 estimate", "58", "burden"),
  row("us_eus", "US upper EUS per year", "290655", "SOURCED", "procedures", "Annual upper endoscopic ultrasounds in the US", "Peery 2022, 2019 estimate", "58", "burden")
)
P <- rbind(P[!P$parameter %in% new$parameter, ], new)
P[P == ""] <- NA
data.table::fwrite(P, "inputs/parameters.csv", na = "")
R <- read.csv("inputs/references.csv", colClasses = "character")
add <- data.frame(ref = c("57", "58", "59", "60"),
  citation = c("Jonas DE, Russell LB, Sandler RS, Chou J, Pignone M. Value of patient time invested in the colonoscopy screening process: time requirements for colonoscopy study. Med Decis Making. 2008;28(1):56-65.",
               "Peery AF, Crockett SD, Murphy CC, et al. Burden and cost of gastrointestinal, liver, and pancreatic diseases in the United States: update 2021. Gastroenterology. 2022;162(2):621-644.",
               "Peery AF, Murphy CC, Anderson C, et al. Burden and cost of gastrointestinal, liver, and pancreatic diseases in the United States: update 2024. Gastroenterology. 2025;168(5):1000-1024.",
               "US Census Bureau. Annual estimates of the resident population by single year of age and sex for the United States, Vintage 2024."),
  url = c("https://doi.org/10.1177/0272989X07309643", "https://doi.org/10.1053/j.gastro.2021.10.017",
          "https://doi.org/10.1053/j.gastro.2024.12.029", "https://www.census.gov/programs-surveys/popest.html"),
  pmid = c("18263561", "34678215", "39920892", ""))
R <- rbind(R[!R$ref %in% add$ref, ], add)
R[R == ""] <- NA
data.table::fwrite(R, "inputs/references.csv", na = "")

# ---- scenarios: the new base case (MISCAN-ESO, within-year onset, Bayesian EGD rule) -------
S <- read.csv("inputs/scenarios.csv", colClasses = "character", check.names = FALSE)
S <- S[!(S$parameter == "miscan_eac" & S$value == "1"), ]                    # now the base case
S$scenario[S$scenario == "MISCAN-ESO natural history for esophageal adenocarcinoma"] <- NA
S <- S[!is.na(S$scenario), ]
fixed_rule <- S$scenario %in% c("No EGD hold (EGD at every positive)", "EGD hold 3 years")
S <- rbind(S, data.frame(scenario = unique(S$scenario[fixed_rule]), parameter = "egd_rule", value = "0", refs = NA))
site <- unique(S$scenario[S$parameter == "site_sojourn"])                    # need the original EAC history
S <- rbind(S, data.frame(scenario = site, parameter = "miscan_eac", value = "0", refs = NA))
add <- data.frame(
  scenario = c("Original natural history (one sojourn for all upper GI cancer)", "Fixed 2-year EGD hold (original rule)",
               "Fixed 2-year EGD hold (original rule)", "Bayesian EGD rule, threshold 0.3", "Bayesian EGD rule, threshold 0.5",
               "Yearly onset timing (original cohort timing)"),
  parameter = c("miscan_eac", "egd_rule", "egd_hold_years", "egd_threshold", "egd_threshold", "within_year_onset"),
  value = c("0", "0", "2", "0.3", "0.5", "0"), refs = c("7", "7, 43", "7, 43", NA, NA, NA))
S <- rbind(S[!S$scenario %in% add$scenario, ], add)
S[S == ""] <- NA
data.table::fwrite(S, "inputs/scenarios.csv", na = "")
