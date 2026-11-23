# =============================================================
# Data for the visual stories (viz/). Run after Rscript run.R:
#   Rscript tools/viz_data.R          -> viz/galleri_viz_data.json
# Base case, per 100,000, by single year of age for every arm, plus the
# workup flows, PSA draws, scenarios, one-way SA and subgroup grid.
# =============================================================
suppressPackageStartupMessages({ library(ceakit); library(jsonlite) })
proj <- load_project(".")
E <- proj$env; p <- param_list(proj$params); S <- proj$settings; D <- proj$data
C <- E$cohort_rows(D); N <- 1e5
wk <- E$workup_summary(p)

ev_keys <- c("tests", "tp", "fp", "fp_other_ca", "tn", "fn", "fn_symp", "found_ugi", "found_other",
             "overdx_ugi", "overdx_other", "investigated", "free_inv", "egd", "ent", "pet", "eus", "samp", "incid",
             "incid_free", "perf", "perf_deaths", "proc_deaths", "missed_inv", "late_screen")
acc_keys <- c("clin_dx", "clin_late", "ca_death", "oc_death", "ly", "symp_perf")
by_age <- list()
for (arm in E$ARMS) {
  tot <- NULL
  for (i in seq_len(nrow(C))) {
    row <- C[i, ]; prep <- E$prepare_row(row, D, p, S)
    s <- E$simulate_arm(arm, row, prep$curve, prep$cal, p, S, wk, prep$V, prep$q)
    n <- nrow(s$ev_t)
    m <- cbind(s$ev_t[, ev_keys], s$acc_t[, acc_keys],
               alive = rowSums(s$trace[seq_len(n), !E$STATES %in% E$DEAD]),
               undiagnosed = s$trace[seq_len(n), "pre"] + s$trace[seq_len(n), "pre_m"],
               dead_cancer = s$trace[seq_len(n) + 1, "dead_cancer"]) * row$share * N
    tot <- if (is.null(tot)) m else tot + m
  }
  df <- as.data.frame(tot); df$age <- C$start_age[1] + seq_len(nrow(df)) - 1
  by_age[[arm]] <- df
}

# workup flows for everyone worked up (base case, lifetime)
d <- E$workup_paths(p); s_egd <- 1 - p$p_perf_egd * p$p_death_perf
inv <- sum(by_age$Galleri$investigated)
flow <- function(from, to, v) data.frame(source = from, target = to, value = v)
pp <- function(leaf) sum(d$prob[d$leaf %in% leaf])
neg <- inv * s_egd * (1 - p$p_egd_pos)
links <- rbind(
  flow("Galleri +ve", "EGD", inv),
  flow("EGD", "EGD +ve (true positive)", inv * s_egd * p$p_egd_pos),
  flow("EGD", "EGD -ve (false positive)", neg),
  flow("EGD", "Fatal perforation", inv * (1 - s_egd)),
  flow("EGD -ve (false positive)", "ENT", neg * p$p_ent_route),
  flow("EGD -ve (false positive)", "PET CT", neg * (1 - p$p_ent_route)),
  flow("ENT", "ENT +ve: cancer", neg * p$p_ent_route * p$p_ent_pos),
  flow("ENT", "ENT -ve", neg * p$p_ent_route * (1 - p$p_ent_pos)),
  flow("PET CT", "Pancreatobiliary +ve", neg * (1 - p$p_ent_route) * p$p_pet_pb),
  flow("PET CT", "Incidental finding", neg * (1 - p$p_ent_route) * p$p_pet_incid),
  flow("PET CT", "PET -ve", neg * (1 - p$p_ent_route) * (1 - p$p_pet_pb - p$p_pet_incid)),
  flow("Pancreatobiliary +ve", "EUS sampling +ve: cancer", neg * (1 - p$p_ent_route) * p$p_pet_pb * p$p_eus_pos),
  flow("Pancreatobiliary +ve", "EUS sampling -ve", neg * (1 - p$p_ent_route) * p$p_pet_pb * p$p_eus_neg),
  flow("Pancreatobiliary +ve", "No sampling", neg * (1 - p$p_ent_route) * p$p_pet_pb * (1 - p$p_eus_pos - p$p_eus_neg)),
  flow("Incidental finding", "Incidental: cancer", neg * (1 - p$p_ent_route) * p$p_pet_incid * p$p_incid_cancer),
  flow("Incidental finding", "Incidental: benign", neg * (1 - p$p_ent_route) * p$p_pet_incid * (1 - p$p_incid_cancer)))

rl <- function(a) { x <- read.csv(file.path("outputs", a, "results_long.csv")); x[, c("key", "arm", "value", "lower_95", "upper_95")] }
psa_c <- read.csv("outputs/galleri/psa_cost.csv", check.names = FALSE)
psa_e <- read.csv("outputs/galleri/psa_effect.csv", check.names = FALSE)
out <- list(
  meta = list(cohort = "50-year-old non-Hispanic White men", per = N, start_age = C$start_age[1],
              stop_age = S$screen_stop_age, max_age = S$max_age, cost_year = S$cost_year, wtp = S$wtp,
              c_galleri = p$c_galleri, generated = format(Sys.time(), "%Y-%m-%d")),
  by_age = by_age, workup_links = links,
  galleri = rl("galleri"), comparison = rl("comparison"),
  psa = data.frame(d_cost = psa_c$Galleri - psa_c$`Usual care`, d_qaly = psa_e$Galleri - psa_e$`Usual care`),
  scenarios = read.csv("outputs/galleri/scenarios.csv"),
  owsa = read.csv("outputs/galleri/owsa.csv"),
  subgroups = subset(read.csv("outputs/galleri/subgroups.csv"), arm == "Galleri"),
  accuracy = read.csv("outputs/galleri/accuracy_by_group.csv"),
  params = data.frame(parameter = proj$params$parameter, label = proj$params$label, value = proj$params$value))
dir.create("viz", showWarnings = FALSE)
write_json(out, "viz/galleri_viz_data.json", digits = 6, auto_unbox = TRUE, na = "null")
message("wrote viz/galleri_viz_data.json")
