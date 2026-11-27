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
             "tb", "lap", "surv", "egd_pos", "egd_neg", "egd_skipped", "found_pet_ugi", "hidden_new", "missed_visible",
             "incid_free", "perf", "perf_deaths", "proc_deaths", "missed_inv", "late_screen")
acc_keys <- c("clin_dx", "clin_late", "ca_death", "oc_death", "ly", "symp_perf")
by_age <- list()
for (arm in E$ARMS) {
  tot <- NULL
  for (i in seq_len(nrow(C))) {
    row <- C[i, ]; prep <- E$prepare_row(row, D, p, S)
    s <- E$simulate_arm(arm, row, prep$curve, prep$cal, prep$p, S, wk, prep$V, prep$q)
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

# workup flows for everyone worked up (base case, lifetime), from the model's own counts
g <- colSums(by_age$Galleri[, ev_keys])
s_egd <- 1 - p$p_perf_egd * p$p_death_perf
flow <- function(from, to, v) data.frame(source = from, target = to, value = v)
egd_scr <- g[["investigated"]] - g[["egd_skipped"]]
rest <- g[["egd_neg"]] * s_egd + g[["egd_skipped"]]
pb <- g[["pet"]] * p$p_pet_pb; inc <- g[["incid"]]
links <- rbind(
  flow("Galleri +ve", "EGD", egd_scr),
  flow("Galleri +ve", "EGD skipped (hold)", g[["egd_skipped"]]),
  flow("EGD", "EGD +ve: cancer found", g[["egd_pos"]]),
  flow("EGD", "EGD -ve", g[["egd_neg"]] * s_egd),
  flow("EGD", "Fatal perforation", g[["egd_neg"]] * (1 - s_egd)),
  flow("EGD -ve", "Rest of workup", g[["egd_neg"]] * s_egd),
  flow("EGD skipped (hold)", "Rest of workup", g[["egd_skipped"]]),
  flow("Rest of workup", "ENT", g[["ent"]]),
  flow("Rest of workup", "PET CT", g[["pet"]]),
  flow("ENT", "ENT +ve: cancer", g[["ent"]] * p$p_ent_pos),
  flow("ENT", "ENT -ve", g[["ent"]] * (1 - p$p_ent_pos)),
  flow("PET CT", "Pancreatobiliary +ve", pb),
  flow("PET CT", "Incidental finding", inc),
  flow("PET CT", "Upper GI cancer on PET", g[["found_pet_ugi"]]),
  flow("PET CT", "PET -ve", g[["pet"]] - pb - inc - g[["found_pet_ugi"]]),
  flow("Pancreatobiliary +ve", "EUS sample +ve: cancer", pb - g[["lap"]] - g[["surv"]]),
  flow("Pancreatobiliary +ve", "Diagnostic laparoscopy", g[["lap"]]),
  flow("Pancreatobiliary +ve", "Watchful waiting", g[["surv"]]),
  flow("Incidental finding", "Incidental: cancer", inc * p$p_incid_cancer),
  flow("Incidental finding", "Incidental: benign", inc * (1 - p$p_incid_cancer)))
links <- links[links$value > 0, ]

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
