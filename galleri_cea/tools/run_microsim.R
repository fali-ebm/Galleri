# Microsimulation run: hesim natural history + per-person screening overlay.
#   Rscript tools/run_microsim.R [n_people = 1e6] [seed = 2026]
# Writes outputs/microsim/: validation.csv (vs the cohort model), summary.csv,
# distributions.csv, per_person_increments.csv, patients_sample.csv (10,000 people x
# arms), patients.csv.gz (everyone) and microsim_distributions.png.
args <- commandArgs(TRUE)
N <- if (length(args) >= 1) as.numeric(args[1]) else 1e6
seed <- if (length(args) >= 2) as.integer(args[2]) else 2026L
suppressPackageStartupMessages({ library(data.table) })
proj <- ceakit::load_project("."); E <- proj$env; p <- ceakit::param_list(proj$params)
S <- proj$settings; D <- proj$data; arms <- E$ARMS
out_dir <- "outputs/microsim"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
C <- E$kunz_cohort(E$cohort_rows(D), E$kunz_settings(p, S))
t0 <- Sys.time()
INC <- list()
X <- rbindlist(lapply(seq_len(nrow(C)), function(i) {
  n <- round(N * C$share[i]); cat(sprintf("%s: %s people\n", C$risk_group[i], format(n, big.mark = ",")))
  r <- E$ms_row(C[i, ], D, p, S, arms, n, seed + 17L * i); INC[[i]] <<- attr(r, "incidence"); r
}))
# natural-history check: usual-care symptomatic incidence vs the SEER target, 5-year bands
I <- as.data.table(do.call(rbind, INC))[, band := 5 * floor(age / 5)]
I <- I[, .(microsim_per_100k = 1e5 * sum(events) / sum(alive), seer_target_per_100k = 1e5 * sum(target * alive) / sum(alive),
           events = sum(events)), by = band]
fwrite(I, file.path(out_dir, "incidence_check.csv")); print(I)
cat("simulated in", format(round(difftime(Sys.time(), t0, units = "secs"))), "\n")
X[, person := .GRP, by = .(risk_group, person)]
n_people <- uniqueN(X$person)

# ---- validation against the cohort model (same parameters) ------------------------
coh <- rbind(fread("outputs/galleri/results_long.csv")[arm %in% c("Galleri", "Usual care")],
             fread("outputs/comparison/results_long.csv")[arm == "Kunzmann -> EGD"])
cv <- function(k, a) { v <- coh[key == k & arm == a, value]; if (length(v)) v[1] else NA_real_ }
per <- 1e5
PER100K <- c("ly_total", "qaly_total", "daly_total")    # cohort stores these per 100,000; cost per person
mc <- function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)))
X[, clin_net := clin_dx - found_other]   # cohort clin_dx_lifetime is net of other-site cancers found early
map <- list(   # microsim column, cohort key, scale
  c("stage12", "stage12_ugi", per), c("stage34", "stage34_ugi", per), c("screen_stage12", "screen_stage12_ugi", per),
  c("screen_detected", "found_upper_gi", per), c("found_other", "found_other", per),
  c("clin_net", "clin_dx_lifetime", per), c("overdiagnosed", "overdiagnosed_ugi", per),
  c("ca_death", "ca_deaths", per), c("false_pos", "free_investigated", per), c("egd", "egd", per),
  c("egd_skipped", "egd_skipped_hold", per), c("hidden", "hidden_after_egd", per), c("pet", "pet_ct", per),
  c("eus", "eus", per), c("laparoscopy", "laparoscopy", per), c("watchful_waiting", "surveillance_imaging", per),
  c("perforation", "perforations", per), c("ly", "ly_total", 1), c("qaly", "qaly_total", 1),
  c("cost", "cost_total", 1), c("daly", "daly_total", 1), c("qaly_fp", "qaly_fp", per), c("daly_fp", "daly_fp", per))
V <- rbindlist(lapply(arms, function(a) rbindlist(lapply(map, function(m) {
  x <- X[arm == a][[m[1]]]; s <- as.numeric(m[3]); r <- mc(x)
  data.table(arm = a, measure = m[1], cohort_key = m[2], microsim = r[["mean"]] * s,
             mc_95 = 1.96 * r[["se"]] * s, cohort = cv(m[2], a) / if (m[2] %in% PER100K) per else 1)
}))))
V[, diff_pct := 100 * (microsim - cohort) / abs(cohort)]
V[, within_mc := abs(microsim - cohort) <= mc_95]
fwrite(V, file.path(out_dir, "validation.csv"))

# ---- per-person increments (same people, Galleri vs usual care; Kunzmann vs usual care)
W <- dcast(X, person + risk_group ~ arm, value.var = c("qaly", "cost", "ly", "ca_death", "screen_detected", "stage_late", "dx_time", "stage12", "stage34"))
setnames(W, gsub("[^A-Za-z0-9_]+", "_", names(W)))
W[, `:=`(dq_g = qaly_Galleri - qaly_Usual_care, dc_g = cost_Galleri - cost_Usual_care,
         dq_k = qaly_Kunzmann_EGD - qaly_Usual_care, dc_k = cost_Kunzmann_EGD - cost_Usual_care)]
inc <- function(dq, dc) c(qaly = mean(dq), qaly_mc95 = 1.96 * sd(dq) / sqrt(length(dq)),
                         cost = mean(dc), cost_mc95 = 1.96 * sd(dc) / sqrt(length(dc)), icer = mean(dc) / mean(dq),
                         share_gain = mean(dq > 1e-9), share_loss = mean(dq < -1e-9), share_same = mean(abs(dq) <= 1e-9))
s12 <- function(a, b) { d <- as.numeric(a) - as.numeric(b); c(stage12_gained_per100k = 1e5 * mean(d), stage12_gained_mc95 = 1e5 * 1.96 * sd(d) / sqrt(length(d))) }
ICR <- rbind(data.table(comparison = "Galleri vs usual care", t(inc(W$dq_g, W$dc_g)), t(s12(W$stage12_Galleri, W$stage12_Usual_care))),
             data.table(comparison = "Kunzmann -> EGD vs usual care", t(inc(W$dq_k, W$dc_k)), t(s12(W$stage12_Kunzmann_EGD, W$stage12_Usual_care))))
ICR[, cost_per_stage12_gained := 1e5 * cost / stage12_gained_per100k]
ICR[, cohort_icer := c(cv("icer", "Galleri"), coh[key == "icer" & arm == "Kunzmann -> EGD", value][1])]
fwrite(ICR, file.path(out_dir, "summary.csv"))

# ---- distributions across people (Galleri arm unless stated) ----------------------
G <- X[arm == "Galleri"]
bucket <- function(x, top) { x <- pmin(x, top); f <- table(factor(x, levels = 0:top)) / length(x)
  data.table(value = c(as.character(0:(top - 1)), paste0(top, "+")), share = as.numeric(f)) }
dist <- rbindlist(list(
  cbind(measure = "Positive Galleri tests per person", bucket(G$positives, 4)),
  cbind(measure = "False positives per person (cancer-free)", bucket(G$false_pos, 3)),
  cbind(measure = "EGDs per person", bucket(G$egd, 4)),
  cbind(measure = "PET CT scans per person", bucket(G$pet, 3)),
  cbind(measure = "Laparoscopies per person", bucket(G$laparoscopy, 2)),
  cbind(measure = "Positives skipped by the EGD hold, per person", bucket(G$egd_skipped, 2))))
fwrite(dist, file.path(out_dir, "distributions.csv"))
# among people whose upper GI cancer is found by Galleri: lead time vs their own usual-care diagnosis
L <- W[screen_detected_Galleri == TRUE & !is.na(dx_time_Usual_care)]
lead <- L$dx_time_Usual_care - L$dx_time_Galleri
stage <- W[screen_detected_Galleri == TRUE & !is.na(stage_late_Usual_care)]
pp <- data.table(
  measure = c("People simulated", "Upper GI cancer in lifetime (usual care), %", "Galleri finds it, % of those cancers",
              "Lead time, median years (screen-detected, presented in usual care)", "Lead time, IQR",
              "Late stage in usual care -> early when found by Galleri, % of screen-detected",
              "Overdiagnosed, % of screen-detected", "People with any positive test, %",
              "People with >= 2 false positives, %", "Galleri: share of people better off (QALYs), %",
              "Galleri: share worse off, %", "Galleri: largest single-person QALY gain", "Galleri: largest single-person QALY loss",
              "Cost per person, 90th percentile (Galleri)", "Cost per person, 99th percentile (Galleri)"),
  value = c(n_people, 100 * mean(X[arm == "Usual care", cancer]),
            100 * mean(G$screen_detected[X[arm == "Usual care", cancer]]),
            median(lead), paste(round(quantile(lead, c(.25, .75)), 2), collapse = "-"),
            100 * mean(stage$stage_late_Usual_care & !stage$stage_late_Galleri),
            100 * mean(G$overdiagnosed[G$screen_detected]), 100 * mean(G$positives > 0), 100 * mean(G$false_pos >= 2),
            100 * mean(W$dq_g > 1e-9), 100 * mean(W$dq_g < -1e-9), max(W$dq_g), min(W$dq_g),
            quantile(G$cost, .9), quantile(G$cost, .99)))
fwrite(pp, file.path(out_dir, "per_person_facts.csv"))
fwrite(W[, .(person, risk_group, dq_g, dc_g, dq_k, dc_k)], file.path(out_dir, "per_person_increments.csv.gz"))

# ---- per-patient files ----------------------------------------------------------
set.seed(seed); keep <- sample(unique(X$person), min(10000, n_people))
fwrite(X[person %in% keep][order(person, arm)], file.path(out_dir, "patients_sample.csv"))
fwrite(X, file.path(out_dir, "patients.csv.gz"))

# ---- figure ----------------------------------------------------------------------
png(file.path(out_dir, "microsim_distributions.png"), 2400, 1400, res = 200)
op <- par(mfrow = c(1, 3), mar = c(5, 4.5, 3.5, 1), family = "sans")
hist(pmin(lead, 8), breaks = seq(0, 8, 0.25), col = "#199e70", border = "white", main = "Lead time per screen-detected cancer",
     xlab = "Years earlier than in usual care (same person)", ylab = "People")
b <- dist[measure == "Positive Galleri tests per person"][-1]
bp <- barplot(1e5 * b$share, names.arg = b$value, col = "#d95926", border = NA, ylim = c(0, 1.15e5 * max(b$share)),
              main = "People with 1, 2, 3, 4+ positive tests", ylab = "People per 100,000", xlab = "Positives over a lifetime of annual tests")
text(bp, 1e5 * b$share, format(round(1e5 * b$share, 1), big.mark = ","), pos = 3, cex = 0.9)
dq <- W$dq_g[abs(W$dq_g) > 1e-9]
hist(pmax(pmin(dq, 5), -1), breaks = 60, col = "#9085e9", border = "white", main = "QALY change per person affected",
     xlab = "Galleri minus usual care, discounted QALYs (clipped -1 to 5)", ylab = "People")
par(op); dev.off()
cat(sprintf("done: %s people, %s rows\n", format(n_people, big.mark = ","), format(nrow(X), big.mark = ",")))
print(V[arm == "Galleri", .(measure, microsim = round(microsim, 3), mc_95 = round(mc_95, 3), cohort = round(cohort, 3), diff_pct = round(diff_pct, 2))])
print(ICR)
