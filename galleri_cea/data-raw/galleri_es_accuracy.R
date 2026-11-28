# Galleri esophagus/stomach (E/S) signal accuracy from PATHFINDER 2 (Nabavizadeh 2026, Nat Med)
# and NHS-Galleri (Sasieni 2026, NEJM). Reads data-raw/raw/galleri_es_accuracy.csv, prints the
# derived inputs with Wilson 95% CIs, and writes them to inputs/parameters.csv.
Sys.setlocale("LC_CTYPE", "C.UTF-8")
A <- read.csv("data-raw/raw/galleri_es_accuracy.csv", check.names = FALSE)
wilson <- function(x, n, z = 1.96) { p <- x / n; c0 <- p + z^2 / (2 * n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(est = p, lo = (c0 - h) / (1 + z^2 / n), hi = (c0 + h) / (1 + z^2 / n)) }
g <- function(m) { r <- A[A$measure == m, ]; c(x = r$numerator, n = r$denominator) }
show <- function(m) { v <- g(m); w <- wilson(v[["x"]], v[["n"]])
  cat(sprintf("%-62s %5d/%-6d %.4f%% (%.4f-%.4f%%)\n", m, v[["x"]], v[["n"]], 100 * w[1], 100 * w[2], 100 * w[3])); invisible(w) }
cat("PATHFINDER 2, 12-month follow-up\n")
s_ugi <- show("Upper GI cancers with a positive E/S signal (sens_g_ugi)")
show("Upper GI cancers detected by the MCED test (true positive)")
s_oth <- show("Other-site cancers given an E/S CSO1 (sens_g_other)")
fpr <- show("E/S false-positive rate (fpr_g_es)")
ppv_u <- show("E/S signal PPV for upper GI cancer")
show("E/S signal PPV for any cancer"); show("False positives with CSO1 = stomach/esophagus")

# sens_g_ugi: only 9 upper GI cancers, so the trial count updates the previous indirect estimate
# (0.637 = 69.8% episode sensitivity for 12 deadly cancers x 91.3% CSO accuracy, weight n = 60)
# by conjugate beta updating: (0.637*60 + 4) / (60 + 9).
prior <- c(m = 0.637, n = 60); v <- g("Upper GI cancers with a positive E/S signal (sens_g_ugi)")
sens_ugi <- (prior[["m"]] * prior[["n"]] + v[["x"]]) / (prior[["n"]] + v[["n"]])
cat(sprintf("\nsens_g_ugi: prior 0.637 (n 60) + PATHFINDER 2 4/9 -> %.3f (n %d)\n", sens_ugi, prior[["n"]] + v[["n"]]))
# Repeat rounds (NHS-Galleri): positivity and programme sensitivity, rounds 2-3 vs round 1
pos <- A$value[grepl("^Positive tests round", A$measure)]; ps <- A$value[grepl("^Programme sensitivity", A$measure)]
cat(sprintf("NHS-Galleri positivity ratio rounds 2-3 / 1: %.3f; programme sensitivity ratio: %.3f\n",
            mean(pos[2:3]) / pos[1], mean(ps[2:3]) / ps[1]))
# Upper GI stage, NHS-Galleri after 3 rounds (all diagnoses in each arm, not only MCED-detected)
st <- function(s) as.numeric(strsplit(s, " / ")[[1]])
ugi <- function(arm) st(A$note[A$measure == sprintf("Esophagus cancers %s (I/II/III/IV/missing)", arm)]) +
                     st(A$note[A$measure == sprintf("Stomach cancers %s (I/II/III/IV/missing)", arm)])
I <- ugi("intervention"); C <- ugi("control")
cat(sprintf("NHS-Galleri upper GI: stage IV %d vs %d (ratio %.2f); stage I-II %d vs %d; staged late share %.1f%% vs %.1f%%\n",
            I[4], C[4], I[4] / C[4], I[1] + I[2], C[1] + C[2], 100 * sum(I[3:4]) / sum(I[1:4]), 100 * sum(C[3:4]) / sum(C[1:4])))

# ---- PMA Galleri (GRAIL FDA panel executive summary, Sept 2026) -------------------
# The marketed test (PMA "Galleri") was rerun on PATHFINDER 2 and NHS-Galleri samples; the
# published papers used the earlier MCED-V2 classifier. The model uses the PMA test.
val <- function(m) A$value[A$measure == m]
ep_eso <- val("Esophagus episode sensitivity combined (random effects)")      # 0.642 (Table 33)
ep_stom <- val("Stomach episode sensitivity combined (random effects)")       # 0.466
stg <- function(m) as.numeric(strsplit(A$note[A$measure == m], " / ")[[1]])
se <- stg("Galleri-detected esophagus cancers by stage (I/II/III/IV/NA)"); ss <- stg("Galleri-detected stomach cancers by stage (I/II/III/IV/NA)")
early_n <- se[1] + se[2] + ss[1] + ss[2]; staged_n <- sum(se[1:4]) + sum(ss[1:4])
cso <- g("CSO accuracy PATHFINDER 2") + g("CSO accuracy NHS-Galleri")
cso_acc <- cso[["x"]] / cso[["n"]]
fp_all <- g("PATHFINDER 2 false-positive rate (all CSOs)") + g("NHS-Galleri round 1 false-positive rate (all CSOs)")
fpr_all <- fp_all[["x"]] / fp_all[["n"]]
es_share <- g("False positives with CSO1 = stomach/esophagus")                 # 4/114 (MCED-V2; not reported for PMA)
fpr_es <- fpr_all * es_share[["x"]] / es_share[["n"]]
cso_v2 <- val("CSO1 accuracy")                                                 # MCED-V2, PATHFINDER 2
oth <- g("Other-site cancers given an E/S CSO1 (sens_g_other)")
sens_oth <- oth[["x"]] / oth[["n"]] * (1 - cso_acc) / (1 - cso_v2)              # wrong-site E/S calls scale with CSO error
rt_pos <- g("Retest after positive test and negative workup: positive retest"); rt_neg <- g("Retest after positive test and negative workup: negative retest")
p_again <- (rt_pos[["n"]] - rt_pos[["x"]]) / (rt_pos[["n"]] - rt_pos[["x"]] + rt_neg[["n"]] - rt_neg[["x"]])   # cancer-free who retest positive
kappa <- (1 - p_again) / (p_again - fpr_all)            # beta-binomial: P(+ again | + once, no cancer) = (m kappa + 1) / (kappa + 1)
cat(sprintf("\nPMA Galleri: esophagus episode sens %.3f, stomach %.3f, stage I-II among detected %d/%d = %.3f\n", ep_eso, ep_stom, early_n, staged_n, early_n / staged_n))
cat(sprintf("CSO accuracy %d/%d = %.3f; FPR all CSOs %d/%d = %.4f%%; E/S FPR %.5f%%; other-site E/S calls %.4f\n",
            cso[["x"]], cso[["n"]], cso_acc, fp_all[["x"]], fp_all[["n"]], 100 * fpr_all, 100 * fpr_es, sens_oth))
cat(sprintf("Retest (Westgate): cancer-free positive again %d/%d = %.3f -> fp_kappa %.2f\n",
            rt_pos[["n"]] - rt_pos[["x"]], rt_pos[["n"]] - rt_pos[["x"]] + rt_neg[["n"]] - rt_neg[["x"]], p_again, kappa))

# ---- write to inputs/parameters.csv ------------------------------------------------
P <- read.csv("inputs/parameters.csv", colClasses = "character", check.names = FALSE)
f <- function(x, d = 3) formatC(x, digits = d, format = "fg")
row <- function(parameter, label, value, status, dist, n, low, high, unit, description, source, refs)
  data.frame(parameter, label, value, status, dist, se = "", n, group = "", low, high, unit, description, source,
             analysis = "both", refs, price_value = "", price_year = "", check.names = FALSE)
new <- rbind(
  row("sens_ep_eso", "Galleri episode sensitivity, esophagus", f(ep_eso), "SOURCED", "beta", "28", "0.453", "0.796", "proportion",
      "12-month episode sensitivity for esophageal cancer (PMA Galleri). The model solves stage-specific sensitivity from this (galleri_stage_sens)",
      "Random-effects pooled PATHFINDER 2 (1/2) and NHS-Galleri round 1 (17/26): 64.2% (45.3-79.6%) (GRAIL FDA summary 2026, Table 33)", "62"),
  row("sens_ep_stom", "Galleri episode sensitivity, stomach", f(ep_stom), "SOURCED", "beta", "11", "0.200", "0.753", "proportion",
      "12-month episode sensitivity for gastric cancer (PMA Galleri)",
      "Pooled PATHFINDER 2 (1/4) and NHS-Galleri round 1 (4/7): 46.6% (20.0-75.3%) (GRAIL FDA summary 2026, Table 33)", "62"),
  row("screen_early_share", "Galleri-detected upper GI cancers at stage I-II", f(early_n / staged_n), "SOURCED", "beta", as.character(staged_n), "0.07", "0.39", "proportion",
      "Share of Galleri-detected esophageal and gastric cancers at stage I-II; sets the early / late sensitivity ratio",
      sprintf("Esophagus 2 of 17 staged (I 0, II 2, III 10, IV 5), stomach 2 of 5 (II 2, III 2, IV 1): %d/%d (GRAIL FDA summary 2026, Figure 12). Range: Wilson 95%% CI", early_n, staged_n), "62"),
  row("cso_acc", "Galleri CSO accuracy (true positives)", f(cso_acc), "SOURCED", "beta", as.character(cso[["n"]]), "0.89", "0.95", "proportion",
      "Share of true positives given the right cancer signal origin; an upper GI cancer is sent to EGD only with a stomach/esophagus CSO",
      sprintf("PATHFINDER 2 100/106 and NHS-Galleri 238/259 = %d/%d (GRAIL FDA summary 2026, Tables 19 and 31)", cso[["x"]], cso[["n"]]), "62"))
P <- rbind(P[!P$parameter %in% c(new$parameter, "sens_g_ugi"), ], new)
put <- function(k, value, status, n, low, high, source, refs) {
  i <- P$parameter == k; P$value[i] <<- value; P$status[i] <<- status; P$n[i] <<- n
  P$low[i] <<- low; P$high[i] <<- high; P$source[i] <<- source; P$refs[i] <<- refs }
put("fpr_g_es", f(fpr_es), "SOURCED", "49000", f(fpr_all * 0.0137, 2), f(fpr_all * 0.0868, 2),
    sprintf("PMA Galleri false-positive rate, all CSOs, PATHFINDER 2 + NHS-Galleri: %d/%d = %.3f%% (GRAIL FDA summary 2026, Tables 13, 25) x share of false positives with a stomach/esophagus CSO 4/114 (PATHFINDER 2, MCED-V2; Suppl. Table S7). Range: Wilson CI of 4/114", fp_all[["x"]], fp_all[["n"]], 100 * fpr_all),
    "1, 62")
put("sens_g_other", f(sens_oth), "SOURCED", "740", f(sens_oth * 0.39, 2), f(sens_oth * 2.54, 2),
    sprintf("PATHFINDER 2 (MCED-V2): 4 of 431 other-site cancers got a stomach/esophagus CSO (Ext. Data Fig. 4), scaled by the PMA test's lower CSO error ((1 - %.3f) / (1 - %.3f)). Range: Wilson CI of 4/431", cso_acc, cso_v2),
    "1, 62")
put("fp_kappa", f(kappa), "SOURCED", "", "", "",
    sprintf("Person-level persistence of false positives (beta-binomial concentration). Real-world retests after a positive test and a negative workup: %d of %d people without cancer tested positive again (Westgate et al., cited in GRAIL FDA summary 2026): (1 - %.3f) / (%.3f - %.4f) = %.2f. Retest within 6 months, all CSOs; 1e6 = no persistence (scenario)",
            rt_pos[["n"]] - rt_pos[["x"]], rt_pos[["n"]] - rt_pos[["x"]] + rt_neg[["n"]] - rt_neg[["x"]], p_again, p_again, fpr_all, kappa), "62")
P$dist[P$parameter == "fp_kappa"] <- "fixed"
i <- P$parameter == "p_late_screen"
P$label[i] <- "P(late | screen-detected), other-site cancers"
P$description[i] <- "Share of other-site cancers found during the upper GI workup that are late stage. Upper GI cancers take their stage at detection from the natural history and stage-specific Galleri sensitivity"
i <- P$parameter == "sens_g_rep_ratio"
P$source[i] <- paste0(sub(" Also: .*", "", P$source[i]), sprintf(". Also: NHS-Galleri positivity 1.03%% / 0.80%% / 0.90%% by round (ratio %.2f) and programme sensitivity 43.1%% / 32.2%% / 27.0%% (ratio %.2f; Sasieni 2026)", mean(pos[2:3]) / pos[1], mean(ps[2:3]) / ps[1]))
P$refs[i] <- "2, 61"
i <- P$parameter == "p_egd_pos"
P$value[i] <- f(ppv_u[["est"]] * 0.887, 2); P$n[i] <- "12"
P$source[i] <- "PATHFINDER 2: 4 of 12 E/S-signal positives had an upper GI cancer (8 of 12 had a cancer somewhere; Ext. Data Fig. 4, Suppl. Table S7) x EGD sensitivity 88.7% (Menon 2014)"
P$refs[i] <- "1, 6"
P[P == ""] <- NA
data.table::fwrite(P, "inputs/parameters.csv", na = "")
R <- read.csv("inputs/references.csv", colClasses = "character")
add <- data.frame(ref = c("61", "62"),
  citation = c("Sasieni P, Johnson P, Round T, et al. Effect of screening with multicancer early-detection test on late-stage cancer diagnosis. N Engl J Med. 2026;395(14):1374-1386.",
               "GRAIL. Galleri multi-cancer early detection test: executive summary for the FDA Molecular and Clinical Genetics Panel, September 23, 2026 (incl. clinician follow-up guide citing Westgate C, et al., early real-world experience with positive MCED test cases and negative initial diagnostic work-up, poster)."),
  url = c("https://doi.org/10.1056/NEJMoa2505723", "https://www.fda.gov/advisory-committees"), pmid = NA)
R <- rbind(R[!R$ref %in% add$ref, ], add); R[R == ""] <- NA
data.table::fwrite(R, "inputs/references.csv", na = "")

# ---- scenarios ------------------------------------------------------------------------
S <- read.csv("inputs/scenarios.csv", colClasses = "character", check.names = FALSE)
S <- S[!S$parameter %in% "sens_g_ugi" & !grepl("^Galleri E/S accuracy as before", S$scenario), ]
add <- data.frame(
  scenario = c("Galleri sensitivity at lower 95% CI (esophagus 45%, stomach 20%)", "Galleri sensitivity at lower 95% CI (esophagus 45%, stomach 20%)",
               "Galleri sensitivity at upper 95% CI (esophagus 80%, stomach 75%)", "Galleri sensitivity at upper 95% CI (esophagus 80%, stomach 75%)",
               "Galleri finds more early cancers (stage I-II share 39%, upper CI)", "Galleri finds fewer early cancers (stage I-II share 7%, lower CI)",
               "Galleri sensitivity equal at every stage (no stage gradient)",
               "False positives independent between tests (no persistence)",
               "Earlier test version (MCED-V2): E/S false-positive rate 0.0127%", "Earlier test version (MCED-V2): E/S false-positive rate 0.0127%",
               "Repeat-round sensitivity at NHS-Galleri programme-sensitivity ratio (0.69)"),
  parameter = c("sens_ep_eso", "sens_ep_stom", "sens_ep_eso", "sens_ep_stom", "screen_early_share", "screen_early_share",
                "screen_early_share", "fp_kappa", "fpr_g_es", "sens_g_other", "sens_g_rep_ratio"),
  value = c("0.453", "0.200", "0.796", "0.753", "0.39", "0.07", "-1", "1e6", "0.000127", "0.00928", f(mean(ps[2:3]) / ps[1], 2)),
  refs = c("62", "62", "62", "62", "62", "62", NA, NA, "1", "1", "61"))
S <- rbind(S[!S$scenario %in% add$scenario, ], add); S[S == ""] <- NA
data.table::fwrite(S, "inputs/scenarios.csv", na = "")
