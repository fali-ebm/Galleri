# Build viz/galleri_vs_usual_care.html from viz/galleri_viz_data.json.
#   Rscript viz/build_compare.R       (run from galleri_cea/, after Rscript tools/viz_data.R)
# Shares the burden page's style: the <style> block of burden_template.html is copied in.
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))   # UTF-8 text (dashes, arrows) whatever the system locale
suppressPackageStartupMessages(library(jsonlite))
here <- "viz"
D <- fromJSON(file.path(here, "galleri_viz_data.json"), simplifyVector = FALSE)
R <- list(); for (r in D$comparison) R[[paste(r$key, r$arm, sep = "|")]] <- r
v <- function(k, a = "Galleri") R[[paste(k, a, sep = "|")]]$value
UC <- D$by_age$`Usual care`; G <- D$by_age$Galleri
sc <- list(); for (r in D$scenarios) if (r$arm == "Galleri") sc[[r$scenario]] <- r
base <- sc[["Base case"]]$icer; half <- sc[["Galleri test half price"]]$icer
# The ICER is linear in the test price: two points (list price, half price) give the line.
slope <- (base - half) / (D$meta$c_galleri / 2); icpt <- base - slope * D$meta$c_galleri
d_cost <- v("cost_total") - v("cost_total", "Usual care")                 # per person, discounted
ugi_uc <- sum(vapply(UC, function(r) r$ca_death, 0)); ugi_g <- sum(vapply(G, function(r) r$ca_death, 0))
K <- "Kunzmann -> EGD"
data <- list(
  per = D$meta$per, price = D$meta$c_galleri, cost_year = D$meta$cost_year, wtp = D$meta$wtp,
  icer = base, icer_lo = R[["icer|Galleri"]]$lower_95, icer_hi = R[["icer|Galleri"]]$upper_95,
  d_cost_total = d_cost * 1e5, qaly_net = v("qaly_net"), ly_gained = v("ly_gained"), late_averted = v("late_averted"),
  deaths_averted = v("ca_deaths_averted"), ugi_deaths_uc = ugi_uc, ugi_deaths_g = ugi_g,
  other_deaths_averted = v("ca_deaths_averted") - (ugi_uc - ugi_g), tests = v("galleri_tests"),
  overdx = v("overdiagnosed"), false_alarms = v("free_investigated"),
  price_icpt = icpt, price_slope = slope,
  cum = Map(function(a, b) list(age = a$age, uc = b$ca_death, g = a$ca_death), G, UC),
  psa = D$psa,
  scenarios = unname(lapply(names(sc), function(k) list(s = k, icer = sc[[k]]$icer, dq = sc[[k]]$inc_effect))),
  subgroups = lapply(D$subgroups, function(r) list(age = r$start_age, sex = r$sex, race = r$race_eth, label = r$race_label,
    icer = r$icer_vs_usual_care, dq = r$inc_qaly, averted = r$cancer_deaths_averted)),
  kunz = list(egd = v("egd", K), perf_deaths = v("perforation_deaths", K), qaly_net = v("qaly_net", K),
              averted = v("ca_deaths_averted", K), cost = (v("cost_total", K) - v("cost_total", "Usual care")) * 1e5))

# ---- 07 microsimulation: one person at a time (outputs/microsim, Rscript tools/run_microsim.R) ----
csv <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
MS <- csv("outputs/microsim/summary.csv"); MV <- csv("outputs/microsim/validation.csv")
MF <- csv("outputs/microsim/per_person_facts.csv"); mf <- function(m) as.numeric(MF$value[MF$measure == m])
GV <- MS[MS$comparison == "Galleri vs usual care", ]
inc <- readr::read_csv("outputs/microsim/per_person_increments.csv.gz", col_select = "dq_g", show_col_types = FALSE, progress = FALSE)$dq_g
# symmetric, log-spaced bins of each person's QALY change (people with any change only)
edges <- c(0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)
bin_of <- function(v) findInterval(abs(v), edges, rightmost.closed = TRUE)
hist_side <- function(v) { b <- tabulate(bin_of(v), nbins = length(edges) - 1); data.frame(lo = head(edges, -1), hi = edges[-1], n = b) }
ch <- inc[inc != 0]
mv <- function(m) { r <- MV[MV$arm == "Galleri" & MV$measure == m, ]; list(cohort = r$cohort, microsim = r$microsim) }
ms <- list(
  n = length(inc), gain = GV$share_gain, loss = GV$share_loss, same = GV$share_same,
  icer = GV$icer, icer_cohort = GV$cohort_icer, qaly = GV$qaly, qaly95 = GV$qaly_mc95, cost = GV$cost,
  gains = hist_side(ch[ch > 0]), losses = hist_side(ch[ch < 0]), per = 1e5 / length(inc),
  max_gain = max(inc), max_loss = min(inc),
  lead_med = mf("Lead time, median years (screen-detected, presented in usual care)"),
  shift = mf("Late stage in usual care -> early when found by Galleri, % of screen-detected"),
  found_share = mf("Galleri finds it, % of those cancers"), ever_pos = mf("People with any positive test, %"))
ms$compare <- unname(Map(function(k, l) c(list(label = l), mv(k)), c("screen_detected", "clin_dx", "overdiagnosed", "ca_death", "false_pos"),
  c("Found by screening", "Diagnosed after symptoms", "Overdiagnosed", "Upper GI cancer deaths", "Healthy men worked up")))

# ---- 08 MISCAN-ESO calibration of esophageal adenocarcinoma (outputs/miscan_eac) ----
MR <- csv("outputs/miscan_eac/rates.csv"); MT <- csv("outputs/miscan_eac/targets.csv"); MB <- csv("outputs/miscan_eac/by_age.csv")
ext <- MT[MT$kind == "external validation" & !grepl("per EAC case", MT$target), ]
miscan <- list(
  rates = lapply(split(MR, MR$sex), function(r) list(S = r$sojourn_target_y, L = r$late_share_target, pe = r$mean_years_in_PE,
    pl = r$mean_years_in_PL, hE = r$rate_PE_to_symptoms_early, lam = r$rate_PE_to_PL, hL = r$rate_PL_to_symptoms_late)),
  by_age = lapply(split(MB, MB$sex), function(d) list(age = d$age, seer = 1e5 * d$seer_incidence, model = 1e5 * d$model_incidence)),
  validation = lapply(seq_len(nrow(ext)), function(i) list(
    label = paste(if (grepl("cases", ext$target[i])) "Cases" else "Deaths", if (grepl("women", ext$target[i])) "women" else "men", sep = ", "),
    cisnet = ext$target_value[i], model = ext$model_value[i])),
  fit = lapply(which(MT$kind == "calibration"), function(i) list(target = MT$target[i], value = MT$target_value[i], model = MT$model_value[i])),
  scenarios = unname(lapply(names(sc)[grepl("^MISCAN-ESO|^Base case$", names(sc))], function(k) list(s = k, icer = sc[[k]]$icer))))
data$ms <- ms; data$miscan <- miscan

rd <- function(f) paste(readLines(file.path(here, f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
burden <- rd("burden_template.html")
style <- regmatches(burden, regexpr("(?s)<style>.*?</style>", burden, perl = TRUE))
style <- sub("^<style>", "", sub("</style>$", "", style))
js <- as.character(toJSON(data, auto_unbox = TRUE, digits = NA, null = "null", na = "null"))
html <- rd("compare_template.html")
html <- sub("__CSS__", style, html, fixed = TRUE); html <- sub("__DATA__", js, html, fixed = TRUE)
writeLines(html, file.path(here, "galleri_vs_usual_care.html"), useBytes = TRUE)
cat(sprintf("wrote viz/galleri_vs_usual_care.html %d; break-even price at $100k/QALY $%.2f, at $150k $%.2f\n",
            nchar(html), (1e5 - icpt) / slope, (1.5e5 - icpt) / slope))
