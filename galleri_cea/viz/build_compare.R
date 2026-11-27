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
