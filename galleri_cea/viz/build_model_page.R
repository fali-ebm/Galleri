# Build viz/galleri_model.html: the decision model, step by step, from the model's own files.
#   Rscript viz/build_model_page.R    (run from galleri_cea/, after run.R and tools/run_microsim.R)
#
# Reads inputs/parameters.csv, inputs/settings.csv, outputs/galleri, outputs/comparison and
# outputs/microsim. Every tree edge label is computed from parameters.csv, so the page cannot
# drift from the inputs the model ran with.
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))   # UTF-8 text (dashes, arrows) whatever the system locale

rd <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
PR <- utils::read.csv("inputs/parameters.csv", colClasses = "character", check.names = FALSE)   # values shown as written
P <- setNames(as.list(as.numeric(PR$value)), PR$parameter)
PST <- split(PR, PR$parameter)
SR <- rd("inputs/settings.csv"); S <- setNames(as.list(SR$value), SR$setting)
GL <- rd("outputs/galleri/results_long.csv"); CL <- rd("outputs/comparison/results_long.csv")
SC <- rd("outputs/galleri/scenarios.csv"); SC <- SC[SC$arm == "Galleri", ]
MS <- rd("outputs/microsim/summary.csv"); MV <- rd("outputs/microsim/validation.csv")
MF <- rd("outputs/microsim/per_person_facts.csv"); MD <- rd("outputs/microsim/distributions.csv")

getv <- function(D, arm, key, col = "value") { v <- D[D$arm == arm & D$key == key, col]; if (!length(v)) NA_real_ else as.numeric(v[1]) }
g  <- function(k, a = "Galleri") getv(GL, a, k)
gl <- function(k) c(getv(GL, "Galleri", k, "lower_95"), getv(GL, "Galleri", k, "upper_95"))
KZ <- "Kunzmann -> EGD"
kz <- function(k) getv(CL, KZ, k)
ms <- function(cmp, col) as.numeric(MS[MS$comparison == cmp, col])
mf <- function(m) MF$value[MF$measure == m]
esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE); x <- gsub("\"", "&quot;", x, fixed = TRUE); gsub("'", "&#x27;", x, fixed = TRUE) }

# ---- number formats ----------------------------------------------------------
MINUS <- "−"
fx <- function(x, d) formatC(x, format = "f", digits = d, big.mark = ",")
money <- function(x) {
  a <- abs(x)
  s <- if (a >= 1e6) sprintf("$%.2fM", a / 1e6) else if (a >= 1e4) sprintf("$%.0fk", a / 1e3) else paste0("$", fx(a, 0))
  paste0(if (x < 0) MINUS else "", s)
}
n <- function(x, d = 0) paste0(if (x < 0) MINUS else "", fx(abs(x), d))
pct <- function(x, d = 1) sprintf("%.*f%%", d, 100 * x)
pr <- function(x) {                       # a probability as a short edge label
  if (x >= 0.995 || x == 0) return(sprintf("%.0f", x))
  if (x >= 0.1) { s <- sprintf("%.2f", x); return(if (endsWith(s, "0")) sub("\\.$", "", sub("0+$", "", s)) else s) }
  if (x >= 0.01) return(sub("0+$", "", sprintf("%.3f", x)))
  sub("0+$", "", sprintf("%.4f", x))
}
num <- function(x) format(x, digits = 15)  # plain numbers in prose, as stored

# ---- left-to-right tree layout -----------------------------------------------
# A node: text t, subtitle sub, class cls (root/arm/dec/found/harm/out/dash), edge label el
# (on the edge coming in), children kids, hot = rose edge, dash = dashed edge, w = minimum width.
N <- function(t, sub = "", cls = "dec", el = "", kids = list(), hot = FALSE, dash = FALSE, w = 0)
  list(t = t, sub = sub, cls = cls, el = el, kids = kids, hot = hot, dash = dash, w = w)
tw <- function(t, bold) nchar(t) * (if (bold) 7.2 else 6.3)   # rough text width at 13px bold / 11.5px

tree_svg <- function(root, label, gap = 78, rowh = 58, top = 14, left = 8) {
  # flatten (pre-order) with depth, parent and leaf rows; a parent sits midway between its first and last child
  nodes <- list(); nleaf <- 0
  visit <- function(z, d, parent) {
    i <- length(nodes) + 1; nodes[[i]] <<- c(z[c("t", "sub", "cls", "el", "hot", "dash", "w")], d = d, parent = parent)
    kid_ids <- integer(0)
    for (k in z$kids) kid_ids <- c(kid_ids, visit(k, d + 1, i))
    if (!length(kid_ids)) { nodes[[i]]$row <<- nleaf; nleaf <<- nleaf + 1 }
    else nodes[[i]]$row <<- (nodes[[kid_ids[1]]]$row + nodes[[kid_ids[length(kid_ids)]]]$row) / 2
    nodes[[i]]$kids <<- kid_ids
    i
  }
  visit(root, 0, 0)
  d <- vapply(nodes, `[[`, 0, "d"); maxd <- max(d)
  bw <- vapply(nodes, function(z) max(z$w, 120, tw(z$t, TRUE) + 22, tw(z$sub, FALSE) + 22), 0)
  colw <- vapply(0:maxd, function(k) max(bw[d == k]), 0)
  xs <- cumsum(c(left, colw[-length(colw)] + gap))
  W <- as.integer(xs[maxd + 1] + colw[maxd + 1] + 12); H <- as.integer(top + nleaf * rowh + 4)
  h <- ifelse(vapply(nodes, function(z) nzchar(z$sub), TRUE), 42, 32)
  bx <- xs[d + 1]; by <- top + vapply(nodes, `[[`, 0, "row") * rowh + (42 - h) / 2
  f0 <- function(v) sprintf("%.0f", v)
  edges <- character(0); boxes <- character(0); labels <- character(0)
  for (i in seq_along(nodes)) {
    z <- nodes[[i]]; x1 <- bx[i] + bw[i]; y1 <- by[i] + h[i] / 2
    for (k in z$kids) {
      mx <- xs[d[i] + 1] + colw[d[i] + 1] + 14; y2 <- by[k] + h[k] / 2; x2 <- bx[k] - 3
      cls <- paste0("e", if (nodes[[k]]$hot) " hot" else "", if (nodes[[k]]$dash) " dash" else "")
      edges <- c(edges, sprintf('<path class="e" d="M%s,%s H%s V%s" fill="none"/>', f0(x1), f0(y1), f0(mx), f0(y2)),
                 sprintf('<path class="%s" d="M%s,%s H%s" fill="none" marker-end="url(#ah)"/>', cls, f0(mx), f0(y2), f0(x2)))
      if (nzchar(nodes[[k]]$el))
        labels <- c(labels, sprintf('<text class="el" x="%s" y="%s" text-anchor="middle">%s</text>',
                                    f0((mx + x2) / 2 + 2), f0(y2 - 5), esc(nodes[[k]]$el)))
    }
    c0 <- z$cls; cx <- f0(bx[i] + bw[i] / 2)
    boxes <- c(boxes, sprintf('<rect class="n %s" x="%s" y="%s" width="%s" height="%d" rx="7"/>', c0, f0(bx[i]), f0(by[i]), f0(bw[i]), h[i]))
    if (nzchar(z$sub)) boxes <- c(boxes,
      sprintf('<text class="nt b %s" x="%s" y="%s" text-anchor="middle">%s</text>', c0, cx, f0(by[i] + 17), esc(z$t)),
      sprintf('<text class="nt s %s" x="%s" y="%s" text-anchor="middle">%s</text>', c0, cx, f0(by[i] + 33), esc(z$sub)))
    else boxes <- c(boxes, sprintf('<text class="nt b %s" x="%s" y="%s" text-anchor="middle">%s</text>', c0, cx, f0(by[i] + 21), esc(z$t)))
  }
  svg <- c(sprintf('<svg viewBox="0 0 %d %d" role="img" aria-label="%s" xmlns="http://www.w3.org/2000/svg">', W, H, esc(label)),
           '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>',
           edges, boxes, labels, "</svg>")
  list(svg = paste(svg, collapse = "\n"), W = W)
}

# ---- the trees (labels from parameters) --------------------------------------
hold <- as.integer(P$egd_hold_years); se <- P$sens_egd; fi <- P$f_invisible
kshare <- as.numeric(S$kunz_share_male); kcap <- as.numeric(S$kunz_capture_male)

arms <- N("Cohort eligible", "50-year-old NH White men", "root", kids = list(
  N("1 · Usual care", "", "arm", kids = list(N("Cancer presents with symptoms", "SEER incidence and stage", "out", kids = list(
    N("Early stage", "", "out"), N("Late stage", "", "out"))))),
  N("2 · Kunzmann score, once", "at the start age", "arm", kids = list(
    N("Score ≥ 8", sprintf("%s flagged; sens. %s", pct(kshare), pct(kcap)), "dec", el = pr(kshare), kids = list(
      N("One EGD", "start age only", "dec", kids = list(
        N("Cancer found", "", "found", el = paste("sens", pr(se))),
        N("No cancer found", "→ usual care", "out"))))),
    N("Score < 8", "→ usual care", "out", el = pr(1 - kshare)))),
  N("3 · Galleri, every year", paste("to age", S$screen_stop_age), "arm", kids = list(
    N("Galleri +ve (E/S signal)", "→ figure 2", "dec", kids = list(
      N("Cancer present", sprintf("sens %s, then ×%s", pct(P$sens_g_ugi), num(P$sens_g_rep_ratio)), "found"),
      N("No upper GI cancer", sprintf("FP rate %s per test", pct(P$fpr_g_es, 3)), "harm", hot = TRUE))),
    N("Galleri −ve", "", "out", kids = list(
      N("True negative", "retested next year", "out"),
      N("False negative", "symptoms later or next round", "harm", hot = TRUE)))))))

egd <- N("Galleri +ve", "E/S signal", "root", kids = list(
  N(sprintf("Negative EGD in last %d y?", hold), "the EGD hold", "dec", kids = list(
    N("Yes: EGD skipped", "go straight to figure 3", "dash"),
    N("No: EGD", sprintf("perforation %s, %s fatal", pct(P$p_perf_egd, 2), pct(P$p_death_perf, 1)), "dec", kids = list(
      N("Cancer present", "", "dec", kids = list(
        N("Found at EGD", "", "found", el = pr(se)),
        N("Missed: visible lesion", "→ figure 3; stays preclinical", "harm", el = pr((1 - se) * (1 - fi)), hot = TRUE),
        N("Missed: not yet visible", sprintf("visible after %s y", num(P$t_visible)), "harm", el = pr((1 - se) * fi), hot = TRUE))),
      N("No upper GI cancer", "EGD −ve", "out", kids = list(
        N(sprintf("Hold starts (%d y)", hold), "→ figure 3", "out")))))))))

pet_ugi <- P$p_pet_ugi
rest <- N("EGD −ve or skipped", "", "root", kids = list(
  N("ENT evaluation", "", "dec", el = pr(P$p_ent_route), kids = list(
    N("ENT +ve", "cancer found", "found", el = pr(P$p_ent_pos)),
    N("ENT −ve", paste("missed cancer", pct(P$p_miss)), "out", el = pr(1 - P$p_ent_pos)))),
  N("PET CT", "", "dec", el = pr(1 - P$p_ent_route), kids = list(
    N("Pancreatobiliary +ve", "→ figure 4", "dec", el = pr(P$p_pet_pb)),
    N("Incidental finding", sprintf("cancer %s; rest benign", pct(P$p_incid_cancer, 0)), "harm", el = pr(P$p_pet_incid)),
    N("Upper GI cancer seen", "targeted EGD; if one is visible", "found", el = paste0(pr(pet_ugi), "*"), dash = TRUE),
    N("PET −ve", paste("missed cancer", pct(P$p_miss)), "out", el = pr(1 - P$p_pet_pb - P$p_pet_incid))))))

ps <- P$p_eus_samp; M <- P$p_pb_malig; est <- P$eus_tissue_sens; tl <- P$p_tb_lap
tb <- function() N("Tumor board", "", "dec", kids = list(
  N("Diagnostic laparoscopy", "", "harm", el = pr(tl)), N("Watchful waiting", "surveillance imaging", "out", el = pr(1 - tl))))
pl2 <- tl * P$r_lap_after_eus; pe2 <- P$p_reeus_after_eus
pl3 <- tl * P$r_lap_after_lap; pe3 <- P$p_reeus_after_lap
pb1 <- N("First finding: EUS", paste("malignant", pct(M, 0)), "root", kids = list(
  N("FNA/FNB taken", "coin toss", "dec", el = pr(ps), kids = list(
    N("Tissue +ve", "cancer found", "found", el = paste0(pr(M), "×", pr(est))),
    N("Tissue −ve", "", "out", el = pr(1 - M * est), kids = list(tb())))),
  N("No sample", "", "out", el = pr(1 - ps), kids = list(tb()))))
rep_branch <- function(lbl, pl, pe, r) N(lbl, "", "dec", kids = list(
  N("Laparoscopy", sprintf("%s × %s", pr(tl), pr(r)), "harm", el = pr(pl)),
  N("Repeat EUS", "", "dec", el = pr(pe)),
  N("Watchful waiting", "", "out", el = pr(1 - pl - pe))))
pb2 <- N("Repeat finding", "personalised decision", "root", kids = list(
  rep_branch("After a negative EUS", pl2, pe2, P$r_lap_after_eus),
  rep_branch("After a negative laparoscopy", pl3, pe3, P$r_lap_after_lap)))

nh <- N("Well", "no upper GI cancer", "out", kids = list(
  N("Preclinical cancer", sprintf("mean sojourn %s y", num(P$sojourn_years)), "dec", el = "onset", kids = list(
    N("Screen-detected", sprintf("late %s; lead time", pct(P$p_late_screen, 0)), "found", kids = list(
      N("Would-be symptom time", "stage-specific survival", "out"))),
    N("Not yet visible at EGD", "after a +ve test", "dash", dash = TRUE),
    N("Symptoms: early", "SEER stage mix", "out", kids = list(N("Cancer death", paste(pct(P$m_early), "/ y"), "harm"))),
    N("Symptoms: late", "", "out", kids = list(N("Cancer death", paste(pct(P$m_late), "/ y"), "harm"))),
    N("Other-cause death first", "overdiagnosed if found", "harm"))),
  N("Other-cause death", "NCHS 2024, cause-deleted", "out")))

# EAC natural history in the MISCAN-ESO scenario (rates from outputs/miscan_eac/rates.csv)
ER <- rd("outputs/miscan_eac/rates.csv"); erm <- ER[ER$sex == "male", ]
rt <- function(x) sprintf("%.3f / y", x)
eac <- N("Preclinical EAC, early", sprintf("mean %.1f y (men)", erm$mean_years_in_PE), "dec", kids = list(
  N("Symptoms, early stage", sprintf("cancer death %s / y", pct(P$m_early_eac)), "out", el = rt(erm$rate_PE_to_symptoms_early)),
  N("Preclinical EAC, late", sprintf("mean %.1f y", erm$mean_years_in_PL), "dec", el = rt(erm$rate_PE_to_PL), kids = list(
    N("Symptoms, late stage", sprintf("cancer death %s / y", pct(P$m_late_eac)), "harm", el = rt(erm$rate_PL_to_symptoms_late)))),
  N("Screen-detected early", "lead time = remaining sojourn", "found", el = "Galleri / EGD", dash = TRUE)))

figs <- list(arms = tree_svg(arms, "Three-arm decision tree"), eac = tree_svg(eac, "Esophageal adenocarcinoma natural history, MISCAN-ESO scenario"), egd = tree_svg(egd, "Galleri positive: the EGD step"),
  rest = tree_svg(rest, "After a negative or skipped EGD: ENT or PET CT"),
  pb1 = tree_svg(pb1, "First pancreatobiliary finding"), pb2 = tree_svg(pb2, "Repeat pancreatobiliary finding"),
  nh = tree_svg(nh, "Natural history states"))
fig <- function(k) sprintf('<figure class="canvas" style="--minw:%dpx">%s</figure>', figs[[k]]$W, figs[[k]]$svg)

# ---- numbers -----------------------------------------------------------------
icer <- g("icer"); ui <- gl("icer")
GV <- "Galleri vs usual care"
ms_icer <- ms(GV, "icer")
kp <- list(
  c(money(icer), sprintf("Galleri cost per QALY gained, cohort model (95%% UI %s–%s)", money(ui[1]), money(ui[2]))),
  c(money(ms_icer), "Galleri cost per QALY gained, microsimulation of 1M people"),
  c(n(g("ca_deaths_averted")), "cancer deaths averted per 100,000"),
  c(n(g("overdiagnosed")), sprintf("cancers overdiagnosed (%s of those found)", pct(g("overdiagnosis_share")))),
  c(n(g("free_investigated")), "cancer-free people worked up after a Galleri +ve"),
  c("Dominated", "one-time Kunzmann → EGD (fewer QALYs, more cost)"))
kpis <- paste(vapply(kp, function(z) sprintf('<div class="kpi"><span class="v">%s</span><span class="l">%s</span></div>', z[1], esc(z[2])), ""), collapse = "\n")

row3 <- function(lbl, a, b) sprintf("<tr><td>%s</td><td class='num'>%s</td><td class='num'>%s</td></tr>", lbl, a, b)
kz_inc_cost <- if (!is.na(kz("inc_cost"))) kz("inc_cost") else kz("cost_total") - getv(CL, "Usual care", "cost_total")
res_rows <- paste0(
  row3("Cost per QALY gained", sprintf("%s (95%% UI %s–%s)", money(icer), money(ui[1]), money(ui[2])), "dominated"),
  row3("Incremental cost per person", money(g("inc_cost")), money(kz_inc_cost)),
  row3("Net QALYs", paste0("+", n(g("qaly_net"))), n(kz("qaly_net"))),
  row3("Cancers found by screening", sprintf("%s (%s upper GI)", n(g("screen_detected")), n(g("found_upper_gi"))), n(kz("screen_detected"))),
  row3("Cancer deaths averted", n(g("ca_deaths_averted")), n(kz("ca_deaths_averted"), 1)),
  row3("Late-stage diagnoses averted", n(g("late_averted")), n(kz("late_averted"), 1)),
  row3("Overdiagnosed cancers", sprintf("%s (%s)", n(g("overdiagnosed")), pct(g("overdiagnosis_share"))), n(kz("overdiagnosed"), 1)),
  row3("EGDs", n(g("egd")), n(kz("egd"))),
  row3("EGD perforations", n(g("perforations"), 2), n(kz("perforations"), 1)),
  row3("Positives inside the EGD hold", n(g("egd_skipped_hold")), "—"),
  row3("Cancers missed after +ve: visible / not yet visible", sprintf("%s / %s", n(g("missed_visible_after_egd")), n(g("hidden_after_egd"))), "—"),
  row3("PET CT · EUS · laparoscopy · watchful waiting",
       sprintf("%s · %s · %s · %s", n(g("pet_ct")), n(g("eus")), n(g("laparoscopy"), 1), n(g("surveillance_imaging"))), "—"),
  row3("False-positive burden", sprintf("%s people; %s days each; %s QALYs, %s DALYs", n(g("free_investigated")), n(g("mean_days_fp")),
       n(g("qaly_fp"), 1), n(g("daly_fp"), 1)), "—"))

# microsimulation vs cohort (Galleri arm)
vrows <- list(c("Upper GI cancers found by screening", "screen_detected", 0), c("Symptomatic diagnoses", "clin_dx", 0),
  c("Overdiagnosed (upper GI)", "overdiagnosed", 0), c("Cancer deaths", "ca_death", 0),
  c("Cancer-free people worked up", "false_pos", 0), c("EGDs", "egd", 0), c("PET CT", "pet", 0), c("EUS", "eus", 0),
  c("Laparoscopies", "laparoscopy", 1), c("QALYs per person", "qaly", 4), c("Lifetime cost per person", "cost", -1))
vt <- paste(vapply(vrows, function(z) {
  r <- MV[MV$arm == "Galleri" & MV$measure == z[2], ]; a <- as.numeric(r$microsim); b <- as.numeric(r$cohort); d <- as.integer(z[3])
  f <- if (d == -1) money else function(x) n(x, d)
  diff <- if (b != 0) sprintf("%+.1f%%", 100 * (a - b) / b) else "—"
  sprintf("<tr><td>%s</td><td class='num'>%s</td><td class='num'>%s</td><td class='num'>%s</td></tr>", z[1], f(b), f(a), diff)
}, ""), collapse = "")

nf <- function(m) as.numeric(mf(m))
facts <- list(
  c("Better off / worse off / unaffected", sprintf("%s / %s / %s", pct(ms(GV, "share_gain"), 2), pct(ms(GV, "share_loss"), 2), pct(ms(GV, "share_same"), 2))),
  c("Upper GI cancers found by Galleri", sprintf("%.0f%% of lifetime cancers", nf("Galleri finds it, % of those cancers"))),
  c("Median lead time", sprintf("%.1f years (IQR %s)", nf("Lead time, median years (screen-detected, presented in usual care)"), mf("Lead time, IQR"))),
  c("Shifted late → early", sprintf("%.0f%% of screen-detected", nf("Late stage in usual care -> early when found by Galleri, % of screen-detected"))),
  c("Ever positive", sprintf("%.2f%% of people", nf("People with any positive test, %"))),
  c("Largest single-person change", sprintf("+%.1f / %.1f QALYs", nf("Galleri: largest single-person QALY gain"), nf("Galleri: largest single-person QALY loss"))),
  c("Lifetime cost, 99th percentile", money(nf("Cost per person, 99th percentile (Galleri)"))))
facts <- paste(vapply(facts, function(z) sprintf("<tr><td>%s</td><td class='num'>%s</td></tr>", z[1], z[2]), ""), collapse = "")
pos <- MD[MD$measure == "Positive Galleri tests per person" & MD$value != "0", ]
posbar <- paste(sprintf("<tr><td>%s positive%s</td><td class='num'>%s</td></tr>", pos$value, ifelse(pos$value != "1", "s", ""),
                        vapply(1e5 * as.numeric(pos$share), n, "", d = 1)), collapse = "")

# scenarios, sorted by ICER
sc <- SC[!is.na(SC$icer), c("scenario", "icer")]; sc <- sc[order(sc$icer), ]
scrows <- paste0(paste(sprintf("<tr%s><td>%s</td><td class='num'>%s</td></tr>", ifelse(sc$scenario == "Base case", ' class="base"', ""),
                               esc(sc$scenario), vapply(sc$icer, money, "")), collapse = ""),
                 sprintf("<tr><td>Microsimulation (continuous time, 1M people)</td><td class='num'>%s</td></tr>", money(ms_icer)))

# inputs that need review
rv_keys <- c("p_pet_ugi", "p_lap_cx", "t_lap", "c_lap", "c_surv", "p_eus_samp", "p_tb_lap", "r_lap_after_eus", "r_lap_after_lap",
             "p_reeus_after_eus", "p_reeus_after_lap", "p_pb_malig", "t_visible", "sens_g_ugi", "fpr_g_es", "fp_kappa", "eac_dwell_ratio")
rv <- paste(vapply(rv_keys, function(k) {
  r <- PST[[k]]; st <- r$status
  sprintf("<tr class='%s'><td>%s</td><td class='num'>%s</td><td><span class='pill %s'>%s</span></td><td><code>%s</code></td></tr>",
          if (st == "DUMMY") "flag-review" else "", esc(r$label), esc(r$value),
          if (st == "DUMMY") "gap" else if (st == "ASSUMPTION") "data" else "ok", st, k)
}, ""), collapse = "")

# MISCAN-ESO targets and scenario results
ET <- rd("outputs/miscan_eac/targets.csv"); ET <- ET[ET$kind != "out of scope" & !grepl("by age", ET$target), ]
fmt_t <- function(x, tgt) if (is.na(x)) "\u2014" else if (grepl("share|per EAC case", tgt)) pct(x) else if (grepl("sojourn", tgt)) sprintf("%.1f y", x) else n(x)
eac_targets <- paste(vapply(seq_len(nrow(ET)), function(i) sprintf("<tr><td>%s</td><td>%s</td><td class='num'>%s</td><td class='num'>%s</td></tr>",
  esc(sub(", NH White", "", ET$target[i])), ET$kind[i], fmt_t(ET$target_value[i], ET$target[i]), fmt_t(ET$model_value[i], ET$target[i])), ""), collapse = "")
es <- SC[grepl("^MISCAN-ESO|^Base case$", SC$scenario), c("scenario", "icer")]
eac_scen <- paste(sprintf("<tr%s><td>%s</td><td class='num'>%s</td></tr>", ifelse(es$scenario == "Base case", ' class="base"', ""),
                          esc(es$scenario), vapply(es$icer, money, "")), collapse = "")

today <- sub("^0", "", format(Sys.Date(), "%d %b %Y"))
CSS <- paste0(paste(readLines("viz/model_page.css", warn = FALSE, encoding = "UTF-8"), collapse = "\n"), "\n")
kappa <- formatC(P$fp_kappa, format = "d", big.mark = ",")

tpl <- paste(readLines("viz/model_page_template.html", warn = FALSE, encoding = "UTF-8"), collapse = "\n")
fill <- list(CSS = CSS, TODAY = today, STOP = S$screen_stop_age, KPIS = kpis, RES_ROWS = res_rows,
  FIG_ARMS = fig("arms"), FIG_EGD = fig("egd"), FIG_REST = fig("rest"), FIG_PB1 = fig("pb1"), FIG_PB2 = fig("pb2"), FIG_NH = fig("nh"),
  KSHARE = pct(kshare), KCAP = pct(kcap), SENS_G = pct(P$sens_g_ugi), REP = num(P$sens_g_rep_ratio),
  HOLD = hold, SE = pct(se), FI = pct(fi, 0), TVIS = num(P$t_visible), PETUGI = pr(pet_ugi),
  TL = pct(tl, 0), PL2 = pr(pl2), PE2 = pct(pe2, 0), PL3 = pr(pl3), PE3 = pct(pe3, 0),
  FPR = pct(P$fpr_g_es, 3), KAPPA = kappa,
  MS_ICER = money(ms_icer), MS_COHORT = money(ms(GV, "cohort_icer")), MS_QALY = sprintf("%.4f", ms(GV, "qaly")),
  MS_QALY95 = sprintf("%.4f", ms(GV, "qaly_mc95")), MS_COST = money(ms(GV, "cost")),
  MS_LOSS = pct(ms(GV, "share_loss"), 2), MS_GAIN = pct(ms(GV, "share_gain"), 2),
  VT = vt, FACTS = facts, POSBAR = posbar, SCROWS = scrows, RV = rv,
  FIG_EAC = fig("eac"), EAC_S_M = num(P$sojourn_eac_male), EAC_S_F = num(P$sojourn_eac_female),
  EAC_TARGETS = eac_targets, EAC_SCEN = eac_scen)
for (k in names(fill)) tpl <- gsub(paste0("{{", k, "}}"), as.character(fill[[k]]), tpl, fixed = TRUE)
stopifnot(!grepl("{{", tpl, fixed = TRUE))
writeLines(tpl, "viz/galleri_model.html", useBytes = TRUE)
cat("wrote viz/galleri_model.html", nchar(tpl), "\n")
