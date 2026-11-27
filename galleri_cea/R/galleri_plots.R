# =============================================================
# Galleri project figures. Signature: function(res, p, data) returning a
# ggplot or drawing with base graphics (used by the app and write_cea()).
# =============================================================

STATE_COLS <- c("Well" = "#C9DCEB", "Undiagnosed cancer" = "#F2B701", "Early-stage dx" = "#4C78A8",
                "Late-stage dx" = "#E45756", "Cancer death" = "#3B3B3B", "Other-cause death" = "#9E9E9E")

dummy_note <- function(res) if (any(res$params$status == "DUMMY", na.rm = TRUE)) " (DUMMY inputs)" else ""

plot_workup_tree <- function(res, p, data, cohort = 1000) {
  g <- res$base[, "Galleri"]
  egd_pos <- g[["egd_positive"]] / max(1e-15, g[["egd_positive"]] + g[["egd_negative_workup"]])
  d <- rest_summary(p, 1)$paths
  E <- unique(do.call(rbind, lapply(strsplit(d$path, " > ", fixed = TRUE), function(z)
    data.frame(from = head(z, -1), to = z[-1]))))
  rt <- rest_tree_edges(p)
  E$prob <- rt$prob[match(paste(E$from, E$to), paste(rt$from, rt$to))]
  E <- rbind(data.frame(from = c("galleri", "egd", "egd"), to = c("egd", "egd_pos", "start"), prob = c(NA, egd_pos, 1 - egd_pos)), E)
  lab <- function(id) {
    t <- sub("\\..*$", "", id)
    unname(c(galleri = sprintf("Galleri +ve, esophagus/stomach (n = %s), outside the EGD hold", ceakit::fmt(cohort)),
             egd = "Upper GI endoscopy (EGD)", egd_pos = "EGD +ve: cancer found", start = "EGD -ve: rest of workup",
             ent = "ENT exam", pet = "PET CT", ent_pos = "ENT +ve", ent_neg = "ENT -ve", pb = "Pancreato-biliary +ve",
             incid = "Incidental finding", pet_neg = "PET -ve", eus = "EUS", samp = "FNA/FNB", nosamp = "No sample",
             euspos = "Sample +ve", eusneg = "Sample -ve", tb = "Tumor board", lap = "Laparo-scopy", lappos = "Lap +ve",
             lapneg = "Lap -ve", ww = "Watchful waiting")[t])
  }
  ids <- unique(c(E$from, E$to)); kids <- split(E$to, E$from)
  depth <- setNames(numeric(length(ids)), ids); x <- setNames(numeric(length(ids)), ids); nx <- 0
  place <- function(id, dpt) {
    depth[id] <<- dpt; ch <- kids[[id]]
    if (is.null(ch)) { nx <<- nx + 1; x[id] <<- nx } else { for (c2 in ch) place(c2, dpt + 1); x[id] <<- mean(x[ch]) }
  }
  place("galleri", 0)
  nodes <- data.frame(id = ids, label = vapply(ids, lab, ""), x = x[ids] * 1.55, y = -depth[ids] * 1.15, w = 1.45,
                      stringsAsFactors = FALSE)
  nodes$w[nodes$id == "galleri"] <- 5
  found <- c("egd_pos", "ent_pos", "euspos.a", "lappos.n", "lappos.s", "incid")
  nodes$fill <- ifelse(nodes$id %in% found, "#9CC3E6", ifelse(startsWith(nodes$id, "ww") | startsWith(nodes$id, "tb"), "#D9D2E9", "#FCCD2A"))
  ceakit::draw_tree(nodes, E,
    sprintf("Galleri +ve workup (first pancreatobiliary finding)%s", dummy_note(res)),
    sprintf(paste("Blue = cancer can be found; purple = tumor board and watchful waiting. Inside the %s-year hold after a negative EGD,\n",
                  "positives skip the EGD and start at 'rest of workup'. After a prior negative EUS or laparoscopy,\n",
                  "a repeat pancreatobiliary finding leads to laparoscopy %s%% / %s%%, repeat EUS %s%% / %s%%, otherwise watchful waiting."),
            p$egd_hold_years, round(100 * p$p_tb_lap * p$r_lap_after_eus, 1), round(100 * p$p_tb_lap * p$r_lap_after_lap, 1),
            round(100 * p$p_reeus_after_eus), round(100 * p$p_reeus_after_lap)),
    wrap = 7, cex = 0.62, pcex = 0.55)
}

plot_strategy_tree <- function(res, p, data) {
  b <- res$base; a <- res$arms
  cb <- calibration_table(p, data, res$settings, a[1])
  pl_clin <- sum(cb$seer_late_share * cb$seer * cb$share) / sum(cb$seer * cb$share)
  per <- function(x, y) if (y > 0) x / y else 0
  co <- cohort_rows(data)
  pk <- sum(co$share[co$kunzmann_scoped == 1]) / sum(co$share)               # share in the Kunzmann high-risk group
  k <- b[, a[2]]; g <- b[, a[3]]
  pk_pos <- per(k[["egd_positive"]], k[["investigated"]])                     # per screening EGD
  g_pos <- per(g[["screen_tp"]] + g[["screen_fp"]], g[["galleri_tests"]])    # per Galleri test
  g_tn <- per(g[["screen_tn"]], g[["screen_tn"]] + g[["screen_fn"]])
  g_egd_pos <- per(g[["egd_positive"]], g[["egd_positive"]] + g[["egd_negative_workup"]])
  cols <- ceakit::arm_colours(a)
  nodes <- data.frame(
    id = c("root", "s1", "s2", "s3", "uc", "uc_late", "uc_early", "k_hi", "k_lo", "k_egd", "k_pos", "k_neg",
           "g_pos", "g_neg", "g_egd", "g_tp", "g_fpw", "g_perf", "g_tn", "g_fn", "mk"),
    label = c("Cohort eligible for upper GI cancer screening", a,
              "Cancers present with symptoms (SEER-calibrated)", "Late stage", "Early stage",
              "Kunzmann high-risk group(s)", "Other groups: usual care", "Upper GI endoscopy",
              "EGD +ve: screen-detected", "EGD -ve",
              "Galleri +ve (esophagus/ stomach signal)", "Galleri -ve",
              "Upper GI endoscopy (EGD)", "EGD +ve: true positive", "EGD -ve: false positive -> further workup",
              "Perforation", "True negative: stays well, rescreened next round",
              "False negative: symptomatic cancer later (or found at a later round)",
              "Every end state -> lifetime Markov model: other-cause death, new cancers, symptomatic presentation, cancer death"),
    x = c(12.5, 3, 12, 22.5, 3, 1.5, 4.5, 10.4, 13.6, 10.4, 8.9, 11.9, 19.6, 25.4, 19.6, 16.6, 19.6, 22.4, 24.0, 26.8, 13.6),
    y = c(9, 7.6, 7.6, 7.6, 6.2, 4.8, 4.8, 6.2, 6.2, 4.8, 3.4, 3.4, 6.2, 6.2, 4.8, 3.4, 3.4, 3.4, 4.8, 4.8, 1.9),
    w = c(5.2, 3.4, 3.4, 3.4, 3.4, 2.0, 2.0, 2.9, 2.9, 2.9, 2.8, 2.8, 2.9, 2.9, 2.9, 2.6, 2.8, 2.2, 2.6, 2.6, 10.0),
    stringsAsFactors = FALSE)
  arm_ids <- c("s1", "s2", "s3")
  nodes$fill <- ifelse(nodes$id == "root", "#3B3B3B", ifelse(nodes$id == "mk", "#E8E8E8",
                ifelse(nodes$id %in% c("g_tp", "k_pos"), "#9CC3E6", ifelse(nodes$id == "g_perf", "#F4B183", "#FCCD2A"))))
  nodes$fill[match(arm_ids, nodes$id)] <- unname(cols[a])
  nodes$tcol <- ifelse(nodes$id %in% c("root", arm_ids), "white", "grey10")
  edges <- data.frame(
    from = c("root", "root", "root", "s1", "uc", "uc", "s2", "s2", "k_hi", "k_egd", "k_egd",
             "s3", "s3", "g_pos", "g_egd", "g_egd", "g_egd", "g_neg", "g_neg"),
    to = c("s1", "s2", "s3", "uc", "uc_late", "uc_early", "k_hi", "k_lo", "k_egd", "k_pos", "k_neg",
           "g_pos", "g_neg", "g_egd", "g_tp", "g_fpw", "g_perf", "g_tn", "g_fn"),
    prob = c(NA, NA, NA, NA, pl_clin, 1 - pl_clin, pk, 1 - pk, NA, pk_pos, 1 - pk_pos,
             g_pos, 1 - g_pos, NA, g_egd_pos, 1 - g_egd_pos, p$p_perf_egd, g_tn, 1 - g_tn))
  ceakit::draw_tree(nodes, edges, paste0("Screening strategies compared (per screening round, averaged over rounds)", dummy_note(res)),
                    "Galleri +ve/-ve per test; TP/FP/TN/FN against preclinical upper GI cancer. Perforation can occur with either EGD result (see the Galleri workup tree).",
                    xlim = c(0, 28.4), ylim = c(1.2, 9.5), wrap = 7.2, cex = 0.72, pcex = 0.62)
}

plot_balance <- function(res, p, data, arm = utils::tail(res$arms, 1)) {
  keys <- c(screen_detected = "Benefit", late_averted = "Benefit", ca_deaths_averted = "Benefit",
            free_investigated = "Harm", invasive_free = "Harm", overdiagnosed = "Harm", missed_investigated = "Harm",
            sae_total = "Harm", proc_deaths = "Harm", radiation_cancers = "Harm")
  d <- res$long[res$long$arm == arm & res$long$key %in% names(keys), ]
  d$type <- keys[d$key]
  d$metric <- sub("^- ", "", d$metric); d$metric <- paste0(toupper(substr(d$metric, 1, 1)), substring(d$metric, 2))
  d$metric <- factor(d$metric, rev(d$metric))
  q <- res$long[res$long$arm == arm & res$long$key %in% c("qaly_net", "daly_averted"), ]
  ggplot2::ggplot(d, ggplot2::aes(value, metric, fill = type)) + ggplot2::geom_col(width = 0.7) +
    { if (!all(is.na(d$lower_95))) ggplot2::geom_errorbar(ggplot2::aes(xmin = lower_95, xmax = upper_95),
                                                         width = 0.2, colour = "grey30", orientation = "y") } +
    ggplot2::geom_text(ggplot2::aes(x = pmax(value, dplyr::coalesce(upper_95, value)),
                                    label = scales::comma(value, accuracy = ifelse(value < 1, 0.01, 1))),
                       hjust = -0.2, size = 3.2) +
    ggplot2::scale_fill_manual(values = c(Benefit = "#4C78A8", Harm = "#F58518"), name = NULL) +
    ggplot2::scale_x_continuous(labels = scales::comma, expand = ggplot2::expansion(c(0, 0.15))) +
    ggplot2::labs(y = NULL, x = res$analysis$unit, title = paste0("Benefits and harms: ", arm, dummy_note(res)),
                  subtitle = sprintf("Net QALYs %s | net DALYs averted %s | vs %s",
                                     scales::comma(q$value[q$key == "qaly_net"], 0.1),
                                     scales::comma(q$value[q$key == "daly_averted"], 0.1), res$arms[1])) +
    ceakit::theme_cea()
}

plot_icons <- function(res, p, data, arm = "Galleri") {
  m <- res$base[, arm]; n <- m[["investigated"]]
  if (n <= 0) return(NULL)
  g <- c("Cancer found: upper GI" = m[["found_upper_gi"]], "Cancer found: other site" = m[["found_other"]],
         "Cancer missed despite workup" = m[["missed_investigated"]],
         "No cancer: false-positive Galleri" = m[["free_investigated"]]) / n * 1000
  f <- floor(g); r <- 1000 - sum(f)
  if (r > 0) { o <- order(g - f, decreasing = TRUE)[seq_len(r)]; f[o] <- f[o] + 1 }
  d <- tibble::tibble(i = 0:999, x = i %% 50, y = -(i %/% 50), group = factor(rep(names(g), f), names(g)))
  ggplot2::ggplot(d, ggplot2::aes(x, y, fill = group)) + ggplot2::geom_tile(width = 0.85, height = 0.85) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_manual(values = c("#1F4E79", "#6FA8DC", "#E8590C", "#C9C9C9"),
                               labels = sprintf("%s (%d)", names(g), f), name = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::labs(title = paste0("What happens to every 1,000 people with a Galleri +ve esophagus/stomach signal", dummy_note(res))) +
    ggplot2::theme_void(base_size = 12) + ggplot2::theme(legend.position = "bottom")
}

plot_trace <- function(res, p, data) {
  d <- trace_table(p, data, res$settings, res$arms)
  d$arm <- factor(d$arm, res$arms)
  ggplot2::ggplot(d, ggplot2::aes(cycle, prop, fill = state)) + ggplot2::geom_area(colour = NA) +
    ggplot2::facet_wrap(~arm) + ggplot2::scale_fill_manual(values = STATE_COLS, name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Years after first screen", y = "Share of cohort",
                  title = "Markov state occupancy over the lifetime (competing risks)") + ceakit::theme_cea()
}

plot_cancer_states <- function(res, p, data) {
  d <- trace_table(p, data, res$settings, res$arms)
  d <- d[d$state %in% c("Undiagnosed cancer", "Early-stage dx", "Late-stage dx", "Cancer death"), ]
  d$arm <- factor(d$arm, res$arms)
  ggplot2::ggplot(d, ggplot2::aes(cycle, prop * res$analysis$cohort_n, colour = arm)) +
    ggplot2::geom_line(linewidth = 1) + ggplot2::facet_wrap(~state, scales = "free_y") +
    ggplot2::scale_colour_manual(values = ceakit::arm_colours(res$arms), name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(x = "Years after first screen", y = res$analysis$unit, title = "Cancer states over time by arm") +
    ceakit::theme_cea()
}

# Subgroup sensitivity grid: start age x sex x race/ethnicity (deterministic)
plot_subgroups <- function(res, p, data) {
  d <- subgroup_grid(p, data, res$settings, res$arms, res$analysis$cohort_n)
  d <- d[d$arm != res$arms[1], ]
  icer <- d$icer_vs_usual_care
  d$band <- dplyr::case_when(
    d$start_age > res$settings$screen_stop_age ~ "Not screened (past stop age)",
    d$inc_qaly <= 0 ~ "Dominated (fewer QALYs)",
    icer < 1e5 ~ "< $100k", icer < 5e5 ~ "$100k-$500k", icer < 1e6 ~ "$500k-$1M",
    icer < 5e6 ~ "$1M-$5M", TRUE ~ "> $5M")
  d$band <- factor(d$band, c("< $100k", "$100k-$500k", "$500k-$1M", "$1M-$5M", "> $5M",
                             "Dominated (fewer QALYs)", "Not screened (past stop age)"))
  d$label <- ifelse(d$band %in% c("Dominated (fewer QALYs)", "Not screened (past stop age)"), "",
                    ifelse(icer < 1e6, scales::dollar(icer, scale = 1e-3, accuracy = 1, suffix = "k"),
                           scales::dollar(icer, scale = 1e-6, accuracy = 0.1, suffix = "M")))
  d$race_label <- factor(d$race_label, rev(unname(RACE_LABELS)))
  d$sex <- factor(d$sex, c("male", "female"), c("Men", "Women"))
  ggplot2::ggplot(d, ggplot2::aes(factor(start_age), race_label, fill = band)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5) +
    ggplot2::facet_grid(arm ~ sex) +
    ggplot2::scale_fill_manual(values = c("< $100k" = "#1A9850", "$100k-$500k" = "#91CF60",
                                          "$500k-$1M" = "#D9EF8B", "$1M-$5M" = "#FEE08B", "> $5M" = "#FC8D59",
                                          "Dominated (fewer QALYs)" = "#D73027",
                                          "Not screened (past stop age)" = "grey85"),
                               drop = FALSE, name = "Cost per QALY vs usual care") +
    ggplot2::labs(x = "Age at first screen", y = NULL,
                  title = "Subgroup sensitivity analysis: cost per QALY gained vs usual care",
                  subtitle = sprintf("Annual screening to age %d; natural history calibrated to SEER for each subgroup",
                                     res$settings$screen_stop_age)) +
    ceakit::theme_cea() + ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))
}

# Calibration: modelled symptomatic incidence vs SEER target, by age
plot_calibration <- function(res, p, data) {
  d <- calibration_table(p, data, res$settings, res$arms)
  d <- d |> dplyr::group_by(arm, age) |>
    dplyr::summarise(seer = sum(seer * share), model = sum(model_clinical * share), .groups = "drop")
  d$arm <- factor(d$arm, res$arms)
  ggplot2::ggplot(d, ggplot2::aes(age)) +
    ggplot2::geom_point(ggplot2::aes(y = seer * 1e5, shape = "SEER 2013-2017 (target)"), colour = "grey30", size = 1.6) +
    ggplot2::geom_line(ggplot2::aes(y = model * 1e5, colour = arm), linewidth = 1) +
    ggplot2::scale_colour_manual(values = ceakit::arm_colours(res$arms), name = "Model: symptomatic diagnoses") +
    ggplot2::scale_shape_manual(values = 16, name = NULL) +
    ggplot2::labs(x = "Age", y = "Upper GI cancers per 100,000 per year",
                  title = "Calibration to SEER: esophagus + stomach incidence",
                  subtitle = "Usual care reproduces SEER; screening arms divert cancers to screen detection") +
    ceakit::theme_cea()
}
