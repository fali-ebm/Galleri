# =============================================================
# Galleri project figures. Signature: function(res, p, data) returning a
# ggplot or drawing with base graphics (used by the app and write_cea()).
# =============================================================

STATE_COLS <- c("Well" = "#C9DCEB", "Undiagnosed cancer" = "#F2B701", "Early-stage dx" = "#4C78A8",
                "Late-stage dx" = "#E45756", "Cancer death" = "#3B3B3B", "Other-cause death" = "#9E9E9E")

dummy_note <- function(res) if (any(res$params$status == "DUMMY", na.rm = TRUE)) " (DUMMY inputs)" else ""

plot_workup_tree <- function(res, p, data, cohort = 1000) {
  paths <- workup_paths(p)
  nodes <- data.frame(
    id = c("galleri", "egd", "egd_neg", "egd_pos", "ent", "pet", "ent_neg", "ent_pos", "pb_pos", "incid",
           "pet_neg", "tb1", "eus", "tb2", "eus_pos", "eus_neg", "no_eus", "tb3"),
    label = c(sprintf("Galleri +ve for esophagus/stomach (n = %s)", ceakit::fmt(cohort)), "Upper GI endoscopy",
              "EGD -ve", "EGD +ve", "ENT evaluation", "PET CT", "ENT -ve", "ENT +ve", "Pancreatobiliary +ve",
              "Incidental non-GI/non-PB +ve", "PET -ve", "Tumor board discussion/workup complete", "EUS",
              "Tumor board discussion/workup complete", "EUS sampling +ve", "EUS sampling -ve",
              "No EUS sampling", "Tumor board discussion/workup complete"),
    x = c(10.1, 10.1, 7.75, 12.75, 3.45, 12.0, 1.6, 5.4, 10.9, 14.75, 18.4, 3.65, 10.9, 15.9, 8.55, 12.0, 15.5, 13.5),
    y = c(9, 8, 7, 7, 6, 6, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 3, 2),
    w = c(4.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.8, 3.2, 5.8, 3.2, 5.8, 3.2, 3.2, 3.2, 5.8),
    stringsAsFactors = FALSE)
  n_end <- stats::setNames(paths$prob * cohort, paths$leaf)
  term <- nodes$id %in% names(n_end)
  nodes$label[term] <- sprintf("%s [%s]", nodes$label[term], ceakit::fmt(n_end[nodes$id[term]]))
  nodes$fill <- ifelse(nodes$id %in% c("egd_pos", "ent_pos", "eus_pos", "incid"), "#9CC3E6", "#FCCD2A")
  edges <- data.frame(
    from = c("galleri", "egd", "egd", "egd_neg", "egd_neg", "ent", "ent", "ent_neg", "pet", "pet", "pet",
             "pb_pos", "incid", "pet_neg", "eus", "eus", "eus", "eus_neg", "no_eus"),
    to = c("egd", "egd_neg", "egd_pos", "ent", "pet", "ent_neg", "ent_pos", "tb1", "pb_pos", "incid", "pet_neg",
           "eus", "tb2", "tb2", "eus_pos", "eus_neg", "no_eus", "tb3", "tb3"),
    prob = c(NA, 1 - p$p_egd_pos, p$p_egd_pos, p$p_ent_route, 1 - p$p_ent_route, 1 - p$p_ent_pos, p$p_ent_pos, NA,
             p$p_pet_pb, p$p_pet_incid, 1 - p$p_pet_pb - p$p_pet_incid, NA, NA, NA,
             p$p_eus_pos, p$p_eus_neg, 1 - p$p_eus_pos - p$p_eus_neg, NA, NA))
  f <- sum(paths$prob * paths$found); m <- sum(paths$prob * paths$missed)
  ceakit::draw_tree(nodes, edges,
    sprintf("Galleri +ve workup: branch probabilities and [people per %s]%s", ceakit::fmt(cohort), dummy_note(res)),
    sprintf("Blue = cancer found. Found %s | missed despite workup %s | cancer-free worked up %s. Each end state enters the lifetime Markov model.",
            ceakit::fmt(f * cohort), ceakit::fmt(m * cohort), ceakit::fmt((1 - f - m) * cohort)),
    xlim = c(-0.2, 20.2), ylim = c(1.5, 9.5))
}

plot_strategy_tree <- function(res, p, data) {
  b <- res$base; a <- res$arms
  pk <- b["investigated", a[2]]; pg <- b["investigated", a[4]]
  pk_pos <- if (pk > 0) b["screen_detected", a[2]] / pk else 0
  pg_pos <- if (pg > 0) b["screen_detected", a[3]] / pg else 0
  cols <- ceakit::arm_colours(a)
  nodes <- data.frame(
    id = c("root", "s1", "s2", "s3", "s4", "uc", "uc_late", "uc_early", "k_hi", "k_lo", "k_egd", "k_pos", "k_neg",
           "g3_pos", "g3_neg", "g3_egd", "g3_ep", "g3_en", "g4_pos", "g4_neg", "g4_wk", "mk"),
    label = c("Cohort eligible for upper GI cancer screening", a,
              "Prevalent cancers stay undiagnosed until symptoms", "Late stage", "Early stage",
              "Kunzmann high-risk group(s)", "Other groups: usual care", "Upper GI endoscopy",
              "EGD +ve: screen-detected", "EGD -ve", "Galleri +ve (esophagus/ stomach signal)",
              "Galleri -ve: usual care", "Upper GI endoscopy", "EGD +ve: screen-detected", "EGD -ve: stop",
              "Galleri +ve (esophagus/ stomach signal)", "Galleri -ve: usual care", "Full workup (see workup tree)",
              "Every end state -> lifetime Markov model: other-cause death, new cancers, symptomatic presentation, cancer death"),
    x = c(12, 2, 8, 14.4, 20.6, 2, 1.0, 3.2, 6.4, 9.6, 6.4, 4.9, 7.9, 12.8, 16.0, 12.8, 11.3, 14.3, 19.0, 22.2, 19.0, 19.5),
    y = c(9, 7.6, 7.6, 7.6, 7.6, 6.2, 4.8, 4.8, 6.2, 6.2, 4.8, 3.4, 3.4, 6.2, 6.2, 4.8, 3.4, 3.4, 6.2, 6.2, 4.8, 3.3),
    w = c(5.2, 3.4, 3.4, 3.4, 3.4, 3.4, 2.0, 2.0, 2.9, 2.9, 2.9, 2.8, 2.8, 2.9, 2.9, 2.9, 2.8, 2.8, 2.9, 2.9, 2.9, 7.6),
    stringsAsFactors = FALSE)
  arm_ids <- c("s1", "s2", "s3", "s4")
  nodes$fill <- ifelse(nodes$id == "root", "#3B3B3B", ifelse(nodes$id == "mk", "#E8E8E8", "#FCCD2A"))
  nodes$fill[match(arm_ids, nodes$id)] <- unname(cols[a])
  nodes$tcol <- ifelse(nodes$id %in% c("root", arm_ids), "white", "grey10")
  edges <- data.frame(
    from = c("root", "root", "root", "root", "s1", "uc", "uc", "s2", "s2", "k_hi", "k_egd", "k_egd",
             "s3", "s3", "g3_pos", "g3_egd", "g3_egd", "s4", "s4", "g4_pos"),
    to = c("s1", "s2", "s3", "s4", "uc", "uc_late", "uc_early", "k_hi", "k_lo", "k_egd", "k_pos", "k_neg",
           "g3_pos", "g3_neg", "g3_egd", "g3_ep", "g3_en", "g4_pos", "g4_neg", "g4_wk"),
    prob = c(NA, NA, NA, NA, NA, p$p_late_clin, 1 - p$p_late_clin, pk, 1 - pk, NA, pk_pos, 1 - pk_pos,
             pg, 1 - pg, NA, pg_pos, 1 - pg_pos, pg, 1 - pg, NA))
  ceakit::draw_tree(nodes, edges, paste0("Screening strategies compared (cohort-weighted)", dummy_note(res)),
                    xlim = c(0, 24), ylim = c(2.6, 9.5), wrap = 7.2, cex = 0.76, pcex = 0.66)
}

plot_balance <- function(res, p, data, arm = utils::tail(res$arms, 1)) {
  keys <- c(screen_detected = "Benefit", late_averted = "Benefit", ca_deaths_averted = "Benefit",
            free_investigated = "Harm", invasive_free = "Harm", missed_investigated = "Harm",
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

plot_icons <- function(res, p, data, arm = "Full workup") {
  m <- res$base[, arm]
  g <- c("Cancer found: upper GI" = m[["found_upper_gi"]], "Cancer found: other site" = m[["found_other"]],
         "Cancer missed despite workup" = m[["missed_investigated"]],
         "No cancer: false-positive Galleri" = m[["cancer_free"]]) * 1000
  f <- floor(g); r <- 1000 - sum(f)
  if (r > 0) { o <- order(g - f, decreasing = TRUE)[seq_len(r)]; f[o] <- f[o] + 1 }
  d <- tibble::tibble(i = 0:999, x = i %% 50, y = -(i %/% 50), group = factor(rep(names(g), f), names(g)))
  ggplot2::ggplot(d, ggplot2::aes(x, y, fill = group)) + ggplot2::geom_tile(width = 0.85, height = 0.85) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_manual(values = c("#1F4E79", "#6FA8DC", "#E8590C", "#C9C9C9"),
                               labels = sprintf("%s (%d)", names(g), f), name = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::labs(title = paste0("What happens to 1,000 people with a Galleri +ve esophagus/stomach signal", dummy_note(res))) +
    ggplot2::theme_void(base_size = 12) + ggplot2::theme(legend.position = "bottom")
}

plot_trace <- function(res, p, data) {
  d <- trace_table(p, data, res$settings, res$analysis$name)
  d$arm <- factor(d$arm, res$arms)
  ggplot2::ggplot(d, ggplot2::aes(cycle, prop, fill = state)) + ggplot2::geom_area(colour = NA) +
    ggplot2::facet_wrap(~arm) + ggplot2::scale_fill_manual(values = STATE_COLS, name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Years after screening", y = "Share of cohort",
                  title = "Markov state occupancy over the lifetime (competing risks)") + ceakit::theme_cea()
}

plot_cancer_states <- function(res, p, data) {
  d <- trace_table(p, data, res$settings, res$analysis$name)
  d <- d[d$state %in% c("Undiagnosed cancer", "Early-stage dx", "Late-stage dx", "Cancer death"), ]
  d$arm <- factor(d$arm, res$arms)
  ggplot2::ggplot(d, ggplot2::aes(cycle, prop * res$analysis$cohort_n, colour = arm)) +
    ggplot2::geom_line(linewidth = 1) + ggplot2::facet_wrap(~state, scales = "free_y") +
    ggplot2::scale_colour_manual(values = ceakit::arm_colours(res$arms), name = NULL) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(x = "Years after screening", y = res$analysis$unit, title = "Cancer states over time by arm") +
    ceakit::theme_cea()
}

plot_subgroups <- function(res, p, data) {
  d <- subgroup_table(p, data, res$settings, res$analysis$name, res$analysis$cohort_n)
  d <- d[d$arm != res$arms[1], ]
  d$subgroup <- sprintf("%s, %s, %s", d$risk_group, d$sex, d$start_age)
  d <- tidyr::pivot_longer(d, c(cancers_found, false_positives, cancer_deaths_averted, net_qalys),
                           names_to = "metric", values_to = "value")
  d$metric <- factor(d$metric, c("cancers_found", "false_positives", "cancer_deaths_averted", "net_qalys"),
                     c("Cancers found", "False positives", "Cancer deaths averted", "Net QALYs"))
  ggplot2::ggplot(d, ggplot2::aes(value, subgroup, fill = factor(arm, res$arms))) +
    ggplot2::geom_col(position = "dodge", width = 0.75) + ggplot2::facet_wrap(~metric, scales = "free_x") +
    ggplot2::scale_fill_manual(values = ceakit::arm_colours(res$arms), name = NULL) +
    ggplot2::labs(x = res$analysis$unit, y = NULL, title = "Results by subgroup (cohort row), vs reference") +
    ceakit::theme_cea()
}
