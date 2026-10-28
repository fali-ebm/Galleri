# =============================================================
# PROJECT MODEL - the only file with project-specific logic.
#
# Contract (read by ceakit::load_project):
#   analyses        list of ceakit::cea_analysis() objects       (required)
#   settings        ceakit::default_settings(...)                  (optional)
#   load_data(dir)  returns named list of data tables              (optional;
#                   default = every table in inputs/ except parameters,
#                   counts and scenarios)
#   validate_data(data, settings) returns list(errors, warnings)  (optional)
#   project_title   shown in the app                               (optional)
#
# EXAMPLE MODEL (replace): one-off screening vs usual care.
#   Time 0 decision tree: test -> positive (confirmatory test) / negative.
#   Lifetime Markov: well, undiagnosed disease, treated disease, dead.
#   Competing risk: age-specific background mortality from inputs/life_table.csv.
# =============================================================

project_title <- "{{title}}"

settings <- ceakit::default_settings(
  wtp = 50000, n_sim = 1000, seed = 2026,
  start_age = 60, max_age = 100,            # project-specific settings
  disc_cost = 0.03, disc_qaly = 0.03)

STATES <- c("well", "undiagnosed", "treated", "dead")
ARMS <- c("Usual care", "Screening")

validate_data <- function(data, settings) {
  e <- character(0)
  lt <- data$life_table
  if (is.null(lt)) return(list(errors = "inputs/life_table.csv is missing", warnings = character(0)))
  if (any(lt$mort < 0 | lt$mort > 1)) e <- c(e, "life_table: mort must be between 0 and 1")
  ok <- tryCatch({ ceakit::lookup_by_age(lt, list(sex = "all"), settings$start_age:(settings$max_age - 1)); TRUE },
                 error = function(err) { e <<- c(e, paste("life_table:", conditionMessage(err))); FALSE })
  list(errors = e, warnings = character(0))
}

# ---- Decision tree at time 0 (screening arm) ------------------------------
screen_tree <- function(p) {
  edges <- data.frame(
    from = c("screen", "screen", "diseased", "diseased", "healthy", "healthy"),
    to   = c("diseased", "healthy", "true_pos", "false_neg", "false_pos", "true_neg"),
    prob = c(p$p_prev, 1 - p$p_prev, p$sens, 1 - p$sens, 1 - p$spec, p$spec))
  nodes <- data.frame(id = c("screen", "true_pos", "false_pos"),
                      cost = c(p$c_test, p$c_confirm, p$c_confirm),
                      confirm = c(0, 1, 1))
  ceakit::tree_paths(edges, nodes)
}

# ---- Markov model ----------------------------------------------------------
markov_parts <- function(p, data, settings) {
  ages <- settings$start_age:(settings$max_age - 1)
  lt <- ceakit::lookup_by_age(data$life_table, list(sex = "all"), ages)
  transition <- function(t) {
    q <- lt$mort[t]; s <- 1 - q
    P <- matrix(0, 4, 4, dimnames = list(STATES, STATES))
    P["well", ] <- c(s, 0, 0, q)
    P["undiagnosed", ] <- c(0, s * (1 - p$m_untreated) * (1 - p$p_dx_clin),
                            s * (1 - p$m_untreated) * p$p_dx_clin, q + s * p$m_untreated)
    P["treated", ] <- c(0, 0, s * (1 - p$m_treated), q + s * p$m_treated)
    P["dead", "dead"] <- 1
    P
  }
  rewards <- function(t, P) {
    alive <- ceakit::alive_fraction(P, "dead")
    cbind(qaly = alive * c(p$u_well, p$u_untreated, p$u_treated, 0),
          cost = alive * c(0, 0, p$c_tx_annual, 0) +
                 c(0, P["undiagnosed", "treated"] * p$c_dx_clin, 0, 0),
          ly = alive)
  }
  list(n = length(ages), transition = transition, rewards = rewards)
}

run_base <- function(p, data, settings) {
  mk <- markov_parts(p, data, settings)
  V <- ceakit::markov_values(mk$n, mk$transition, mk$rewards,
                             disc = c(cost = settings$disc_cost, qaly = settings$disc_qaly))
  tree <- screen_tree(p)
  tp <- tree$prob[tree$leaf == "true_pos"]; fn <- tree$prob[tree$leaf == "false_neg"]
  fp <- tree$prob[tree$leaf == "false_pos"]
  init <- rbind("Usual care" = c(1 - p$p_prev, p$p_prev, 0, 0),
                "Screening"  = c(1 - p$p_prev, fn, tp, 0))
  colnames(init) <- STATES
  life <- init %*% V
  cost0 <- c(0, sum(tree$prob * tree$cost))
  qaly0 <- c(0, fp * p$disu_fp)                       # anxiety/harm of false positives
  m <- rbind(cost = cost0 + life[, "cost"], qaly = life[, "qaly"] - qaly0, ly = life[, "ly"],
             detected = c(0, tp), false_pos = c(0, fp), tests = c(0, 1))
  colnames(m) <- ARMS
  m
}

analyses <- list(
  base = ceakit::cea_analysis(
    name = "base", title = "Screening vs usual care", arms = ARMS, run = run_base,
    cost = "cost", effect = "qaly", effect_label = "QALY",
    required = c("p_prev", "sens", "spec", "c_test", "c_confirm", "p_dx_clin", "m_untreated",
                 "m_treated", "c_tx_annual", "c_dx_clin", "u_well", "u_untreated", "u_treated", "disu_fp"),
    cohort_n = 1000, unit = "per 1,000 screened",
    metrics = tibble::tribble(
      ~key, ~section, ~label, ~type, ~digits,
      "cost", "Outcomes", "Lifetime cost per person ($)", "person", 2,
      "qaly", "Outcomes", "QALYs", "count", 1,
      "ly", "Outcomes", "Life years", "count", 1,
      "tests", "Screening", "Screening tests", "count", 0,
      "detected", "Screening", "Cases detected by screening", "count", 1,
      "false_pos", "Screening", "False positives", "count", 1),
    key_metrics = c("detected", "false_pos", "qaly", "inc_nmb"),
    plots = list(
      decision_tree = function(res, p, data) {
        t <- screen_tree(p)
        nodes <- data.frame(id = c("screen", "diseased", "healthy", t$leaf),
          label = c("Screened", "Has disease", "No disease", "True positive", "False negative",
                    "False positive", "True negative"),
          x = c(5, 2.5, 7.5, 1.3, 3.7, 6.3, 8.7), y = c(3, 2, 2, 1, 1, 1, 1), w = 2)
        edges <- data.frame(from = c("screen", "screen", "diseased", "diseased", "healthy", "healthy"),
                            to = c("diseased", "healthy", t$leaf),
                            prob = c(p$p_prev, 1 - p$p_prev, p$sens, 1 - p$sens, 1 - p$spec, p$spec))
        ceakit::draw_tree(nodes, edges, "Screening decision tree", xlim = c(0, 10), ylim = c(0.4, 3.6))
      },
      markov_trace = function(res, p, data) {
        mk <- markov_parts(p, data, settings)
        tr <- ceakit::markov_trace(mk$n, mk$transition, c(well = 1 - p$p_prev, undiagnosed = p$p_prev,
                                                         treated = 0, dead = 0))
        d <- tidyr::pivot_longer(tr, -cycle, names_to = "state", values_to = "prop")
        ggplot2::ggplot(d, ggplot2::aes(cycle, prop, fill = state)) + ggplot2::geom_area() +
          ggplot2::labs(x = "Years", y = "Share of cohort", title = "Markov trace (usual care)") +
          ceakit::theme_cea()
      }))
)
