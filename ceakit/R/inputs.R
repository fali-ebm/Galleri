#' Read a CSV or Excel table
#'
#' @param path File path (`.csv`, `.xlsx` or `.xls`).
#' @param name Optional original file name (used to detect the type of
#'   uploaded files whose temporary path has no extension).
#' @return A tibble.
#' @export
read_table <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) {
    if (tools::file_ext(path) == "") { tmp <- paste0(path, ".", ext); file.copy(path, tmp, TRUE); path <- tmp }
    return(tibble::as_tibble(readxl::read_excel(path)))
  }
  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA"), guess_max = 10000)
}

PARAM_REQUIRED_COLS <- c("parameter", "value", "dist")
PARAM_OPTIONAL_COLS <- c("label", "status", "se", "n", "group", "low", "high", "unit",
                         "description", "source")
DISTRIBUTIONS <- c("fixed", "beta", "gamma", "normal", "lognormal", "dirichlet")

#' Tidy a parameter table
#'
#' Coerces column types and adds any missing optional columns. The parameter
#' table has one row per model input: `parameter`, `value`, `dist`, and
#' optionally `label`, `status` (DUMMY / ASSUMPTION / SOURCED / DATA), `se`, `n`, `group`
#' (dirichlet groups), `low`, `high` (one-way SA range), `unit`, `description`,
#' `source`.
#' @param df Data frame.
#' @return A tibble.
#' @export
tidy_parameters <- function(df) {
  df <- tibble::as_tibble(df)
  for (col in setdiff(PARAM_OPTIONAL_COLS, names(df))) df[[col]] <- NA
  num <- c("value", "se", "n", "low", "high")
  chr <- c("parameter", "label", "status", "dist", "group", "unit", "description", "source")
  df[num] <- lapply(df[num], function(x) suppressWarnings(as.numeric(x)))
  df[chr] <- lapply(df[chr], as.character)
  df$dist <- tolower(trimws(df$dist))
  df$status <- toupper(trimws(df$status))
  df
}

#' Read a parameter table from CSV or Excel
#' @param path File path.
#' @return A tidy parameter tibble.
#' @export
read_parameters <- function(path) tidy_parameters(read_table(path))

#' Replace probabilities with estimates from observed counts
#'
#' `counts` has columns `branch`, `outcome`, `parameter`, `n`. Each branch lists
#' all of its mutually exclusive outcomes; rows naming a `parameter` set that
#' parameter. A branch is used only when every outcome has a count. Estimates
#' use a Jeffreys prior (+0.5 per outcome) and become Beta (2 outcomes) or
#' Dirichlet (3+ outcomes) distributions sized by the observed total.
#' @param params Parameter table.
#' @param counts Counts table (or `NULL`).
#' @return Updated parameter table.
#' @export
apply_counts <- function(params, counts) {
  if (is.null(counts) || !nrow(counts)) return(params)
  counts$n <- suppressWarnings(as.numeric(counts$n))
  if (all(is.na(counts$n))) return(params)
  for (b in unique(counts$branch)) {
    rows <- counts[counts$branch == b, ]
    if (all(is.na(rows$n))) next
    if (any(is.na(rows$n))) stop("counts: branch '", b, "' is partly filled; fill every outcome or none.", call. = FALSE)
    if (sum(rows$n) <= 0 || any(rows$n < 0)) stop("counts: branch '", b, "' needs non-negative counts with a positive total.", call. = FALSE)
    k <- nrow(rows); tot <- sum(rows$n) + 0.5 * k
    for (i in which(!is.na(rows$parameter) & rows$parameter != "")) {
      j <- match(rows$parameter[i], params$parameter)
      if (is.na(j)) stop("counts: unknown parameter ", rows$parameter[i], call. = FALSE)
      params$value[j] <- (rows$n[i] + 0.5) / tot
      params$n[j] <- tot; params$se[j] <- NA
      params$dist[j] <- if (k == 2) "beta" else "dirichlet"
      params$group[j] <- if (k == 2) NA else paste0("counts_", b)
      params$status[j] <- "DATA"
      params$source[j] <- sprintf("observed counts: %g of %g", rows$n[i], sum(rows$n))
    }
  }
  params
}

#' Validate a parameter table
#'
#' @param params Parameter table.
#' @param required Character vector of parameter names the model needs.
#' @param probability Regular expression for parameter names that must lie
#'   in 0-1 (in addition to beta / dirichlet parameters).
#' @return `list(errors, warnings)`; never stops.
#' @export
validate_parameters <- function(params, required = NULL,
                                probability = "^(p_|prob_|sens_|spec_|fpr_)") {
  e <- character(0); w <- character(0); P <- params
  miss <- setdiff(PARAM_REQUIRED_COLS, names(P))
  if (length(miss)) return(list(errors = paste("parameters: missing columns", paste(miss, collapse = ", ")),
                                warnings = w))
  m <- setdiff(required, P$parameter)
  if (length(m)) e <- c(e, paste("parameters: missing rows", paste(m, collapse = ", ")))
  d <- unique(P$parameter[duplicated(P$parameter)])
  if (length(d)) e <- c(e, paste("parameters: duplicated", paste(d, collapse = ", ")))
  na <- P$parameter[is.na(P$value)]
  if (length(na)) e <- c(e, paste("parameters: no value for", paste(na, collapse = ", ")))
  bd <- P$parameter[!P$dist %in% DISTRIBUTIONS]
  if (length(bd)) e <- c(e, paste("parameters: unknown dist for", paste(bd, collapse = ", "),
                                  "- use", paste(DISTRIBUTIONS, collapse = ", ")))
  is_p <- P$dist %in% c("beta", "dirichlet") | grepl(probability, P$parameter)
  bp <- P$parameter[is_p & !is.na(P$value) & (P$value < 0 | P$value > 1)]
  if (length(bp)) e <- c(e, paste("parameters: outside 0-1:", paste(bp, collapse = ", ")))
  ns <- P$parameter[P$dist %in% c("gamma", "normal", "lognormal") & is.na(P$se)]
  if (length(ns)) e <- c(e, paste("parameters: se needed for", paste(ns, collapse = ", ")))
  nb <- P$parameter[P$dist == "beta" & is.na(P$se) & is.na(P$n)]
  if (length(nb)) e <- c(e, paste("parameters: se or n needed for", paste(nb, collapse = ", ")))
  neg <- P$parameter[P$dist %in% c("gamma", "lognormal") & !is.na(P$value) & P$value < 0]
  if (length(neg)) e <- c(e, paste("parameters: gamma/lognormal must be >= 0:", paste(neg, collapse = ", ")))
  dir <- P[P$dist == "dirichlet", ]
  if (nrow(dir)) {
    if (any(is.na(dir$group) | is.na(dir$n))) e <- c(e, "parameters: every dirichlet parameter needs group and n")
    for (g in unique(stats::na.omit(dir$group))) {
      s <- sum(dir$value[dir$group %in% g])
      if (s > 1 + 1e-9) e <- c(e, sprintf("parameters: dirichlet group '%s' sums to %.3f (> 1)", g, s))
    }
  }
  dummy <- P$parameter[!is.na(P$status) & P$status == "DUMMY"]
  if (length(dummy)) w <- c(w, sprintf("%d parameters are still DUMMY placeholders: %s",
                                       length(dummy), paste(dummy, collapse = ", ")))
  list(errors = e, warnings = w)
}

#' Convert a parameter table to a named list of values
#' @param params Parameter table.
#' @export
param_list <- function(params) stats::setNames(as.list(params$value), params$parameter)

draw_one <- function(dist, m, se, n) {
  if (dist %in% c("fixed", "dirichlet") || is.na(m) || m == 0 || (dist == "beta" && m >= 1)) return(m)
  switch(dist,
    beta = {
      ab <- if (!is.na(n)) c(m * n, (1 - m) * n) else {
        k <- m * (1 - m) / se^2 - 1
        if (k <= 0) stop("beta se too large for mean ", m, call. = FALSE)
        c(m * k, (1 - m) * k) }
      stats::rbeta(1, ab[1], ab[2]) },
    gamma     = stats::rgamma(1, shape = (m / se)^2, scale = se^2 / m),
    normal    = stats::rnorm(1, m, se),
    lognormal = { s2 <- log(1 + (se / m)^2); stats::rlnorm(1, log(m) - s2 / 2, sqrt(s2)) })
}

#' Draw one set of parameter values for the PSA
#'
#' Beta: `n` (effective sample size) or `se`. Gamma, normal, lognormal: `se`.
#' Dirichlet: parameters sharing a `group` and `n`; the remainder of the group
#' is the implied last outcome.
#' @param params Parameter table.
#' @param constraints Optional function(list) -> list applied after sampling
#'   (e.g. to keep utilities ordered).
#' @return Named list of sampled values.
#' @export
sample_parameters <- function(params, constraints = NULL) {
  q <- param_list(params)
  for (i in seq_len(nrow(params)))
    q[[params$parameter[i]]] <- draw_one(params$dist[i], params$value[i], params$se[i], params$n[i])
  for (g in unique(stats::na.omit(params$group[params$dist == "dirichlet"]))) {
    r <- params[params$dist == "dirichlet" & params$group %in% g, ]
    m <- c(r$value, max(0, 1 - sum(r$value)))
    x <- stats::rgamma(length(m), shape = pmax(m * r$n[1], 1e-9)); x <- x / sum(x)
    q[r$parameter] <- as.list(x[seq_len(nrow(r))])
  }
  if (!is.null(constraints)) q <- constraints(q)
  q
}

#' Ranges for one-way sensitivity analysis
#' @param params Parameter table.
#' @param pct Fractional range used when `low` / `high` are blank.
#' @export
owsa_ranges <- function(params, pct = 0.25) {
  v <- params[params$dist != "fixed", ]
  hi <- ifelse(is.na(v$high), v$value * (1 + pct), v$high)
  hi <- ifelse(v$dist %in% c("beta", "dirichlet"), pmin(hi, 1), hi)
  tibble::tibble(parameter = v$parameter,
                 low = ifelse(is.na(v$low), v$value * (1 - pct), v$low), high = hi)
}

#' Look up age-specific values from a stratified table
#'
#' For tables with `age_from`, `age_to` and stratifying columns (e.g.
#' `risk_group`, `sex`). Each stratum value falls back to `fallback`
#' (default `"all"`) when no specific rows exist.
#' @param table Data frame with `age_from`, `age_to` and strata columns.
#' @param strata Named list, e.g. `list(risk_group = "high", sex = "male")`.
#' @param ages Integer vector of ages.
#' @param fallback Value used as wildcard.
#' @return Tibble with `age` and the remaining columns, one row per age.
#' @export
lookup_by_age <- function(table, strata, ages, fallback = "all") {
  keys <- names(strata)
  opts <- lapply(strata, function(v) unique(c(v, fallback)))
  combos <- expand.grid(opts, stringsAsFactors = FALSE)
  combos <- combos[order(rowSums(combos == fallback)), , drop = FALSE]   # most specific first
  rows <- NULL
  for (i in seq_len(nrow(combos))) {
    keep <- Reduce(`&`, lapply(keys, function(k) table[[k]] == combos[i, k]))
    if (any(keep)) { rows <- table[keep, ]; break }
  }
  if (is.null(rows)) stop("no rows for ", paste(keys, unlist(strata), sep = " = ", collapse = ", "),
                          call. = FALSE)
  idx <- vapply(ages, function(a) { k <- which(rows$age_from <= a & rows$age_to >= a); if (length(k)) k[1] else NA_integer_ },
                integer(1))
  if (anyNA(idx)) stop(sprintf("rows for %s do not cover ages %s",
                               paste(unlist(strata), collapse = "/"),
                               paste(range(ages[is.na(idx)]), collapse = "-")), call. = FALSE)
  tibble::tibble(age = ages) |>
    dplyr::bind_cols(rows[idx, setdiff(names(rows), c(keys, "age_from", "age_to")), drop = FALSE])
}
