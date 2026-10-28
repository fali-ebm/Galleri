#' Evaluate a decision tree
#'
#' @param edges Data frame with `from`, `to`, `prob`. The root is the node that
#'   never appears in `to`. Children of every node must sum to 1.
#' @param nodes Optional data frame with `id` and numeric attribute columns
#'   (e.g. procedure counts, costs). Attributes are summed along each path.
#' @param tol Tolerance for the sum-to-one check.
#' @return Tibble with one row per leaf: `leaf`, `prob`, `path` and the summed
#'   node attributes.
#' @export
tree_paths <- function(edges, nodes = NULL, tol = 1e-8) {
  edges <- as.data.frame(edges)
  check_tree(edges, tol)
  root <- setdiff(edges$from, edges$to)
  if (length(root) != 1) stop("tree must have exactly one root; found: ", paste(root, collapse = ", "), call. = FALSE)
  attrs <- if (is.null(nodes)) character(0) else setdiff(names(nodes), "id")
  node_val <- function(id) {
    if (!length(attrs)) return(numeric(0))
    r <- nodes[nodes$id == id, attrs, drop = FALSE]
    if (!nrow(r)) stats::setNames(rep(0, length(attrs)), attrs) else unlist(r[1, ])
  }
  out <- list()
  walk <- function(id, p, path, acc) {
    kids <- edges[edges$from == id, ]
    if (!nrow(kids)) { out[[length(out) + 1]] <<- c(list(leaf = id, prob = p, path = path), as.list(acc)); return(invisible()) }
    for (i in seq_len(nrow(kids))) {
      to <- kids$to[i]
      walk(to, p * kids$prob[i], paste(path, to, sep = " > "), acc + node_val(to))
    }
  }
  walk(root, 1, root, node_val(root))
  dplyr::bind_rows(lapply(out, tibble::as_tibble))
}

#' Check a decision tree's branch probabilities
#' @inheritParams tree_paths
#' @export
check_tree <- function(edges, tol = 1e-8) {
  if (any(is.na(edges$prob) | edges$prob < -tol | edges$prob > 1 + tol))
    stop("tree: branch probabilities must be between 0 and 1", call. = FALSE)
  s <- tapply(edges$prob, edges$from, sum)
  bad <- names(s)[abs(s - 1) > tol]
  if (length(bad)) stop("tree: probabilities leaving ", paste(bad, collapse = ", "),
                        " do not sum to 1 (", paste(round(s[bad], 6), collapse = ", "), ")", call. = FALSE)
  invisible(TRUE)
}

#' Expected lifetime values of a cohort Markov model (backward induction)
#'
#' Computes, for every starting state, the expected (discounted) sum of
#' rewards over `n_cycles`. Because results are linear in the starting
#' distribution, outcomes for any arm are `init %*% values`.
#'
#' @param n_cycles Number of cycles.
#' @param transition `function(t)` returning an S x S transition matrix for
#'   cycle `t` (rows = from, columns = to; rows must sum to 1).
#' @param rewards `function(t, P)` returning an S x K matrix of rewards earned
#'   in cycle `t` by someone starting the cycle in each state (include
#'   half-cycle correction and transition-triggered rewards here).
#' @param disc Named numeric vector of annual discount rates, one per reward
#'   column (missing columns are undiscounted).
#' @param cycle_length Cycle length in years (for discounting).
#' @return S x K matrix of expected values from each starting state.
#' @export
markov_values <- function(n_cycles, transition, rewards, disc = NULL, cycle_length = 1) {
  V <- NULL
  for (t in rev(seq_len(n_cycles))) {
    P <- transition(t)
    check_transition(P, t)
    R <- rewards(t, P)
    if (is.null(V)) {
      V <- R * 0
      d <- rep(1, ncol(R)); names(d) <- colnames(R)
      if (!is.null(disc)) for (k in intersect(names(disc), colnames(R))) d[k] <- 1 / (1 + disc[[k]])^cycle_length
    }
    V <- R + sweep(P %*% V, 2, d, "*")
  }
  V
}

#' Forward state occupancy of a cohort Markov model
#' @inheritParams markov_values
#' @param init Starting distribution over states (sums to 1).
#' @return Tibble with `cycle` (0 = start) and one column per state.
#' @export
markov_trace <- function(n_cycles, transition, init) {
  S <- matrix(0, n_cycles + 1, length(init), dimnames = list(NULL, names(init)))
  S[1, ] <- init
  for (t in seq_len(n_cycles)) {
    P <- transition(t); check_transition(P, t)
    S[t + 1, ] <- S[t, ] %*% P
  }
  tibble::as_tibble(S) |> dplyr::mutate(cycle = 0:n_cycles, .before = 1)
}

#' Check a transition matrix
#' @param P Transition matrix.
#' @param t Cycle (for messages).
#' @param tol Tolerance.
#' @export
check_transition <- function(P, t = NA, tol = 1e-8) {
  if (any(P < -tol | P > 1 + tol, na.rm = TRUE) || anyNA(P))
    stop("transition matrix (cycle ", t, ") has values outside 0-1 or NA", call. = FALSE)
  rs <- rowSums(P)
  if (any(abs(rs - 1) > tol))
    stop("transition matrix (cycle ", t, ") rows do not sum to 1: ",
         paste(rownames(P)[abs(rs - 1) > tol], collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

#' Half-cycle-corrected share of a cycle spent alive
#' @param P Transition matrix.
#' @param dead Names of dead states.
#' @return Numeric vector, one per state (0 for dead states).
#' @export
alive_fraction <- function(P, dead) {
  die <- rowSums(P[, dead, drop = FALSE])
  ifelse(rownames(P) %in% dead, 0, 1 - 0.5 * die)
}

#' Check that a Markov trace conserves the cohort
#' @param trace Output of [markov_trace()].
#' @param tol Tolerance.
#' @export
check_trace <- function(trace, tol = 1e-8) {
  s <- rowSums(as.matrix(trace[, setdiff(names(trace), "cycle")]))
  if (any(abs(s - 1) > tol)) stop("trace does not sum to 1 in every cycle", call. = FALSE)
  invisible(TRUE)
}
