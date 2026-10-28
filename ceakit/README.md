# ceakit

A toolkit for decision modelling and cost-effectiveness analysis (CEA) in R.
Every project uses the same input format, engines, uncertainty analyses,
outputs, report and app. Only the model structure changes between projects.

## Install

```r
# from a local file
install.packages("ceakit_0.1.0.tar.gz", repos = NULL, type = "source")
# or from GitHub
remotes::install_github("fali-ebm/Galleri", subdir = "ceakit")
```

## Start a new project

```r
ceakit::new_project("~/projects/my_cea", title = "My CEA")
```

Then fill in `inputs/parameters.csv`, write `R/model.R` and run `run.R`.
See the project README for details.

## What's in the box

| Area | Functions |
|---|---|
| Inputs | `read_parameters`, `tidy_parameters`, `validate_parameters`, `apply_counts`, `sample_parameters`, `owsa_ranges`, `lookup_by_age` (stratified baseline/competing risks with `all` fallbacks) |
| Engines | `tree_paths` (decision trees, attributes summed along paths), `markov_values` (backward induction, any number of reward channels, per-channel discounting), `markov_trace`, `alive_fraction` (half-cycle), `check_tree`, `check_transition`, `check_trace` |
| Analysis | `cea_analysis`, `default_settings`, `evaluate`, `run_cea` (base, PSA, frontier, CEAC, EVPI, EVPPI, OWSA), `run_owsa`, `run_scenarios` |
| CEA methods | `icer_frontier` (strong + extended dominance), `ceac`, `evpi`, `evppi` (uses voi when installed, built-in GAM otherwise; Strong et al. 2014) |
| Other packages | `to_dampack`, `to_bcea`, `to_voi`, `psa_matrices`. Package tests check that ICERs, CEAC and EVPI match dampack and BCEA exactly, and that EVPPI agrees with voi |
| Outputs | `results_wide`, `results_long`, `write_cea` (CSV, Excel, PNG), `plot_ce_plane`, `plot_ceac`, `plot_evpi`, `plot_evppi`, `plot_tornado`, `plot_metrics`, `draw_tree` |
| Projects | `new_project`, `load_project`, `check_project`, `run_project`, `cea_app`, `run_app` |

## The project contract

A project's `R/model.R` defines `analyses`: a list of
`cea_analysis(name, arms, run, ...)`. Each `run(p, data, settings)` returns
a matrix of per-person outcomes:
- one column per arm;
- one row per metric, including `cost` and the effect (e.g. `qaly`).

Everything else is handled by the toolkit: sampling, PSA, incremental
analysis, VOI, sensitivity analyses, tables, figures, report and app.
