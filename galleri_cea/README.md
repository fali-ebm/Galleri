# Galleri upper GI screening

Decision model and cost-effectiveness analysis built with **ceakit**.

## This project

Two analyses, both a decision tree at time 0 followed by a lifetime Markov
model with competing risks:

- **`galleri_positive`**: people with a Galleri +ve esophagus/stomach
  signal. It compares no workup, EGD only and the full flowchart workup
  (EGD, then ENT or PET CT, then EUS ± sampling, then tumor board). The
  results are reported per 1,000 Galleri +ve.
- **`screening`**: a screening population. It compares usual care, a
  Kunzmann score ≥ 8 followed by EGD, Galleri followed by EGD only, and
  Galleri followed by the full workup. The results are reported per
  100,000 screened.

Screen-detected cancers enter lead-time states (`lead_early`, `lead_late`)
and have no cancer deaths until they would have presented with symptoms.
Stage-specific survival is measured from clinical diagnosis, so screening
is credited only for the stage shift, not for the lead time. The test suite
checks this: detection without a stage shift leaves cancer deaths
unchanged.

**Code:**
- `R/galleri_engine.R`: the tree (`ceakit::tree_paths`), the Markov model
  (`ceakit::markov_values`), the time-0 event ledgers and the derived
  metrics.
- `R/galleri_plots.R`: project figures.
- `R/model.R`: the ceakit contract (settings, data loading and validation,
  metric definitions, analyses).

**Inputs:**
- `parameters.csv`: all inputs, with their distributions.
- `baseline_risks.csv`: competing and baseline risks by risk group × sex ×
  age band (other-cause mortality, incidence, prevalence, utility,
  reference life expectancy for DALYs).
- `cohort.csv`: the cohort mix, including which rows the Kunzmann arm
  scopes.
- `counts.csv`: optional observed workup counts.
- `scenarios.csv`: deterministic scenarios.

`data-raw/make_us_inputs.R` builds `baseline_risks.csv` and `cohort.csv`
from published US data: the SSA 2022 life table, SEER incidence, Census 2024
population and the Kunzmann score capture. `inputs/SOURCES.md` lists the
source or reasoning for every input. Inputs are marked `SOURCED`,
`ASSUMPTION` or `DUMMY`. The workup branch probabilities are still `DUMMY`
until local counts are added to `inputs/counts.csv`.

## Workflow

1. **Inputs**: edit `inputs/parameters.csv`. Each row has a value, a PSA
   distribution (`dist` with `se` or `n`), and a `status` and `source`.
   - Change `DUMMY` to `SOURCED`, `ASSUMPTION` or `DATA` as evidence comes in.
   - Other tables in `inputs/` (e.g. `life_table.csv`) are data the model
     reads.
   - `inputs/scenarios.csv` lists deterministic scenarios.
   - Optional `inputs/counts.csv` (columns `branch, outcome, parameter, n`)
     turns observed counts into probabilities.
2. **Model**: `R/model.R` is the only project-specific code. It defines
   `analyses` (arms plus a `run()` function), `settings`, and
   `validate_data()`.
3. **Run**: `Rscript run.R`, or `targets::tar_make()`. Each run:
   - checks the inputs;
   - runs the base case, PSA, incremental analysis, CEAC, EVPI, EVPPI,
     one-way SA and scenarios;
   - writes everything to `outputs/`;
   - runs the verification tests;
   - renders `outputs/report.html` (CHEERS-structured).
4. **Explore or share**: `shiny::runApp()` opens the interactive app.

## Reproducibility

- `git init` and commit `inputs/`, `R/`, `report.Rmd` and `tests/`.
  `outputs/` is ignored.
- `renv::init()` once to lock package versions, and `renv::snapshot()`
  after updating packages.
- Add a test in `tests/testthat/` for every bug you fix.

## Other packages

To use other packages, call `ceakit::psa_matrices(res)`. It returns cost,
effect and parameter matrices for `dampack::make_psa_obj()`,
`BCEA::bcea()` and the `voi` package.
