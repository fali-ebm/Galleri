# Galleri upper GI screening

Decision model and cost-effectiveness analysis built with **ceakit**.

## This project

An annual-cycle cohort model over a lifetime horizon (to age 100), with
three arms and **annual screening from the start age to 85**:

1. **Usual care**: no screening. Cancers present with symptoms.
2. **Kunzmann -> EGD**: the score is applied once; people with a score ≥ 8
   have an EGD every year. Everyone else follows usual care.
3. **Galleri**: a test every year. A -ve result (or a signal pointing
   elsewhere) follows usual care. A +ve esophagus/stomach signal goes to
   the full flowchart workup (EGD, then ENT or PET CT, then EUS ± sampling,
   then tumor board).

**Calibrated to SEER.** For each subgroup (race/ethnicity × sex × Kunzmann
group) the model solves the yearly onset of preclinical upper GI cancer so
that, without screening, symptomatic diagnoses reproduce SEER 2013-2017
age-specific incidence (esophagus + stomach) exactly, and their stage
reproduces SEER's stage distribution (esophagus and stomach separately, by
race and sex, weighted by the subgroup's site mix at each age). That calibrated
natural history is the baseline risk in every arm, so screening can only
move cancers from symptomatic to screen-detected. The fit is checked in
the test suite and written to `outputs/<analysis>/calibration.csv`.

**Base case:** 50-year-old non-Hispanic White men (`inputs/cohort.csv`).
**Sensitivity analysis:** every combination of age at first screen (50-90
by 5 years), sex and race/ethnicity (all; non-Hispanic White, Black,
Asian/Pacific Islander, American Indian/Alaska Native; Hispanic), each
re-calibrated to its own SEER incidence and life table
(`outputs/<analysis>/subgroups.csv` and `subgroups.png`). Change the grid,
the stop age or the interval in `settings` in `R/model.R`.

Two analyses, both per 100,000 people:

- **`galleri`**: usual care vs Galleri, giving the benefits and harms of
  Galleri, the workup tree and what happens to each 1,000 people worked up.
- **`comparison`**: all three arms, with the incremental analysis that
  ranks them.

In the Galleri arm, Galleri accuracy (sensitivity and false-positive rate)
sets how many people test +ve each round and how many of them have cancer.
The workup inputs set the procedures, costs, harms and days for each
person worked up, and the share of those cancers the workup finds.

Screen-detected cancers enter lead-time states (`lead_early`, `lead_late`)
and have no cancer deaths until they would have presented with symptoms.
Stage-specific survival is measured from clinical diagnosis, so screening
is credited only for the stage shift, not for the lead time. The test suite
checks this: detection without a stage shift leaves cancer deaths
unchanged.

**Code:**
- `R/galleri_engine.R`: the workup tree (`ceakit::tree_paths`), SEER
  calibration, the annual screening simulation, the subgroup grid and the derived
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
from the raw published tables in `data-raw/raw/`: SEER 2013-2017
incidence by age, sex and race/ethnicity, and NCHS 2024 life tables by
race/ethnicity and sex. `inputs/SOURCES.md` lists the source or reasoning
for every input. Inputs are marked `SOURCED` or `ASSUMPTION`. The workup
branch probabilities come from the literature until local counts are
added to `inputs/counts.csv`.

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
