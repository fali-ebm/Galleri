# Galleri upper GI screening

Decision model and cost-effectiveness analysis built with **ceakit**.

## This project

An annual-cycle cohort model over a lifetime horizon (to age 100), with
three arms and **annual screening from the start age to 85**:

1. **Usual care**: no screening. Cancers present with symptoms.
2. **Kunzmann -> EGD**: the score is applied once; people with a score ≥ 8
   have an EGD every year. Everyone else follows usual care.
3. **Galleri**: a test every year, scored against preclinical upper GI
   cancer as a 2×2 table:
   - **Galleri +ve** (esophagus/stomach signal) → EGD. **EGD +ve** is a true
     positive. **EGD -ve** is treated as a false positive for upper GI
     cancer and goes on to the rest of the workup (ENT or PET CT, then EUS ±
     sampling, then tumor board). Any EGD can **perforate**
     (`p_perf_egd`); a fatal perforation (`p_death_perf`) ends the workup.
   - **Galleri -ve** → usual care. **True negatives** stay well and are
     rescreened next round. **False negatives** keep their undiagnosed
     cancer: they present with symptoms later (reported as
     `fn_symptomatic_1y` for the first year) or are found at a later round.
   - Outputs report true positives, false positives (and how many of those
     had another cancer), true negatives and false negatives, summed over all
     rounds. The Kunzmann arm reports the same table for its screening EGD.

**Overdiagnosis.** A screen-detected cancer is overdiagnosed if, left alone,
it would never have caused symptoms within the lifetime horizon (death from
another cause first). For each age the model computes that probability from
the calibrated natural history (`overdx_prob`), so overdiagnosis is reported
for upper GI and other-site cancers in both screening arms. A test checks it
equals the excess lifetime diagnoses versus usual care exactly.

**Other model features.**
- Galleri sensitivity is lower in repeat rounds than in the first round
  (`sens_g_rep_ratio`, NHS-Galleri).
- Other-cause mortality is cause-deleted: SEER upper GI cancer deaths are
  removed from the NCHS life tables.
- Costs are in 2025 USD (CPI medical care; `data-raw/inflate_costs.R`).
- The Kunzmann flagged share and capture differ by sex (59.7% of men, 3.5%
  of women flagged).
- Incidental PET findings are followed up (`c_incid_fu`), and symptomatic
  diagnoses carry the EGD perforation risk.
- Scenarios: missed cancers stay missed by Galleri (`p_fn_persist`), and the
  Galleri price shared across all cancers it detects (`galleri_cost_share`).

**Editing inputs in Excel.** `Galleri_model_inputs.xlsx` holds every input
(parameters, settings, cohort, scenarios, workup counts, references and the
raw SEER/NCHS/CPI tables). Edit the blue cells, save, then:

```
Rscript import_inputs_xlsx.R --check   # validate and list every change, write nothing
Rscript import_inputs_xlsx.R           # back up, write the CSVs, re-inflate costs, rebuild the baseline
Rscript run.R
```

The importer validates on a scratch copy (types, ranges, keys, complete raw
tables, then a model load and calibration check) before touching anything,
and backs up replaced files to `data-raw/import_backup/`. After editing the
CSVs directly, refresh the workbook with `python3 tools/export_inputs_xlsx.py`.

**False-positive burden (QALYs and DALYs).** Cancer-free people worked up after a
Galleri positive lose the whole-pathway QALYs measured by Matza 2024 (0.079-0.111
by where the workup ends), plus persisting worry (17% for 6 months), complications,
deaths and radiation cancers. DALYs use the GBD anxiety weight over the days of
uncertainty. Every component is reported per 100,000. `accuracy_by_group.csv` gives
the 2×2, PPV and NPV by attained age for each sex × race/ethnicity group.

**Workup counts and costs of false positives.** For cancer-free people worked up after a Galleri
+ve, the model reports EGDs, ENT, PET CT, EUS and sampling, perforations,
complication deaths, PET CT dose and radiation cancers, days in workup,
workup and complication cost, and the QALYs lost. These follow the workup
paths cancer-free people actually take (each path weighted by its
cancer-free share). In the Kunzmann arm, every cancer-free person screened
carries the EGD burden.

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
the stop age or the interval in `inputs/settings.csv` (defaults in `R/model.R`).

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
