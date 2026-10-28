# Galleri decision model and the ceakit toolkit

This repository holds two things:

| Folder | What it is |
|---|---|
| [`ceakit/`](ceakit) | An R package with reusable tools for decision modelling and cost-effectiveness analysis: input tables and validation, decision-tree and Markov engines, PSA, one-way SA, scenarios, ICER frontier, CEAC, EVPI/EVPPI, tables, plots, a Shiny app, and a project template. |
| [`galleri_cea/`](galleri_cea) | The Galleri upper GI model, built with ceakit. It has two analyses: the Galleri +ve workup (no workup vs EGD only vs full workup) and the screening comparison (usual care vs Kunzmann score ≥ 8 → EGD vs Galleri → EGD only vs Galleri → full workup). |

> **Inputs are placeholders.** Every parameter except the Kunzmann score
> capture (29.5% flagged, 77.4% of EAC cases; Kunzmann et al., *Clin
> Gastroenterol Hepatol* 2018) is marked `DUMMY`. Results are illustrative
> until real inputs are added.

## Install the toolkit

```r
install.packages("remotes")
remotes::install_github("fali-ebm/Galleri", subdir = "ceakit")
# optional, for EVPPI via voi and hand-off to other packages:
install.packages(c("voi", "dampack", "BCEA", "targets", "shiny", "bslib", "DT", "markdown"))
```

## Run the Galleri model

```r
setwd("galleri_cea")
source("run.R")          # checks inputs, runs both analyses, tests, report -> outputs/
# or: targets::tar_make()  (re-runs only what changed)
shiny::runApp()          # interactive app
```

## Start a new model

```r
ceakit::new_project("~/projects/my_model", title = "My model")
```

See [`ceakit/README.md`](ceakit/README.md) for the toolkit and the project
contract, and [`galleri_cea/README.md`](galleri_cea/README.md) for the
Galleri model.

## Validation

- `R CMD check` passes with status OK. There are 43 package test
  expectations and 19 Galleri verification tests.
- ICERs, dominance, CEAC and EVPI match `dampack` and `BCEA` exactly.
- EVPPI uses `voi` when it is installed; the built-in fallback agrees with
  `voi` to within 2%.
