# Visual pages

Three self-contained pages built in R from the model outputs. The two story pages use
D3 7.9 and d3-sankey 0.12.3 from CDNs; the decision-model page is static SVG.

| Page | Built by | Shows |
|---|---|---|
| `galleri_model.html` | `build_model_page.R` + `model_page_template.html` + `model_page.css` | The decision model step by step: the three arms, the EGD step (hold, two kinds of miss), ENT/PET CT, pancreatobiliary first and repeat findings, natural history, cohort vs microsimulation, scenarios, inputs to review. Every edge label is read from `inputs/parameters.csv`. |
| `galleri_burden.html` | `build_burden.R` + `burden_template.html` | The burden of yearly Galleri for 100,000 men aged 50-85: every test and positive by age, fates of the positives, the workup flow, one false alarm, overdiagnosis by age, false negatives. |
| `galleri_vs_usual_care.html` | `build_compare.R` + `compare_template.html` | Galleri vs usual care: deaths prevented, cumulative deaths, 1,000 PSA runs, price threshold, scenarios, 108 subgroups, the microsimulation person by person (who gains, who loses, cohort vs microsim), the MISCAN-ESO recalibration (natural history, fit by age, CISNET validation, scenarios), the Kunzmann option. Needs `outputs/microsim` and `outputs/miscan_eac`. |

Rebuild (from `galleri_cea/`, after `Rscript run.R` and `Rscript tools/run_microsim.R`), or run
everything at once with `Rscript make_all.R`:

```
Rscript tools/viz_data.R          # viz/galleri_viz_data.json (per-age events, flows, PSA, scenarios, subgroups)
Rscript viz/build_burden.R
Rscript viz/build_compare.R
Rscript viz/build_model_page.R
```

Colours (validated for colour-vision deficiency on the dark surface): aqua = cancer found / benefit, orange = false alarm / cost, violet = overdiagnosed; a hollow ring marks a missed cancer.
