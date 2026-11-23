# Visual stories

Two self-contained pages built from the model outputs (D3 7.9 + d3-sankey 0.12.3 from CDNs).

| Page | Built by | Shows |
|---|---|---|
| `galleri_burden.html` | `build_burden.py` + `burden_template.html` | The burden of yearly Galleri for 100,000 men aged 50-85: every test and positive by age, fates of the positives, the workup flow, one false alarm, overdiagnosis by age, false negatives. |
| `galleri_vs_usual_care.html` | `build_compare.py` + `compare_template.html` | Galleri vs usual care: deaths prevented, cumulative deaths, 1,000 PSA runs, price threshold, scenarios, 108 subgroups, the Kunzmann option. |

Rebuild after `Rscript run.R`:

```
Rscript tools/viz_data.R          # viz/galleri_viz_data.json (per-age events, flows, PSA, scenarios, subgroups)
python3 viz/build_burden.py
python3 viz/build_compare.py
```

Colours (validated for colour-vision deficiency on the dark surface): aqua = cancer found / benefit, orange = false alarm / cost, violet = overdiagnosed; a hollow ring marks a missed cancer.
