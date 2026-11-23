# Input sources (US, payer perspective, adults 50-79)

The `status` column in each table means:

- **SOURCED**: taken from, or derived arithmetically from, the source cited.
- **ASSUMPTION**: no direct source. The value and the reasoning are given so a
  reviewer can challenge it.
- **DUMMY**: placeholder with no source. None are left. The workup branch
  probabilities now come from the literature (below) until local counts
  are entered in `inputs/counts.csv`, which then overrides them.

`baseline_risks.csv` and `cohort.csv` are built by
`data-raw/make_us_inputs.R`. Edit the script, not the CSVs.

## Galleri test performance (`parameters.csv`)

| Parameter | Value | Derivation |
|---|---|---|
| `sens_g_ugi` | 0.637 (range 0.40-0.85) | PATHFINDER 2 episode sensitivity for the 12 deadly cancers (69.8%) × cancer signal origin (CSO) accuracy (91.3%). The range runs from CCGA3 stage I-II sensitivity (esophagus 12/25, stomach 4/12) to CCGA3 all-stage sensitivity (esophagus 84/99, stomach 20/29). **ASSUMPTION**: the 12-cancer figure applies to esophagus and stomach. |
| `sens_g_other` | 0.0068 | Overall episode sensitivity 39.3% × CSO error 8.7% × an assumed 20% of wrong CSOs that are called esophagus/stomach. **ASSUMPTION** |
| `fpr_g_es` | 0.00036 | False-positive rate 0.36% (specificity 99.64%) × an assumed 10% of false positives with an esophagus/stomach CSO. **ASSUMPTION** |
| Workup yield | from your counts | Galleri accuracy sets how many people test +ve and how many of them have cancer. Your workup counts set the procedures per person worked up and the share of those cancers the workup finds: found ÷ (found + missed despite workup). |
| `p_late_screen` | 0.47 | 53% of PATHFINDER 2 screen-detected cancers were stage I-II. |

## Galleri +ve workup (literature values until local counts are added)

Searched in PubMed and Consensus in November 2026. No study reports a
step-by-step workup for an esophagus/stomach cancer signal, so each branch
is derived from the closest published anchor. Your `counts.csv` replaces
all of these.

| Parameter | Value | Derivation |
|---|---|---|
| `p_egd_pos` | 0.50 | PPV of the stomach/esophagus signal 64.6% (NHS-Galleri, 3 rounds) × first signal-origin prediction correct 87% (real-world, 111,080 tests) × EGD sensitivity 88.7%. Cross-check: the model's own Galleri accuracy inputs give a PPV of 67.6% among E/S positives, close to the 64.6% observed. |
| `p_ent_route` | 0.30 | **ASSUMPTION**: a local pathway choice, so there are no published data. |
| `p_ent_pos` | 0.10 | **ASSUMPTION**: only case series exist (3 early oropharyngeal cancers found after a Galleri head and neck signal). |
| `p_pet_pb`, `p_pet_incid` | 0.087, 0.159 | **ASSUMPTION**, anchored on PATHFINDER 2: PET/CT after a negative targeted workup found cancer in 7/63 (11.1%). Half of that yield is assumed to come from pancreatobiliary findings (70% malignant × EUS sensitivity 91%) and half from incidental findings (one-third malignant). |
| `p_eus_pos`, `p_eus_neg` | 0.64, 0.26 | **ASSUMPTION**: 70% of PET pancreatobiliary findings malignant × EUS-guided tissue acquisition sensitivity 91% (newer-generation needles 91-94%); 90% of EUS include sampling. |
| `p_miss` | 0.034 | PATHFINDER 2: 4 of about 117 positives with no cancer at diagnostic resolution had cancer diagnosed within 12 months. This matches the study's ratio of cancers found to missed among positives (173 vs 4). Scenario: 17.8% of NHS-Galleri round-1 false positives had cancer diagnosed over the next 2 rounds (includes new cancers). |
| `p_incid_cancer` | 0.33 | Incidental FDG PET/CT findings are malignant about one-third of the time (thyroid 30.8%, meta-analysis of 5,151 nodules; breast 33.5%; prostate 21.3%). Scenario: 19.4% (nationwide Dutch thyroid cohort). |

## Natural history

| Parameter | Value | Derivation |
|---|---|---|
| `sojourn_years` | 2.37 (SE 0.27) | Mean sojourn time of gastric cancer, 2.37 years (95% CI 1.92-2.96). The same value is used for esophageal cancer (**ASSUMPTION**). |
| `p_late_clin` | 0.703 | SEER share of staged cases that are regional or distant: esophagus 71/90, stomach 58/90. Weighted by 2013-17 incidence at ages 50-79 (esophagus 40.5%). |
| `m_early` | 0.079 | SEER 5-year relative survival for localized disease (esophagus 48.6%, stomach 78.1%; weighted 66.1%), converted to a constant annual risk: 1 − S^(1/5). |
| `m_late` | 0.286 | The same conversion for regional and distant disease (esophagus 29.1% / 5.3%, stomach 39.0% / 8.1%; stage-weighted 18.6%). A constant annual risk understates deaths in the first years. |
| `sens_egd` | 0.887 | 11.3% of upper GI cancers were missed at an EGD within the previous 3 years. |
| `inc_ugi` | by age and sex | SEER 21 age-specific incidence, 2013-2017: esophagus (Table 8.7) + stomach (Table 24.7). |
| `prev_ugi` | inc × 2.37 | Preclinical prevalence = incidence × mean sojourn time (steady state). |
| `prev_other` | (all-site − UGI) × 2.37 | All-site incidence (Table 2.7) minus esophagus and stomach. **ASSUMPTION**: same sojourn time, and every other cancer counts as findable by the downstream workup (generous; this matters only through `sens_g_other`). |
| `mort_other` | qx by single year | SSA 2022 period life table (2025 Trustees Report). This is **all-cause** mortality, which slightly overstates competing mortality because cancer deaths are included. |
| Kunzmann groups | RR 2.62 (high) / 0.32 (low) | 29.5% of people score ≥ 8, and they contain 77.4% of esophageal adenocarcinoma (EAC). The RR is applied to the EAC share of esophageal cancer only (63.1%, SEER 2004-15). Stomach risk is unchanged. **ASSUMPTION**: the 29.5% share is the same in every age-sex group (in reality the score rises with age and with male sex). |
| `cohort.csv` | shares | Census Vintage 2024 resident population ages 50-79 by 5-year band and sex, entered at the band midpoint. Both analyses use this cohort. |

## Harms

| Parameter | Value | Derivation |
|---|---|---|
| `sae_egd` | 0.0002 | Perforation 1/2,500-1/11,000. |
| `sae_eus` | 0.0004 | EUS perforation 0.02-0.06%. |
| `sae_eus_samp` | 0.0062 | EUS-FNA bleeding (≈0.18%) + pancreatitis (0.44%). |
| `p_death_sae` | 0.1 | **ASSUMPTION**. A scenario sets it to 0.5, which reproduces the quoted EGD mortality of about 1/10,000. |
| `sae_ent` | 0.0001 | **ASSUMPTION** |
| `msv_pet` | 14.0 (SE 1.3) | PET-CT with low-dose CT. A scenario uses 24.4 mSv (with diagnostic CT). |
| `risk_per_sv` | 0.041 | ICRP 103 nominal risk coefficient for adults. |

## Costs (USD, Medicare rates as the payer proxy)

| Parameter | Value | Derivation |
|---|---|---|
| `c_galleri` | 949 | List price (not covered by Medicare). |
| `c_egd` | 1,001 | CPT 43239, 2024: OPPS $864.59 + physician $136.48. |
| `c_eus` | 2,007 | CPT 43237, 2024: OPPS $1,814.88 + physician $192.07. |
| `c_eus_samp` | 150 | CPT 43238 adds $35.62 over 43237; pathology is assumed. **ASSUMPTION** |
| `c_pet` | 1,459 | CPT 78815, OPPS, July 2025 (facility component only). |
| `c_ent` | 127 | CPT 31575, non-facility, 2026 (excludes the office visit). |
| `c_tx_first_early` | 45,063 | Initial-phase cost of cancer care (USD 2023). |
| `c_fu_early` | 7,141 | Continuing-phase cost. |
| `c_eol` | 121,746 | Last-year-of-life cost. |
| `c_tx_first_late` | 156,982 | First-year cost of stage IV stomach cancer. |
| `c_fu_late` | 21,423 | 3 × the continuing-phase cost. **ASSUMPTION** |
| `c_kunz`, `c_tb`, `c_symp_dx`, `c_sae` | | **ASSUMPTION** |

Costs from different years are not yet inflated to a common year. To do that,
use the medical care component of the CPI.

## Utilities and DALYs

| Parameter | Value | Derivation |
|---|---|---|
| `utility_well` | 0.816 / 0.815 / 0.824 / 0.811 | US EQ-5D-5L norms for ages 45-54 / 55-64 / 65-74 / 75+. |
| `u_mult_early` | 0.91 | Esophageal cancer EQ-5D: stage I 0.66, II/III 0.84, averaged, ÷ the US norm of about 0.82. |
| `u_mult_late` | 0.82 | Stage IV 0.67 ÷ 0.82. |
| `dw_early` | 0.049 | GBD 2019 controlled phase. The diagnosis and primary therapy weight, 0.288, is the upper bound. |
| `dw_late` | 0.451 | GBD 2019 metastatic phase. |
| `daly_ref_le` | 38.10 … 7.05 | WHO GHE 2019 standard years of life lost per death by 5-year age group. The 85+ value is used for all ages 85-100. |
| `d_proc`, `d_sae`, `loss_rad_ca` | | **ASSUMPTION** |

## References

- Neal RD, et al. Performance of a multi-cancer early detection test in the randomized controlled NHS-Galleri trial. *Nat Med* 2026;32:3770-80. doi:10.1038/s41591-026-04652-8
- Nabavizadeh N, et al. Performance and safety of a multi-cancer early detection test: the PATHFINDER 2 study. *Nat Med* 2026;32:3760-9. doi:10.1038/s41591-026-04618-w
- Matrana M, et al. Real-world data and clinical experience from over 100,000 multi-cancer early detection tests. *Nat Commun* 2025;16:9625. doi:10.1038/s41467-025-64094-7
- Talwar A, et al. The role of multicancer early detection tests in head and neck squamous cell carcinoma: a case series. *Head Neck* 2026;48:E126-31. doi:10.1002/hed.70337
- Endoscopic ultrasound-guided Franseen fine-needle biopsy for solid pancreatic lesions: systematic review and meta-analysis. *Dig Dis Sci* 2026. doi:10.1007/s10620-026-09970-y
- de Leijer JF, et al. Focal thyroid incidentalomas on 18F-FDG PET/CT: systematic review and meta-analysis. *Front Endocrinol* 2021.
- Coerts HI, et al. Malignancy risk of thyroid incidentalomas detected by [18F]FDG PET/CT: nationwide study in the Netherlands. *Eur J Nucl Med Mol Imaging* 2026.
- Iacovitti CM, et al. Breast incidental uptake on PET/CT: updated systematic review and meta-analysis. *Pharmaceuticals* 2025.
- Martinello C, et al. Prostatic incidental uptake on PET/CT: updated systematic review and meta-analysis. *Front Oncol* 2026.

- GRAIL. PATHFINDER 2 results fact sheet (2026). Galleri performance in 35,878 adults aged 50+.
- Schrag D, et al. Blood-based tests for multicancer early detection (PATHFINDER). *Lancet* 2023;402:1251-60.
- Klein EA, et al. Clinical validation of a targeted methylation-based multi-cancer early detection test (CCGA3). *Ann Oncol* 2021;32:1167-77.
- Menon S, Trudgill N. How commonly is upper gastrointestinal cancer missed at endoscopy? A meta-analysis. *Endosc Int Open* 2014;2:E46-50.
- Bae JM, et al. Sojourn time of gastric cancer. 2014.
- SEER Cancer Statistics Review 1975-2017, Tables 2.7, 8.7 and 24.7, and SEER stage distribution and 5-year relative survival (seer.cancer.gov).
- Then EO, et al. Esophageal cancer: an updated SEER database analysis. *World J Oncol* 2020;11:55-64.
- Kunzmann AT, et al. Model for identifying individuals at risk for esophageal adenocarcinoma. *Clin Gastroenterol Hepatol* 2018.
- Social Security Administration. Period Life Table, 2022 (2025 Trustees Report). https://www.ssa.gov/oact/STATS/table4c6_2022_TR2025.html
- US Census Bureau. NC-EST2024-SYASEXN, resident population by single year of age and sex, July 1, 2024.
- ASGE Standards of Practice Committee. Adverse events of upper GI endoscopy. *Gastrointest Endosc* 2012;76:707-18.
- Forbes N, et al. Adverse events associated with EUS and EUS-guided procedures (ASGE). *Gastrointest Endosc* 2022;95:16-26.
- Quinn B, et al. Radiation dosimetry of 18F-FDG PET/CT. *J Nucl Med Technol* 2016.
- ICRP Publication 103 (2007).
- Mariotto AB, et al. Medical care costs associated with cancer in the US. *JNCI* 2026.
- Reddy SR, et al. Cost of stomach cancer care in the US. 2022.
- Petric et al. Health-related quality of life in esophageal cancer. *Ann Esophagus* 2025.
- Jiang R, et al. United States EQ-5D-5L population norms. *Qual Life Res* 2021.
- GBD 2019 Disability Weights (IHME).
- WHO. Methods and data sources for global burden of disease estimates 2000-2019 (WHO/DDI/DNA/GHE/2020.3), Table 2.1.
- Medicare physician fee schedule and OPPS rates as compiled in the Cook Medical 2024 reimbursement guide (CPT 43237-43239), CMS OPPS Addendum B July 2025 (78815), and CMS PFS 2026 (31575).
