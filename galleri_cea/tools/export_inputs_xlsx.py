"""Export every model input to one Excel workbook (the mirror of import_inputs_xlsx.R).

    python3 tools/export_inputs_xlsx.py [output.xlsx]     (run from galleri_cea/)

Reads inputs/*.csv and data-raw/raw/*.csv only. Default output: Galleri_model_inputs.xlsx.
Formulas (Kunzmann RRs, upper GI totals, late shares) get values when the file is
opened in Excel or LibreOffice; the importer never relies on them.
"""
import csv, os, sys
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table, TableStyleInfo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "Galleri_model_inputs.xlsx")
rd = lambda f: list(csv.DictReader(open(f"{ROOT}/{f}", newline="", encoding="utf-8")))
REFS = [(int(r["ref"]), r["citation"], r["url"], r["pmid"]) for r in rd("inputs/references.csv")]
NONE_REF = "none (assumption)"

F = "Arial"
H_FILL = PatternFill("solid", fgColor="1F3A5F")
H_FONT = Font(name=F, bold=True, color="FFFFFF", size=10)
BASE = Font(name=F, size=10)
BLUE = Font(name=F, size=10, color="0000FF")
BOLD = Font(name=F, size=10, bold=True)
SEC_FILL = PatternFill("solid", fgColor="E3E8EF")
STATUS_FILL = {"SOURCED": "D9EAD3", "ASSUMPTION": "FCE5CD", "DATA": "CFE2F3", "DUMMY": "F4CCCC"}
YELLOW = PatternFill("solid", fgColor="FFFF00")
thin = Side(style="thin", color="D5DCE5")

wb = Workbook()

def sheet(title, headers, rows, widths, fmts=None, value_cols=(), wrap_cols=(), freeze="A2", table=True):
    ws = wb.create_sheet(title)
    ws.append(headers)
    for c in ws[1]:
        c.font = H_FONT; c.fill = H_FILL; c.alignment = Alignment(vertical="center", wrap_text=True)
    for r in rows:
        ws.append(r)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = BLUE if c.column in value_cols else BASE
            c.alignment = Alignment(vertical="top", wrap_text=c.column in wrap_cols)
            if fmts and c.column in fmts: c.number_format = fmts[c.column]
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = freeze
    ws.row_dimensions[1].height = 30
    if table and len(rows):
        ref = f"A1:{get_column_letter(len(headers))}{len(rows) + 1}"
        t = Table(displayName=title.replace(" ", "_").replace("-", "_"), ref=ref)
        t.tableStyleInfo = TableStyleInfo(name="TableStyleLight1", showRowStripes=True)
        ws.add_table(t)
    return ws

def num(x):
    try:
        v = float(x); return int(v) if v.is_integer() and "." not in x and "e" not in x.lower() else v
    except (TypeError, ValueError):
        return x if x != "" else None

# ---------------------------------------------------------------- README
IMPORT_STEPS = [
    ("How to update the model from this workbook", "h"),
    ("1. Edit blue cells (and add or delete rows on Parameters, Scenarios, References and the raw-data sheets). Black cells are formulas or labels; the Baseline (calibration) sheet is rebuilt from the raw-data sheets, so edits there are ignored.", ""),
    ("2. Save the workbook as galleri_cea/Galleri_model_inputs.xlsx (keep the sheet names and column headers).", ""),
    ("3. From galleri_cea/, preview: Rscript import_inputs_xlsx.R --check   (validates and lists every change, writes nothing).", ""),
    ("4. Apply: Rscript import_inputs_xlsx.R   (validates on a scratch copy, backs up the old CSVs, writes the new ones, rebuilds the SEER/NCHS baseline).", ""),
    ("5. Run the model: Rscript run.R", ""),
]
NO_IMPORT = [("Editing this workbook does not yet change the model: edit the CSVs in inputs/ (or data-raw/raw/, then run Rscript data-raw/make_us_inputs.R).", "")]
ws = wb.active; ws.title = "README"
lines = [
    ("Galleri upper GI screening model: all inputs", "title"),
    ("Source: github.com/fali-ebm/Galleri, folder galleri_cea. Exported from the files the model reads (tools/export_inputs_xlsx.py).", ""),
    ("", ""),
] + (IMPORT_STEPS if os.path.exists(os.path.join(ROOT, 'import_inputs_xlsx.R')) else NO_IMPORT) + [
    ("", ""),
    ("Model in one paragraph", "h"),
    ("Annual-cycle cohort model, lifetime horizon (to age 100), US payer perspective, 3% discounting of costs and QALYs. Arms: usual care; Kunzmann score >= 8 -> EGD every year; Galleri every year: +ve -> EGD (EGD +ve = true positive; EGD -ve = false positive -> rest of the workup; any EGD can perforate), -ve -> true negative (stays well) or false negative (symptomatic cancer later or found at a later round), then usual care. Screening runs from the start age to 85. Natural history is calibrated to SEER for every subgroup: onset of preclinical upper GI cancer is solved so that, without screening, symptomatic diagnoses reproduce SEER incidence by age, sex and race/ethnicity, and their stage reproduces SEER stage at diagnosis. Base case: 50-year-old non-Hispanic White men.", ""),
    ("", ""),
    ("Sheets", "h"),
    ("Parameters - every scalar input with value, PSA distribution, range, status and source (one row per parameter).", ""),
    ("References - numbered sources; Refs columns elsewhere point to these numbers.", ""),
    ("Settings - analysis settings: willingness to pay, discounting, horizon, screening schedule, Kunzmann split, subgroup grid, PSA size.", ""),
    ("Cohort - base-case cohort rows (start age, sex, race/ethnicity). The Kunzmann shares are linked to Settings.", ""),
    ("Scenarios - deterministic scenario analyses (parameter overrides).", ""),
    ("Workup counts - template for local Galleri +ve workup counts; filled counts override the workup probabilities.", ""),
    ("Baseline (calibration) - life tables and SEER calibration targets by race/ethnicity x sex x single year of age 50-99 (600 rows).", ""),
    ("SEER incidence, SEER race ratios, SEER stage, EAC share, NCHS life tables, EQ-5D norms, WHO standard YLL - the raw published tables the baseline is built from.", ""),
    ("", ""),
    ("Status legend", "h"),
    ("SOURCED", "SOURCED"), ("ASSUMPTION", "ASSUMPTION"), ("DATA", "DATA"), ("DUMMY", "DUMMY"),
    ("", ""),
    ("Conventions", "h"),
    ("Blue text = editable input values; black = formulas, labels or derived values. Probabilities are stored as fractions and shown as %. Costs in USD (Medicare rates as payer proxy; years 2023-2026, not yet inflated to a common year). Rates in the baseline sheet are per person-year; per-100,000 columns are provided for reading.", ""),
    ("Yellow cells on the Workup counts sheet are the ones to fill in (n column).", ""),
    ("Full derivations and references: galleri_cea/inputs/SOURCES.md.", ""),
]
desc = {"SOURCED": "Taken from, or derived arithmetically from, a published source.",
        "ASSUMPTION": "No direct source; value and reasoning documented so it can be challenged.",
        "DATA": "Set from local data (counts.csv).", "DUMMY": "Placeholder with no source (none remain)."}
for i, (txt, kind) in enumerate(lines, 1):
    c = ws.cell(i, 1, txt); c.font = BASE; c.alignment = Alignment(wrap_text=True, vertical="top")
    if kind == "title": c.font = Font(name=F, size=14, bold=True)
    elif kind == "h": c.font = Font(name=F, size=11, bold=True)
    elif kind in STATUS_FILL:
        c.fill = PatternFill("solid", fgColor=STATUS_FILL[kind]); c.font = BOLD
        d = ws.cell(i, 2, desc[kind]); d.font = BASE
ws.column_dimensions["A"].width = 110; ws.column_dimensions["B"].width = 70

# ---------------------------------------------------------------- Parameters
CATS = [
    ("Galleri test performance", ["sens_g_ugi", "sens_g_other", "fpr_g_es"]),
    ("Endoscopy (Kunzmann arm)", ["sens_egd"]),
    ("Galleri +ve workup tree", ["p_egd_pos", "p_ent_route", "p_ent_pos", "p_pet_pb", "p_pet_incid", "p_eus_pos",
                                 "p_eus_neg", "p_miss", "p_incid_cancer"]),
    ("Natural history and stage", ["sojourn_years", "p_late_screen", "or_late_clin", "m_early", "m_late"]),
    ("Procedure harms", ["p_perf_egd", "p_death_perf", "sae_ent", "sae_eus", "sae_eus_samp", "p_death_sae", "msv_pet", "risk_per_sv",
                         "loss_rad_ca"]),
    ("Time to workup steps", ["t_egd", "t_ent", "t_pet", "t_eus", "t_tb"]),
    ("Costs (USD)", ["c_galleri", "c_kunz", "c_egd", "c_ent", "c_pet", "c_eus", "c_eus_samp", "c_tb", "c_symp_dx",
                     "c_tx_first_early", "c_tx_first_late", "c_fu_early", "c_fu_late", "c_eol", "c_sae", "c_perf"]),
    ("Utilities, disutilities and disability weights", ["u_mult_early", "u_mult_late", "d_proc", "d_sae",
                                                        "dw_early", "dw_late"]),
]
P = {r["parameter"]: r for r in rd("inputs/parameters.csv")}
assert sorted(sum((v for _, v in CATS), [])) == sorted(P), set(P) ^ set(sum((v for _, v in CATS), []))
hdr = ["Category", "Parameter", "Label", "Value", "Unit", "Status", "PSA distribution", "SE", "n (beta/Dirichlet)",
       "Dirichlet group", "Low (OWSA)", "High (OWSA)", "Description", "Source / derivation", "Refs (see References)", "Used in"]
ws = wb.create_sheet("Parameters"); ws.append(hdr)
for c in ws[1]:
    c.font = H_FONT; c.fill = H_FILL; c.alignment = Alignment(vertical="center", wrap_text=True)
for cat, keys in CATS:
    ws.append([cat]); r = ws.max_row
    for c in ws[r]: c.fill = SEC_FILL
    ws.cell(r, 1).font = BOLD
    for k in keys:
        p = P[k]
        ws.append([None, k, p["label"], num(p["value"]), p["unit"], p["status"], p["dist"], num(p["se"]), num(p["n"]),
                   p["group"] or None, num(p["low"]), num(p["high"]), p["description"], p["source"], p.get("refs") or NONE_REF, {"both": "Both analyses", "comparison": "Comparison (Kunzmann arm)"}.get(p["analysis"], p["analysis"])])
        r = ws.max_row
        for c in ws[r]:
            c.font = BASE; c.alignment = Alignment(vertical="top", wrap_text=c.column in (3, 13, 14))
        for col in (4, 11, 12):
            ws.cell(r, col).font = BLUE
            u = p["unit"]
            fmt = ("0.0%" if u == "proportion" and float(p["value"] or 0) >= 0.01 else
                   "0.000%" if u in ("proportion", "per procedure", "per year") else
                   "$#,##0" if u.startswith("USD") else "0.00###")
            if u == "per year" and k.startswith("m_"): fmt = "0.0%"
            ws.cell(r, col).number_format = fmt
        ws.cell(r, 6).fill = PatternFill("solid", fgColor=STATUS_FILL.get(p["status"], "FFFFFF"))
widths = [24, 18, 30, 12, 13, 13, 13, 9, 10, 10, 11, 11, 48, 60, 14, 12]
for i, w in enumerate(widths, 1): ws.column_dimensions[get_column_letter(i)].width = w
ws.freeze_panes = "C2"; ws.row_dimensions[1].height = 30
ws.auto_filter.ref = f"A1:P{ws.max_row}"
n0 = ws.max_row + 2
ws.cell(n0, 2, "Notes").font = BOLD
for j, t in enumerate([
    "Value = base case. SE / n define the probabilistic (PSA) distribution: beta and Dirichlet use n as the effective sample size; gamma, normal and lognormal use SE.",
    "Low/High are one-way sensitivity bounds; where blank the model uses +/-25% (owsa_pct) of the base value.",
    "EGD perforation (p_perf_egd, p_death_perf, c_perf) applies to every EGD in both screening arms; a fatal perforation ends the Galleri workup. p_death_sae, c_sae and d_sae now cover the other serious adverse events (ENT, EUS, sampling); d_sae is also the QALY loss of a non-fatal perforation.",
    "Stage at symptomatic diagnosis is no longer a single parameter: it is calibrated to SEER by site, race and sex (Baseline sheet, late_eso / late_stom). or_late_clin = 1 reproduces SEER; it carries the PSA uncertainty.",
    "Workup tree probabilities are literature values (PubMed/Consensus, Nov 2026); counts entered on the Workup counts sheet (inputs/counts.csv) override them."], 1):
    c = ws.cell(n0 + j, 2, t); c.font = BASE
    ws.merge_cells(start_row=n0 + j, start_column=2, end_row=n0 + j, end_column=14)
    c.alignment = Alignment(wrap_text=True, vertical="top"); ws.row_dimensions[n0 + j].height = 28

# ---------------------------------------------------------------- Settings
ST = rd("inputs/settings.csv")
FMT = {"wtp": "$#,##0", "disc_cost": "0.0%", "disc_qaly": "0.0%", "disc_daly": "0.0%", "kunz_share": "0.0%",
       "kunz_capture": "0.0%", "owsa_pct": "0%", "n_sim": "#,##0"}
ws = wb.create_sheet("Settings")
ws.append(["Area", "Setting", "Key", "Value", "Note", "Refs"])
for c in ws[1]:
    c.font = H_FONT; c.fill = H_FILL; c.alignment = Alignment(vertical="center", wrap_text=True)
row_of = {}
for r in ST:
    ws.append([r["area"], r["label"], r["setting"], num(r["value"]), r["note"], r["refs"] or None])
    i = ws.max_row; row_of[r["setting"]] = i
    for c in ws[i]: c.font = BASE; c.alignment = Alignment(vertical="top", wrap_text=c.column in (2, 4, 5))
    ws.cell(i, 4).font = BLUE; ws.cell(i, 4).alignment = Alignment(horizontal="left", vertical="top", wrap_text=True)
    if r["setting"] in FMT: ws.cell(i, 4).number_format = FMT[r["setting"]]
ks, kc = row_of["kunz_share"], row_of["kunz_capture"]
lst = lambda k: f'(LEN(D{row_of[k]})-LEN(SUBSTITUTE(D{row_of[k]},",",""))+1)'
derived = [
    ("Kunzmann score", "Implied relative risk, flagged", f"=D{kc}/D{ks}", "applied to EAC share of esophageal cancer", "8, 9", "0.00"),
    ("Kunzmann score", "Implied relative risk, not flagged", f"=(1-D{kc})/(1-D{ks})", "applied to EAC share of esophageal cancer", "8, 9", "0.00"),
    ("Subgroup grid", "Number of subgroups", f"={lst('subgroup_ages')}*{lst('subgroup_sexes')}*{lst('subgroup_races')}", "deterministic runs, each re-calibrated", None, "0"),
    ("Model structure", "Cycle length", 1, "year (fixed in model code)", None, "0"),
    ("Model structure", "Half-cycle correction", "Yes", "mid-cycle alive fraction (fixed in model code)", None, None),
    ("Economic", "Perspective", "US payer", "Medicare rates as proxy (fixed)", None, None),
    ("Screening", "Galleri arm", "Galleri every round; +ve -> full workup tree", "fixed in model code", "1, 2", None),
    ("Screening", "Kunzmann arm", "Score once at start; EGD every round if score >= 8", "fixed in model code", None, None),
    ("Base case", "Start age, sex, race/ethnicity", "set on the Cohort sheet", "", None, None),
]
ws.append([]); ws.append(["Derived and fixed (not imported)"]); ws.cell(ws.max_row, 1).font = BOLD
for c in ws[ws.max_row]: c.fill = SEC_FILL
for a, l, v, n, rf, fm in derived:
    ws.append([a, l, None, v, n, rf]); i = ws.max_row
    for c in ws[i]: c.font = BASE; c.alignment = Alignment(vertical="top", wrap_text=c.column in (2, 4, 5))
    ws.cell(i, 4).alignment = Alignment(horizontal="left", vertical="top", wrap_text=True)
    if fm: ws.cell(i, 4).number_format = fm
for i, w in enumerate([18, 44, 16, 40, 48, 10], 1): ws.column_dimensions[get_column_letter(i)].width = w
ws.freeze_panes = "A2"; ws.row_dimensions[1].height = 30
ws.cell(ws.max_row + 2, 1, "Only rows with a Key are imported (blue Value cells). Lists are comma-separated.").font = Font(name=F, size=9, italic=True)

# ---------------------------------------------------------------- Cohort
C = rd("inputs/cohort.csv")
ws = sheet("Cohort", list(C[0].keys()), [[num(v) for v in r.values()] for r in C],
           [12, 10, 8, 10, 16, 10, 16, 13, 80], fmts={6: "0.0%"}, value_cols=(2, 3, 4), wrap_cols=(9,))
for i, r in enumerate(C, 2):
    link = {"kunzmann_high": f"=Settings!D{ks}", "kunzmann_low": f"=1-Settings!D{ks}"}.get(r["risk_group"])
    if link: ws.cell(i, 6).value = link
    else: ws.cell(i, 6).font = BLUE
ws.cell(len(C) + 3, 1, "Shares of the Kunzmann groups are linked to Settings (kunz_share) and set from it on import. Ref [8].").font = Font(name=F, size=9, italic=True)

# ---------------------------------------------------------------- Scenarios
SC = rd("inputs/scenarios.csv")
lab = {k: P[k]["label"] for k in P}
rows = [[r["scenario"], r["parameter"], lab.get(r["parameter"], ""), num(P[r["parameter"]]["value"]), num(r["value"]),
         r.get("refs") or NONE_REF] for r in SC]
ws = sheet("Scenarios", ["Scenario", "Parameter", "Label", "Base value", "Scenario value", "Refs"], rows, [58, 16, 32, 13, 14, 16],
           value_cols=(1, 2, 5, 6))

# ---------------------------------------------------------------- Workup counts
CT = rd("inputs/counts.csv")
rows = [[r["branch"], r["outcome"], r["parameter"] or None, num(r["n"]), r["notes"]] for r in CT]
ws = sheet("Workup counts", ["Branch", "Outcome", "Sets parameter", "n (patients)", "What to count"], rows,
           [12, 20, 16, 13, 70], value_cols=(4,), wrap_cols=(5,), table=False)
for i in range(2, len(rows) + 2): ws.cell(i, 4).fill = YELLOW
k = len(rows) + 3
for j, t in enumerate(["Fill the yellow n column with numbers of patients. Fill every row of a branch, or leave the whole branch blank (a partly filled branch stops the run).",
                       "Example: egd positive = 42, egd negative = 58 gives P(EGD +ve) = (42 + 0.5) / (100 + 1) = 0.42 with a beta(42.5, 58.5) PSA distribution.",
                       "Then run Rscript import_inputs_xlsx.R (writes inputs/counts.csv)."]):
    c = ws.cell(k + j, 1, t); c.font = BASE

# ---------------------------------------------------------------- Baseline
B = rd("inputs/baseline_risks.csv")
RL = {"all": "All races/ethnicities", "nhw": "Non-Hispanic White", "nhb": "Non-Hispanic Black",
      "nhapi": "Non-Hispanic Asian/Pacific Islander", "nhaian": "Non-Hispanic American Indian/Alaska Native",
      "hisp": "Hispanic"}
rows = []
for r in B:
    rows.append([r["race_eth"], RL[r["race_eth"]], r["sex"], int(r["age_from"]), float(r["mort_other"]),
                 float(r["inc_eso"]) * 1e5, float(r["inc_stom"]) * 1e5, None, float(r["inc_other"]) * 1e5,
                 float(r["eac_share"]), float(r["late_eso"]), float(r["late_stom"]), float(r["utility_well"]),
                 float(r["daly_ref_le"])])
hdr = ["race_eth", "Race/ethnicity", "Sex", "Age", "Other-cause death prob. qx (NCHS 2024) [30]",
       "SEER esophagus incidence per 100,000 [26]", "SEER stomach incidence per 100,000 [27]",
       "Upper GI target per 100,000 (formula)", "Other cancers incidence per 100,000 [26-28]",
       "Adenocarcinoma share of esophageal [26]", "Late-stage share at dx, esophagus [26]",
       "Late-stage share at dx, stomach [27]", "Utility, general population [17] (EQ-5D norms sheet)",
       "Standard YLL per death [19] (WHO standard YLL sheet)"]
ws = sheet("Baseline (calibration)", hdr, rows, [9, 30, 8, 6, 14, 14, 14, 14, 14, 13, 13, 13, 12, 12],
           fmts={5: "0.000000", 6: "0.00", 7: "0.00", 8: "0.00", 9: "0.0", 10: "0.0%", 11: "0.0%", 12: "0.0%",
                 13: "0.000", 14: "0.00"}, value_cols=(5, 6, 7, 9, 10, 11, 12, 13, 14), freeze="E2")
for i in range(2, len(rows) + 2):
    ws.cell(i, 8).value = f"=F{i}+G{i}"; ws.cell(i, 8).font = BASE
ws.row_dimensions[1].height = 66
ws.cell(len(rows) + 4, 1, "DERIVED SHEET - not imported. Rebuilt from the raw-data sheets by data-raw/make_us_inputs.R on import; edit those sheets instead.").font = Font(name=F, size=10, bold=True, color="C00000")
ws.cell(len(rows) + 3, 1, "Numbers in [ ] refer to the References sheet. Race/ethnicity groups without SEER age-specific rates are scaled by age-adjusted ratios [26-28]; stage for groups other than White/Black uses all races.").font = Font(name=F, size=9, italic=True)

# ---------------------------------------------------------------- Raw source tables
def raw(f): return rd(f"data-raw/raw/{f}")
R1 = raw("seer_csr_2013_2017.csv")
sheet("SEER incidence", ["Site", "Race", "Sex", "Age band (start)", "Rate per 100,000"],
      [[r["site"], r["race"], r["sex"], f'{r["age_band"]}{"+" if r["age_band"] == "85" else "-" + str(int(r["age_band"]) + 4)}',
        float(r["rate"])] for r in R1], [12, 10, 10, 16, 16], fmts={5: "#,##0.0"}, value_cols=(5,))
R2 = raw("seer_csr_age_adjusted_by_race.csv")
sheet("SEER race ratios", ["Site", "Race/ethnicity", "Male, age-adjusted per 100,000", "Female, age-adjusted per 100,000"],
      [[r["site"], r["race"], float(r["male"]), float(r["female"])] for r in R2], [12, 14, 18, 18],
      fmts={3: "#,##0.0", 4: "#,##0.0"}, value_cols=(3, 4))
R3 = raw("seer_stage_by_race_sex.csv")
ws = sheet("SEER stage", ["Site", "Race", "Sex", "Localized %", "Regional %", "Distant %", "Unstaged %",
                          "Late share among staged (formula)"],
           [[r["site"], r["race"], r["sex"], int(r["localized"]), int(r["regional"]), int(r["distant"]),
             int(r["unstaged"]), None] for r in R3], [12, 10, 10, 12, 12, 12, 12, 16], value_cols=(4, 5, 6, 7))
for i in range(2, len(R3) + 2):
    ws.cell(i, 8).value = f"=(E{i}+F{i})/(D{i}+E{i}+F{i})"; ws.cell(i, 8).number_format = "0.0%"; ws.cell(i, 8).font = BASE
R4 = raw("seer_eac_share_by_race.csv")
sheet("EAC share", ["race_eth", "Race/ethnicity", "Adenocarcinoma share of esophageal cancer"],
      [[r["race"], RL[r["race"]], float(r["eac_share"])] for r in R4], [10, 40, 22], fmts={3: "0.0%"}, value_cols=(3,))
R5 = raw("nchs_lifetables_2024_qx.csv")
groups = list(dict.fromkeys((r["race_eth"], r["sex"]) for r in R5))
q = {(r["race_eth"], r["sex"], int(r["age"])): float(r["qx"]) for r in R5}
hdr = ["Age"] + [f"{RL[g]} {s}" for g, s in groups]
rows = [[a] + [q[(g, s, a)] for g, s in groups] for a in range(50, 100)]
ws = sheet("NCHS life tables", hdr, rows, [6] + [16] * len(groups),
           fmts={i: "0.000000" for i in range(2, len(groups) + 2)}, value_cols=tuple(range(2, len(groups) + 2)),
           freeze="B2", table=False)
ws.row_dimensions[1].height = 54
R6 = raw("eq5d_us_norms.csv")
sheet("EQ-5D norms", ["Age from", "Age to", "EQ-5D-5L utility"], [[int(r["age_from"]), int(r["age_to"]), float(r["utility"])] for r in R6],
      [10, 10, 16], fmts={3: "0.000"}, value_cols=(1, 2, 3))
R7 = raw("who_ghe_standard_yll.csv")
sheet("WHO standard YLL", ["Age (band start; 85 = 85+)", "Standard YLL per death"], [[int(r["age"]), float(r["yll"])] for r in R7],
      [22, 20], fmts={2: "0.00"}, value_cols=(1, 2))
src = {"SEER incidence": "[26, 27, 28] SEER Cancer Statistics Review 1975-2017, Tables 2.7 (all sites), 8.7 (esophagus), 24.7 (stomach): age-specific SEER 21 incidence 2013-2017.",
       "SEER race ratios": "[26, 27, 28] SEER CSR 1975-2017 Tables 2.15, 8.15, 24.15: age-adjusted SEER 21 incidence 2013-2017 by race/ethnicity and sex.",
       "SEER stage": "[26, 27] SEER CSR 1975-2017 Tables 8.8, 24.8: stage distribution 2010-2016 (SEER Summary Stage 2000).",
       "EAC share": "[26] SEER CSR 1975-2017 Table 8.22: histology by race, 2013-2017 (White used for NH White, Black for NH Black).",
       "EQ-5D norms": "[17] Jiang R, Janssen MFB, Pickard AS. Qual Life Res 2021;30:803-16 (US EQ-5D-5L norms by age band). Bands must cover ages 50-99 once each.",
       "WHO standard YLL": "[19] WHO GHE 2019 DALY methods, Table 2.1 standard life table: remaining life expectancy at the start of each 5-year band.",
       "NCHS life tables": "[30] Arias E, Xu JQ, Tejada-Vera B. United States Life Tables, 2024. NVSR vol 75 no 5 (Tables 2-3, 5-6, 8-9, 11-12, 14-15, 17-18)."}
for name, s in src.items():
    w = wb[name]; r0 = w.max_row + 2 if name != "NCHS life tables" else w.max_row + 2
    c = w.cell(r0, 1, "Source: " + s); c.font = Font(name=F, size=9, italic=True)

# ---------------------------------------------------------------- References
used = {}
for k, p in P.items():
    for n in [x.strip() for x in (p.get("refs") or "").split(",") if x.strip()]:
        used.setdefault(int(n), []).append(k)
extra = {8: ["Kunzmann split (Settings, Cohort)"], 9: ["Kunzmann RR on EAC share (Settings)"],
         17: ["utility_well (EQ-5D norms)"], 19: ["daly_ref_le (WHO standard YLL)"], 21: ["scenario: p_incid_cancer"],
         26: ["inc_eso, eac_share, late_eso (Baseline); SEER incidence/stage sheets"],
         27: ["inc_stom, late_stom (Baseline); SEER incidence/stage sheets"],
         28: ["inc_other (Baseline); SEER incidence/race-ratio sheets"], 30: ["mort_other (Baseline); NCHS life tables sheet"]}
for n, v in extra.items(): used.setdefault(n, []).extend(v)
ws = wb.create_sheet("References", 1)
ws.append(["#", "Reference", "DOI / URL", "PMID", "Used for (derived, not imported)"])
for c in ws[1]:
    c.font = H_FONT; c.fill = H_FILL; c.alignment = Alignment(vertical="center", wrap_text=True)
for n, cit, url, pmid in REFS:
    ws.append([n, cit, url or None, int(pmid) if pmid else None, ", ".join(used.get(n, []))])
    r = ws.max_row
    for c in ws[r]:
        c.font = BLUE if c.column in (1, 2, 3, 4) else BASE; c.alignment = Alignment(vertical="top", wrap_text=c.column in (2, 3, 5))
    if url:
        ws.cell(r, 3).hyperlink = url; ws.cell(r, 3).font = Font(name=F, size=10, color="0563C1", underline="single")
    if pmid:
        ws.cell(r, 4).hyperlink = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        ws.cell(r, 4).font = Font(name=F, size=10, color="0563C1", underline="single")
for i, w in enumerate([5, 80, 45, 11, 45], 1): ws.column_dimensions[get_column_letter(i)].width = w
ws.freeze_panes = "B2"
r0 = ws.max_row + 2
ws.cell(r0, 2, "Inputs marked 'none (assumption)' have no published source; their reasoning is in the Source / derivation column of the Parameters sheet. Citations and DOIs verified in PubMed (Nov 2026) except where noted.").font = Font(name=F, size=9, italic=True)
ws.cell(r0, 2).alignment = Alignment(wrap_text=True)

# README pointer
rd_ws = wb["README"]
rd_ws.cell(rd_ws.max_row + 2, 1, "References - numbered citations with DOI/PubMed links; the Refs columns on other sheets point to these numbers.").font = BASE

wb.save(OUT)
print("saved", OUT)
