"""Build viz/galleri_burden.html from viz/galleri_viz_data.json (python3 viz/build_burden.py)."""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "galleri_viz_data.json")))
G = [r for r in D["by_age"]["Galleri"] if r["age"] <= D["meta"]["stop_age"]]
R = {r["key"]: r for r in D["galleri"] if r["arm"] == "Galleri"}
P = {r["parameter"]: r["value"] for r in D["params"]}

def largest_remainder(values, total):
    fl = [math.floor(v) for v in values]; rem = total - sum(fl)
    order = sorted(range(len(values)), key=lambda i: values[i] - fl[i], reverse=True)
    for i in order[:rem]: fl[i] += 1
    return fl

# people with a positive result, by fate and age (whole people for the dot field)
groups = {
  "found":   [r["found_ugi"] + r["found_other"] - r["overdx_ugi"] - r["overdx_other"] for r in G],
  "overdx":  [r["overdx_ugi"] + r["overdx_other"] for r in G],
  "missed":  [r["tp"] + r["fp_other_ca"] - r["found_ugi"] - r["found_other"] for r in G],
  "falsealarm": [r["free_inv"] for r in G],
}
tot = {k: round(sum(v)) for k, v in groups.items()}
alloc = {k: largest_remainder(v, tot[k]) for k, v in groups.items()}
ages = [{"age": r["age"], "tests": r["tests"], "alive": r["alive"], "tp": r["tp"], "fn": r["fn"], "fn_symp": r["fn_symp"],
         "found_ugi": r["found_ugi"], "found_other": r["found_other"], "overdx": r["overdx_ugi"] + r["overdx_other"],
         "pos": r["tp"] + r["fp"], **{k: alloc[k][i] for k in alloc}} for i, r in enumerate(G)]
n_fp = R["free_investigated"]["value"]
gfree = {k: R[k]["value"] for k in ["ent_free", "pet_ct_free", "eus_free", "eus_sampling_free", "incidental_followups_free",
                                    "perforations_free", "proc_deaths_free", "dose_sv_free", "radiation_cancers_free"]}
ledger = {
  "n_fp": n_fp, "ent": gfree["ent_free"] / n_fp, "pet": gfree["pet_ct_free"] / n_fp, "eus": gfree["eus_free"] / n_fp,
  "samp": gfree["eus_sampling_free"] / n_fp, "incid": gfree["incidental_followups_free"] / n_fp,
  "days": R["mean_days_fp"]["value"], "cost": R["cost_fp"]["value"] * 1e5 / n_fp, "msv": gfree["dose_sv_free"] * 1000 / n_fp,
  "qaly_days": R["qaly_fp"]["value"] / n_fp * 365.25, "perf_1_in": 1 / P["p_perf_egd"],
  "cost_total": R["cost_fp"]["value"] * 1e5, "qaly_total": R["qaly_fp"]["value"],
}
def ci(k, scale=1, d=0):
    r = R[k]; f = lambda v: f"{v*scale:,.{d}f}"
    return {"v": f(r["value"]), "lo": f(r["lower_95"]), "hi": f(r["upper_95"])}
table = [
  ("Galleri tests over 35 years", ci("galleri_tests")),
  ("People with a positive result, worked up", ci("investigated")),
  ("Cancers found (upper GI and other sites)", ci("screen_detected")),
  ("  of which overdiagnosed", ci("overdiagnosed")),
  ("Cancer-free people worked up (false alarms)", ci("free_investigated")),
  ("PET CT scans in cancer-free people", ci("pet_ct_free")),
  ("EGD perforations, all people worked up", ci("perforations", d=2)),
  ("Radiation-induced cancers, all people worked up", ci("radiation_cancers", d=2)),
  ("False-negative test results", ci("screen_fn")),
  ("  presenting with symptoms within a year", ci("fn_symptomatic_1y")),
  ("QALYs lost to the workup and its harms, everyone worked up", ci("qaly_harms", d=1)),
  ("QALYs lost by cancer-free people worked up", ci("qaly_fp", d=1)),
  ("DALYs incurred by cancer-free people worked up", ci("daly_fp", d=1)),
]
fpc = lambda k: R[k]["value"]
fpb = {"counts": [
    ["Negative EGDs", fpc("egd_free")], ["Negative PET CT scans", fpc("pet_ct_free")], ["Negative ENT exams", fpc("ent_free")],
    ["Negative endoscopic ultrasounds", fpc("eus_free")], ["  of which with a biopsy", fpc("eus_sampling_free")],
    ["Negative diagnostic laparoscopies", fpc("laparoscopy_free")], ["Watchful-waiting imaging", fpc("surveillance_free")],
    ["Follow-ups of a harmless PET finding", fpc("incidental_followups_free")], ["Days waiting for the all-clear", fpc("days_free")],
    ["Still worried 6 months later", fpc("fp_worry_cases")], ["EGD perforations", fpc("perforations_free")],
    ["Other serious complications", fpc("fp_sae_other")], ["Deaths from complications", fpc("proc_deaths_free")],
    ["Radiation-induced cancers, later in life", fpc("radiation_cancers_free")]],
  "comp": [
    ["Workup: procedures and uncertainty", fpc("fp_qaly_pathway"), fpc("fp_daly_anxiety")],
    ["Persisting worry", fpc("fp_qaly_worry"), fpc("fp_daly_worry")],
    ["Radiation-induced cancers", fpc("fp_qaly_radiation"), fpc("fp_daly_radiation")],
    ["Deaths from complications", fpc("fp_qaly_deaths"), fpc("fp_daly_deaths")],
    ["Non-fatal complications", fpc("fp_qaly_complications"), fpc("fp_daly_complications")]],
  "qaly": [fpc("qaly_fp"), R["qaly_fp"]["lower_95"], R["qaly_fp"]["upper_95"]],
  "daly": [fpc("daly_fp"), R["daly_fp"]["lower_95"], R["daly_fp"]["upper_95"]],
  "icer": R["icer"]["value"], "qaly_net": R["qaly_net"]["value"]}
SC = {r["scenario"]: r["icer"] for r in D["scenarios"] if r["arm"] == "Galleri"}
fpb["icer_eq5d"] = SC.get("False positives: per-procedure loss only (EQ-5D finds none)")
acc = [{k: r[k] for k in ["sex", "race_eth", "race_label", "age", "tests", "tp", "fp", "fp_other_cancer", "fp_cancer_free", "tn", "fn",
                          "prevalence_per_100k_tested", "ppv_upper_gi", "ppv_any_cancer", "npv", "cancer_free_per_cancer_found"]}
       for r in D["accuracy"]]
data = {"fpb": fpb, "acc": acc, "ages": ages, "tot": tot, "links": D["workup_links"], "ledger": ledger, "table": table,
        "tests_total": R["galleri_tests"]["value"], "investigated": R["investigated"]["value"],
        "found": R["screen_detected"]["value"], "overdx": R["overdiagnosed"]["value"], "fn": R["screen_fn"]["value"],
        "fn_symp": R["fn_symptomatic_1y"]["value"], "price": P["c_galleri"], "cost_year": D["meta"]["cost_year"]}

html = open(os.path.join(HERE, "burden_template.html")).read().replace("__DATA__", json.dumps(data, separators=(",", ":")))
open(os.path.join(HERE, "galleri_burden.html"), "w").write(html)
print("wrote galleri_burden.html", len(html), tot)
