"""Build viz/galleri_vs_usual_care.html (python3 viz/build_compare.py). Shares the burden piece's style."""
import json, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "galleri_viz_data.json")))
R = {(r["key"], r["arm"]): r for r in D["comparison"]}
v = lambda k, a="Galleri": R[(k, a)]["value"]
UC, G = D["by_age"]["Usual care"], D["by_age"]["Galleri"]
sc = {r["scenario"]: r for r in D["scenarios"] if r["arm"] == "Galleri"}
base, half = sc["Base case"]["icer"], sc["Galleri test half price"]["icer"]
slope = (base - half) / (D["meta"]["c_galleri"] / 2); icpt = base - slope * D["meta"]["c_galleri"]
d_cost = v("cost_total") - v("cost_total", "Usual care")          # per person, discounted
ugi_uc, ugi_g = sum(r["ca_death"] for r in UC), sum(r["ca_death"] for r in G)
data = {
  "per": D["meta"]["per"], "price": D["meta"]["c_galleri"], "cost_year": D["meta"]["cost_year"], "wtp": D["meta"]["wtp"],
  "icer": base, "icer_lo": R[("icer", "Galleri")]["lower_95"], "icer_hi": R[("icer", "Galleri")]["upper_95"],
  "d_cost_total": d_cost * 1e5, "qaly_net": v("qaly_net"), "ly_gained": v("ly_gained"), "late_averted": v("late_averted"),
  "deaths_averted": v("ca_deaths_averted"), "ugi_deaths_uc": ugi_uc, "ugi_deaths_g": ugi_g,
  "other_deaths_averted": v("ca_deaths_averted") - (ugi_uc - ugi_g), "tests": v("galleri_tests"),
  "overdx": v("overdiagnosed"), "false_alarms": v("free_investigated"),
  "price_icpt": icpt, "price_slope": slope,
  "cum": [{"age": a["age"], "uc": a2["ca_death"], "g": a["ca_death"]} for a, a2 in zip(G, UC)],
  "psa": D["psa"],
  "scenarios": [{"s": k, "icer": r["icer"], "dq": r["inc_effect"]} for k, r in sc.items()],
  "subgroups": [{"age": r["start_age"], "sex": r["sex"], "race": r["race_eth"], "label": r["race_label"], "icer": r["icer_vs_usual_care"],
                 "dq": r["inc_qaly"], "averted": r["cancer_deaths_averted"]} for r in D["subgroups"]],
  "kunz": {"egd": v("egd", "Kunzmann -> EGD"), "perf_deaths": v("perforation_deaths", "Kunzmann -> EGD"),
           "qaly_net": v("qaly_net", "Kunzmann -> EGD"), "averted": v("ca_deaths_averted", "Kunzmann -> EGD"),
           "cost": (v("cost_total", "Kunzmann -> EGD") - v("cost_total", "Usual care")) * 1e5},
}
burden = open(os.path.join(HERE, "burden_template.html")).read()
style = re.search(r"<style>(.*?)</style>", burden, re.S).group(1)
html = open(os.path.join(HERE, "compare_template.html")).read().replace("__CSS__", style).replace("__DATA__", json.dumps(data, separators=(",", ":")))
open(os.path.join(HERE, "galleri_vs_usual_care.html"), "w").write(html)
print("wrote galleri_vs_usual_care.html", len(html), round(icpt), round(slope, 1), round((100000 - icpt) / slope, 2), round((150000 - icpt) / slope, 2))
