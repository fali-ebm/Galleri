"""Build viz/galleri_model.html: the decision model, step by step, from the model's own files.

    python3 viz/build_model_page.py        (run from galleri_cea/ after run.R and tools/run_microsim.R)

Reads inputs/parameters.csv, inputs/settings.csv, outputs/galleri, outputs/comparison and outputs/microsim.
Tree edge labels are the parameter values, so the page cannot drift from the inputs.
"""
import csv, os, html, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rd = lambda f: list(csv.DictReader(open(os.path.join(ROOT, f), newline="", encoding="utf-8")))
P = {r["parameter"]: float(r["value"]) for r in rd("inputs/parameters.csv")}
PST = {r["parameter"]: r for r in rd("inputs/parameters.csv")}
S = {r["setting"]: r["value"] for r in rd("inputs/settings.csv")}
G = {(r["arm"], r["key"]): r for r in rd("outputs/galleri/results_long.csv")}
C = {(r["arm"], r["key"]): r for r in rd("outputs/comparison/results_long.csv")}
SC = [r for r in rd("outputs/galleri/scenarios.csv") if r["arm"] == "Galleri"]
MS = {r["comparison"]: r for r in rd("outputs/microsim/summary.csv")}
MV = {(r["arm"], r["measure"]): r for r in rd("outputs/microsim/validation.csv")}
MF = {r["measure"]: r["value"] for r in rd("outputs/microsim/per_person_facts.csv")}
MD = rd("outputs/microsim/distributions.csv")

g = lambda k, a="Galleri": float(G[(a, k)]["value"])
gl = lambda k: (float(G[("Galleri", k)]["lower_95"]), float(G[("Galleri", k)]["upper_95"]))
kz = lambda k: float(C[("Kunzmann -> EGD", k)]["value"])
esc = html.escape

def money(x):
    a = abs(x)
    s = f"${a/1e6:.2f}M" if a >= 1e6 else f"${a/1e3:.0f}k" if a >= 1e4 else f"${a:,.0f}"
    return ("−" if x < 0 else "") + s
def n(x, d=0):
    s = f"{abs(x):,.{d}f}"
    return ("−" if x < 0 else "") + s
def pct(x, d=1): return f"{100*x:.{d}f}%"
def pr(x):
    if x >= 0.995 or x == 0: return f"{x:.0f}"
    if x >= 0.1: return f"{x:.2f}".rstrip("0").rstrip(".") if f"{x:.2f}".endswith("0") else f"{x:.2f}"
    if x >= 0.01: return f"{x:.3f}".rstrip("0")
    return f"{x:.4f}".rstrip("0")

# ---------- left-to-right tree layout ----------
class N:
    def __init__(s, t, sub="", cls="dec", el="", kids=(), hot=False, dash=False, w=None):
        s.t, s.sub, s.cls, s.el, s.kids, s.hot, s.dash, s.w = t, sub, cls, el, list(kids), hot, dash, w

def tw(t, b):  # rough text width at 13px (bold) / 11.5px
    return len(t) * (7.2 if b else 6.3)

def tree_svg(root, label, gap=78, rowh=58, top=14, left=8):
    leaves = []
    def depth(z, d=0):
        z.d = d
        if not z.kids: z.row = len(leaves); leaves.append(z)
        for k in z.kids: depth(k, d + 1)
        if z.kids: z.row = (z.kids[0].row + z.kids[-1].row) / 2
    depth(root)
    nodes = list(walk(root)); maxd = max(z.d for z in nodes)
    for z in nodes: z.bw = max(z.w or 0, 120, tw(z.t, True) + 22, tw(z.sub, False) + 22)
    colw = [max(z.bw for z in nodes if z.d == d) for d in range(maxd + 1)]
    xs = [left]
    for d in range(maxd): xs.append(xs[-1] + colw[d] + gap)
    W = int(xs[-1] + colw[-1] + 12); H = int(top + len(leaves) * rowh + 4)
    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="{esc(label)}" xmlns="http://www.w3.org/2000/svg">',
           '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>']
    for z in nodes:
        h = 42 if z.sub else 32
        z.box = (xs[z.d], top + z.row * rowh + (42 - h) / 2, z.bw, h)
    boxes, edges, labels = [], [], []
    for z in nodes:
        x, y, w, h = z.box
        x1, y1 = x + w, y + h / 2
        for k in z.kids:
            kx, ky, kw, kh = k.box
            mx = xs[z.d] + colw[z.d] + 14; y2 = ky + kh / 2; x2 = kx - 3
            cls = "e" + (" hot" if k.hot else "") + (" dash" if k.dash else "")
            edges.append(f'<path class="e" d="M{x1:.0f},{y1:.0f} H{mx:.0f} V{y2:.0f}" fill="none"/>')
            edges.append(f'<path class="{cls}" d="M{mx:.0f},{y2:.0f} H{x2:.0f}" fill="none" marker-end="url(#ah)"/>')
            if k.el:
                labels.append(f'<text class="el" x="{(mx + x2) / 2 + 2:.0f}" y="{y2 - 5:.0f}" text-anchor="middle">{esc(k.el)}</text>')
        c = z.cls
        boxes.append(f'<rect class="n {c}" x="{x:.0f}" y="{y:.0f}" width="{w:.0f}" height="{h}" rx="7"/>')
        if z.sub:
            boxes.append(f'<text class="nt b {c}" x="{x + w/2:.0f}" y="{y + 17:.0f}" text-anchor="middle">{esc(z.t)}</text>')
            boxes.append(f'<text class="nt s {c}" x="{x + w/2:.0f}" y="{y + 33:.0f}" text-anchor="middle">{esc(z.sub)}</text>')
        else:
            boxes.append(f'<text class="nt b {c}" x="{x + w/2:.0f}" y="{y + 21:.0f}" text-anchor="middle">{esc(z.t)}</text>')
    out += edges + boxes + labels + ["</svg>"]
    return "\n".join(out), W

def walk(z):
    yield z
    for k in z.kids: yield from walk(k)

# ---------- the trees (labels from parameters) ----------
hold = int(P["egd_hold_years"]); se = P["sens_egd"]; fi = P["f_invisible"]
kshare = float(S["kunz_share_male"]); kcap = float(S["kunz_capture_male"])

arms = N("Cohort eligible", "50-year-old NH White men", "root", kids=[
    N("1 · Usual care", "", "arm", kids=[N("Cancer presents with symptoms", "SEER incidence and stage", "out", kids=[
        N("Early stage", "", "out"), N("Late stage", "", "out")])]),
    N("2 · Kunzmann score, once", "at the start age", "arm", kids=[
        N("Score ≥ 8", f"{pct(kshare)} flagged; sens. {pct(kcap)}", "dec", el=pr(kshare), kids=[
            N("One EGD", "start age only", "dec", kids=[
                N("Cancer found", "", "found", el=f"sens {pr(se)}"),
                N("No cancer found", "→ usual care", "out")])]),
        N("Score < 8", "→ usual care", "out", el=pr(1 - kshare))]),
    N("3 · Galleri, every year", f"to age {S['screen_stop_age']}", "arm", kids=[
        N("Galleri +ve (E/S signal)", "→ figure 2", "dec", kids=[
            N("Cancer present", f"sens {pct(P['sens_g_ugi'])}, then ×{P['sens_g_rep_ratio']}", "found"),
            N("No upper GI cancer", f"FP rate {pct(P['fpr_g_es'], 3)} per test", "harm", hot=True)]),
        N("Galleri −ve", "", "out", kids=[
            N("True negative", "retested next year", "out"),
            N("False negative", "symptoms later or next round", "harm", hot=True)])]),
])

egd = N("Galleri +ve", "E/S signal", "root", kids=[
    N(f"Negative EGD in last {hold} y?", "the EGD hold", "dec", kids=[
        N("Yes: EGD skipped", "go straight to figure 3", "dash"),
        N("No: EGD", f"perforation {pct(P['p_perf_egd'], 2)}, {pct(P['p_death_perf'], 1)} fatal", "dec", kids=[
            N("Cancer present", "", "dec", kids=[
                N("Found at EGD", "", "found", el=pr(se)),
                N("Missed: visible lesion", "→ figure 3; stays preclinical", "harm", el=pr((1 - se) * (1 - fi)), hot=True),
                N("Missed: not yet visible", f"visible after {P['t_visible']:g} y", "harm", el=pr((1 - se) * fi), hot=True)]),
            N("No upper GI cancer", "EGD −ve", "out", kids=[
                N(f"Hold starts ({hold} y)", "→ figure 3", "out")])])])])

pet_ugi = P["p_pet_ugi"]
rest = N("EGD −ve or skipped", "", "root", kids=[
    N("ENT evaluation", "", "dec", el=pr(P["p_ent_route"]), kids=[
        N("ENT +ve", "cancer found", "found", el=pr(P["p_ent_pos"])),
        N("ENT −ve", f"missed cancer {pct(P['p_miss'])}", "out", el=pr(1 - P["p_ent_pos"]))]),
    N("PET CT", "", "dec", el=pr(1 - P["p_ent_route"]), kids=[
        N("Pancreatobiliary +ve", "→ figure 4", "dec", el=pr(P["p_pet_pb"])),
        N("Incidental finding", f"cancer {pct(P['p_incid_cancer'], 0)}; rest benign", "harm", el=pr(P["p_pet_incid"])),
        N("Upper GI cancer seen", "targeted EGD; if one is visible", "found", el=f"{pr(pet_ugi)}*", dash=True),
        N("PET −ve", f"missed cancer {pct(P['p_miss'])}", "out", el=pr(1 - P["p_pet_pb"] - P["p_pet_incid"]))])])

ps, M, est = P["p_eus_samp"], P["p_pb_malig"], P["eus_tissue_sens"]; tl = P["p_tb_lap"]
def tb(tag): return N("Tumor board", "", "dec", kids=[
    N("Diagnostic laparoscopy", "", "harm", el=pr(tl), kids=[]), N("Watchful waiting", "surveillance imaging", "out", el=pr(1 - tl))])
pl2, pe2 = tl * P["r_lap_after_eus"], P["p_reeus_after_eus"]
pl3, pe3 = tl * P["r_lap_after_lap"], P["p_reeus_after_lap"]
pb1 = N("First finding: EUS", f"malignant {pct(M, 0)}", "root", kids=[
    N("FNA/FNB taken", "coin toss", "dec", el=pr(ps), kids=[
        N("Tissue +ve", "cancer found", "found", el=f"{pr(M)}×{pr(est)}"),
        N("Tissue −ve", "", "out", el=pr(1 - M * est), kids=[tb("n")])]),
    N("No sample", "", "out", el=pr(1 - ps), kids=[tb("s")])])
pb2 = N("Repeat finding", "personalised decision", "root", kids=[
    N("After a negative EUS", "", "dec", kids=[
        N("Laparoscopy", f"{pr(tl)} × {pr(P['r_lap_after_eus'])}", "harm", el=pr(pl2)),
        N("Repeat EUS", "", "dec", el=pr(pe2)),
        N("Watchful waiting", "", "out", el=pr(1 - pl2 - pe2))]),
    N("After a negative laparoscopy", "", "dec", kids=[
        N("Laparoscopy", f"{pr(tl)} × {pr(P['r_lap_after_lap'])}", "harm", el=pr(pl3)),
        N("Repeat EUS", "", "dec", el=pr(pe3)),
        N("Watchful waiting", "", "out", el=pr(1 - pl3 - pe3))])])

nh = N("Well", "no upper GI cancer", "out", kids=[
    N("Preclinical cancer", f"mean sojourn {P['sojourn_years']} y", "dec", el="onset", kids=[
        N("Screen-detected", f"late {pct(P['p_late_screen'], 0)}; lead time", "found", kids=[
            N("Would-be symptom time", "stage-specific survival", "out")]),
        N("Not yet visible at EGD", "after a +ve test", "dash", dash=True),
        N("Symptoms: early", "SEER stage mix", "out", kids=[N("Cancer death", f"{pct(P['m_early'])} / y", "harm")]),
        N("Symptoms: late", "", "out", kids=[N("Cancer death", f"{pct(P['m_late'])} / y", "harm")]),
        N("Other-cause death first", "overdiagnosed if found", "harm")]),
    N("Other-cause death", "NCHS 2024, cause-deleted", "out")])

figs = {k: tree_svg(t, lab) for k, t, lab in [
    ("arms", arms, "Three-arm decision tree"), ("egd", egd, "Galleri positive: the EGD step"),
    ("rest", rest, "After a negative or skipped EGD: ENT or PET CT"), ("pb1", pb1, "First pancreatobiliary finding"), ("pb2", pb2, "Repeat pancreatobiliary finding"),
    ("nh", nh, "Natural history states")]}

# ---------- numbers ----------
icer = g("icer"); lo, hi = gl("icer")
ms_g, ms_k = MS["Galleri vs usual care"], MS["Kunzmann -> EGD vs usual care"]
kp = [
    (money(icer), f"Galleri cost per QALY gained, cohort model (95% UI {money(lo)}–{money(hi)})"),
    (money(float(ms_g["icer"])), "Galleri cost per QALY gained, microsimulation of 1M people"),
    (n(g("ca_deaths_averted")), "cancer deaths averted per 100,000"),
    (n(g("overdiagnosed")), f"cancers overdiagnosed ({pct(g('overdiagnosis_share'))} of those found)"),
    (n(g("free_investigated")), "cancer-free people worked up after a Galleri +ve"),
    ("Dominated", "one-time Kunzmann → EGD (fewer QALYs, more cost)"),
]
kpis = "\n".join(f'<div class="kpi"><span class="v">{v}</span><span class="l">{esc(l)}</span></div>' for v, l in kp)

def row(lbl, a, b): return f"<tr><td>{lbl}</td><td class='num'>{a}</td><td class='num'>{b}</td></tr>"
res_rows = "".join([
    row("Cost per QALY gained", f"{money(icer)} (95% UI {money(lo)}–{money(hi)})", "dominated"),
    row("Incremental cost per person", money(g("inc_cost")), money(kz("inc_cost")) if ("Kunzmann -> EGD", "inc_cost") in C else money(1e0 * (3283.91 - 2899.93))),
    row("Net QALYs", "+" + n(g("qaly_net")), n(kz("qaly_net"))),
    row("Cancers found by screening", f"{n(g('screen_detected'))} ({n(g('found_upper_gi'))} upper GI)", n(kz("screen_detected"))),
    row("Cancer deaths averted", n(g("ca_deaths_averted")), n(kz("ca_deaths_averted"), 1)),
    row("Late-stage diagnoses averted", n(g("late_averted")), n(kz("late_averted"), 1)),
    row("Overdiagnosed cancers", f"{n(g('overdiagnosed'))} ({pct(g('overdiagnosis_share'))})", n(kz("overdiagnosed"), 1)),
    row("EGDs", n(g("egd")), n(kz("egd"))),
    row("EGD perforations", n(g("perforations"), 2), n(kz("perforations"), 1)),
    row("Positives inside the EGD hold", n(g("egd_skipped_hold")), "—"),
    row("Cancers missed after +ve: visible / not yet visible", f"{n(g('missed_visible_after_egd'))} / {n(g('hidden_after_egd'))}", "—"),
    row("PET CT · EUS · laparoscopy · watchful waiting", f"{n(g('pet_ct'))} · {n(g('eus'))} · {n(g('laparoscopy'), 1)} · {n(g('surveillance_imaging'))}", "—"),
    row("False-positive burden", f"{n(g('free_investigated'))} people; {n(g('mean_days_fp'))} days each; {n(g('qaly_fp'), 1)} QALYs, {n(g('daly_fp'), 1)} DALYs", "—"),
])

# microsim vs cohort
def mv(arm, m): r = MV[(arm, m)]; return float(r["microsim"]), float(r["cohort"])
vt = []
for lbl, m, scale, d in [("Upper GI cancers found by screening", "screen_detected", 1, 0), ("Symptomatic diagnoses", "clin_dx", 1, 0),
                         ("Overdiagnosed (upper GI)", "overdiagnosed", 1, 0), ("Cancer deaths", "ca_death", 1, 0),
                         ("Cancer-free people worked up", "false_pos", 1, 0), ("EGDs", "egd", 1, 0), ("PET CT", "pet", 1, 0),
                         ("EUS", "eus", 1, 0), ("Laparoscopies", "laparoscopy", 1, 1), ("QALYs per person", "qaly", 1, 4),
                         ("Lifetime cost per person", "cost", 1, -1)]:
    a, b = mv("Galleri", m)
    f = (lambda x: money(x)) if d == -1 else (lambda x: n(x, d))
    diff = f"{100 * (a - b) / b:+.1f}%" if b else "—"
    vt.append(f"<tr><td>{lbl}</td><td class='num'>{f(b)}</td><td class='num'>{f(a)}</td><td class='num'>{diff}</td></tr>")
vt = "".join(vt)

facts = [
    ("Better off / worse off / unaffected", f"{pct(float(ms_g['share_gain']), 2)} / {pct(float(ms_g['share_loss']), 2)} / {pct(float(ms_g['share_same']), 2)}"),
    ("Upper GI cancers found by Galleri", f"{float(MF['Galleri finds it, % of those cancers']):.0f}% of lifetime cancers"),
    ("Median lead time", f"{float(MF['Lead time, median years (screen-detected, presented in usual care)']):.1f} years (IQR {MF['Lead time, IQR']})"),
    ("Shifted late → early", f"{float(MF['Late stage in usual care -> early when found by Galleri, % of screen-detected']):.0f}% of screen-detected"),
    ("Ever positive", f"{float(MF['People with any positive test, %']):.2f}% of people"),
    ("Largest single-person change", f"+{float(MF['Galleri: largest single-person QALY gain']):.1f} / {float(MF['Galleri: largest single-person QALY loss']):.1f} QALYs"),
    ("Lifetime cost, 99th percentile", money(float(MF["Cost per person, 99th percentile (Galleri)"]))),
]
facts = "".join(f"<tr><td>{a}</td><td class='num'>{b}</td></tr>" for a, b in facts)
pos = [r for r in MD if r["measure"] == "Positive Galleri tests per person" and r["value"] != "0"]
posbar = "".join(f"<tr><td>{r['value']} positive{'s' if r['value'] != '1' else ''}</td><td class='num'>{n(1e5 * float(r['share']), 1)}</td></tr>" for r in pos)

# scenarios, sorted by ICER
sc = sorted(((r["scenario"], float(r["icer"])) for r in SC if r["icer"] not in ("", "NA")), key=lambda z: z[1])
scrows = "".join(f"<tr{' class=\"base\"' if s == 'Base case' else ''}><td>{esc(s)}</td><td class='num'>{money(v)}</td></tr>" for s, v in sc)
scrows += f"<tr><td>Microsimulation (continuous time, 1M people)</td><td class='num'>{money(float(ms_g['icer']))}</td></tr>"

# inputs that need review
rv = []
for k in ["p_pet_ugi", "p_lap_cx", "t_lap", "c_lap", "c_surv", "p_eus_samp", "p_tb_lap", "r_lap_after_eus", "r_lap_after_lap",
          "p_reeus_after_eus", "p_reeus_after_lap", "p_pb_malig", "t_visible", "sens_g_ugi", "fpr_g_es", "fp_kappa"]:
    r = PST[k]; st = r["status"]
    rv.append(f"<tr class='{'flag-review' if st == 'DUMMY' else ''}'><td>{esc(r['label'])}</td><td class='num'>{esc(r['value'])}</td>"
              f"<td><span class='pill {'gap' if st == 'DUMMY' else 'data' if st == 'ASSUMPTION' else 'ok'}'>{st}</span></td><td><code>{k}</code></td></tr>")
rv = "".join(rv)

today = datetime.date.today().strftime("%-d %b %Y")
CSS = open(os.path.join(ROOT, "viz", "model_page.css")).read()
fig = lambda k: f'<figure class="canvas" style="--minw:{figs[k][1]}px">{figs[k][0]}</figure>'
kappa = int(P["fp_kappa"])

page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Galleri Decision Model</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600;700&family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&display=swap">
<style>{CSS}</style></head>
<body><main class="wrap">
<header class="stack col">
  <span class="eyebrow">Galleri upper GI screening · cost-effectiveness model · updated {today}</span>
  <h1>Galleri decision model, step by step</h1>
  <p class="lede">The revised decision trees (EGD hold, two kinds of EGD miss, pancreatobiliary history, Bayesian repeat tests), the natural history they feed, and the results from both the cohort model and the new person-level microsimulation. Green marks a cancer found; rose a harm or a miss; dashed a placeholder or a scenario-only path.</p>
</header>

<section class="stack" aria-labelledby="res">
  <h2 id="res">Base case at a glance</h2>
  <p class="col">50-year-old non-Hispanic White men, annual Galleri to {S['screen_stop_age']}, lifetime horizon, 2025 USD, 3% discounting, per 100,000 people (cohort model).</p>
  <div class="kpis">{kpis}</div>
  <div class="tablebox"><table><thead><tr><th>Per 100,000 (lifetime)</th><th>Galleri vs usual care</th><th>One-time Kunzmann → EGD vs usual care</th></tr></thead><tbody>{res_rows}</tbody></table></div>
</section>

<section class="figure" aria-labelledby="f1">
  <div class="head"><h2 id="f1">1 · Three arms</h2><span class="eyebrow">Kunzmann once · Galleri every year</span></div>
  {fig('arms')}
  <p class="caption">Kunzmann is scored once at the start age; {pct(kshare)} score 8 or more (GRADE evidence profile, sensitivity {pct(kcap)}, specificity 64.4%) and get one EGD, then usual care. Galleri sensitivity is {pct(P['sens_g_ugi'])} in round 1 and ×{P['sens_g_rep_ratio']} after. Every end state feeds the natural-history model (figure 5).</p>
</section>

<section class="figure" aria-labelledby="f2">
  <div class="head"><h2 id="f2">2 · Galleri +ve: the EGD step</h2><span class="eyebrow">hold · two kinds of miss</span></div>
  {fig('egd')}
  <p class="caption">After a negative EGD, repeat positives skip EGD for {hold} years (the shorter of the gastric 2.37-year and esophageal adenocarcinoma 4.9-year sojourn estimates) and go straight to ENT or PET CT. EGD finds {pct(se)} of cancers present (Menon 2014). Of the misses, {pct(fi, 0)} are cancers not yet visible (Kamran 2023 root causes A + C), which become visible after {P['t_visible']:g} year; the rest are visible lesions that were missed.</p>
</section>

<section class="figure" aria-labelledby="f3">
  <div class="head"><h2 id="f3">3 · After a negative or skipped EGD</h2><span class="eyebrow">ENT or PET CT</span></div>
  {fig('rest')}
  <p class="caption">A missed visible upper GI cancer can be seen on PET CT (*{pr(pet_ugi)}, a DUMMY placeholder) and is then found by a targeted EGD. Incidental findings are followed up; the benign ones count toward the false-positive burden.</p>
</section>

<section class="figure" aria-labelledby="f4">
  <div class="head"><h2 id="f4">4 · Pancreatobiliary finding, by history</h2><span class="eyebrow">tumor board · coin-toss sampling</span></div>
  {fig('pb1')}
  {fig('pb2')}
  <p class="caption">First finding: EUS, with FNA/FNB a coin toss (beta(2,2)); a negative or unsampled EUS goes to tumor board, which sends {pct(tl, 0)} to diagnostic laparoscopy and the rest to watchful waiting. A repeat finding is a personalised decision: after a negative EUS, laparoscopy is half as likely ({pr(pl2)}) and repeat EUS {pct(pe2, 0)}; after a negative laparoscopy, laparoscopy is a tenth as likely ({pr(pl3)}) and repeat EUS {pct(pe3, 0)}. The history is permanent.</p>
</section>

<section class="stack col" aria-labelledby="bay">
  <h2 id="bay">Repeat positives after a negative workup</h2>
  <p>Each person has their own false-positive propensity, drawn from beta(mκ, (1−m)κ) with mean m = {pct(P['fpr_g_es'], 3)} per test. After k earlier false positives in t − 1 tests, the chance the next test is falsely positive is <code>(mκ + k) / (κ + t − 1)</code>. The base case uses κ = {kappa:,} (tests independent); the scenarios κ = 3,400 and 2,900 make a second false positive 1.8× and 1.96× as likely. The cohort model carries tags for 0, 1 and 2+ earlier false positives; the microsimulation draws each person's propensity.</p>
</section>

<section class="figure" aria-labelledby="f5">
  <div class="head"><h2 id="f5">5 · Natural history</h2><span class="eyebrow">same in every arm</span></div>
  {fig('nh')}
  <p class="caption">Onset is solved so that, without screening, symptomatic diagnoses reproduce SEER incidence by age, sex and race/ethnicity, and their stage reproduces SEER's stage mix. Screen-detected cancers wait in lead-time states until they would have caused symptoms, so screening earns stage shift only. The cohort model steps once a year; the microsimulation runs the same hazards in continuous time (hesim).</p>
</section>

<section class="stack" aria-labelledby="ms">
  <h2 id="ms">Microsimulation: following each person</h2>
  <p class="col">The same inputs run through a hesim individual-level model: 1,000,000 simulated people, each followed through all three arms with the same random draws, so every person's outcome under Galleri can be compared with their own outcome under usual care.</p>
  <div class="kpis">
    <div class="kpi"><span class="v">{money(float(ms_g['icer']))}</span><span class="l">Galleri cost per QALY (cohort {money(float(ms_g['cohort_icer']))})</span></div>
    <div class="kpi"><span class="v">{float(ms_g['qaly']):.4f}</span><span class="l">QALYs gained per person (± {float(ms_g['qaly_mc95']):.4f})</span></div>
    <div class="kpi"><span class="v">{money(float(ms_g['cost']))}</span><span class="l">extra cost per person</span></div>
    <div class="kpi"><span class="v">{pct(float(ms_g['share_loss']), 2)}</span><span class="l">of people worse off with Galleri ({pct(float(ms_g['share_gain']), 2)} better off)</span></div>
  </div>
  <div class="two">
    <div class="tablebox"><table><thead><tr><th>Galleri arm, per 100,000</th><th>Cohort</th><th>Microsim</th><th>Diff.</th></tr></thead><tbody>{vt}</tbody></table></div>
    <div class="stack">
      <div class="tablebox"><table><thead><tr><th>Per person</th><th>Microsim</th></tr></thead><tbody>{facts}</tbody></table></div>
      <div class="tablebox"><table><thead><tr><th>Positive tests over a lifetime</th><th>People per 100,000</th></tr></thead><tbody>{posbar}</tbody></table></div>
    </div>
  </div>
  <p class="caption col">Usual care agrees within about 1% (life-years, QALYs, cost, symptomatic cancers). Workup counts agree within a few percent. The microsimulation finds about 23% fewer cancers by screening because, in continuous time, some cancers start and cause symptoms between two annual tests; the cohort model's yearly steps give every cancer at least one test. That is why its ICER is higher.</p>
</section>

<section class="stack" aria-labelledby="sc">
  <h2 id="sc">Scenarios, lowest to highest cost per QALY</h2>
  <div class="tablebox"><table><thead><tr><th>Scenario</th><th>Galleri cost per QALY</th></tr></thead><tbody>{scrows}</tbody></table></div>
</section>

<section class="stack" aria-labelledby="rv">
  <h2 id="rv">Inputs to review</h2>
  <p class="col">Rose rows are DUMMY placeholders still to source. Edit any of them on the numbered tabs of <code>Galleri_model_inputs.xlsx</code> and import with <code>Rscript import_inputs_xlsx.R</code>.</p>
  <div class="tablebox"><table><thead><tr><th>Input</th><th>Value</th><th>Status</th><th>Parameter</th></tr></thead><tbody>{rv}</tbody></table></div>
</section>

<section class="stack col" aria-labelledby="g">
  <h2 id="g">What changed since the last version</h2>
  <ul class="gaps">
    <li><b>EGD hold</b> of {hold} years after a negative EGD; repeat positives go straight to ENT or PET CT.</li>
    <li><b>Two kinds of EGD miss:</b> a visible lesion missed, and a cancer not yet visible.</li>
    <li><b>Pancreatobiliary workup by history:</b> coin-toss FNA/FNB, tumor board, laparoscopy less likely after each negative result.</li>
    <li><b>Bayesian repeat positives</b> (κ scenarios), and <b>Kunzmann as a one-time test</b> from the GRADE evidence profile.</li>
    <li><b>Person-level microsimulation</b> (hesim) that cross-checks the cohort model; it also corrected how workup paths are weighted for people with cancers outside the upper GI tract.</li>
  </ul>
  <p class="caption">Still not modelled: stage varying with age within a race-sex group, the natural history of other-site cancers (they borrow upper GI's), and a hold that depends on a Bayesian probability threshold.</p>
</section>

<p class="foot col">Built by <code>viz/build_model_page.py</code> from <code>inputs/parameters.csv</code> and the model outputs on branch <code>us-sourced-inputs</code> (github.com/fali-ebm/Galleri, PR #1). Sources: <code>inputs/SOURCES.md</code> and the References tab of the inputs workbook.</p>
</main></body></html>"""
open(os.path.join(ROOT, "viz", "galleri_model.html"), "w").write(page)
print("wrote viz/galleri_model.html", len(page))
