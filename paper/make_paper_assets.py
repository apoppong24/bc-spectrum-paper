#!/usr/bin/env python3
"""Generate the LaTeX table snippets and the two remaining figures for the
paper from results/*.json, so no number is transcribed by hand.

    python paper/make_paper_assets.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
RES = ROOT / "results"
OUT = ROOT / "paper"
(OUT / "tables").mkdir(parents=True, exist_ok=True)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

from quarkonium import pdg  # noqa: E402

fits = {p.stem[4:]: json.loads(p.read_text()) for p in RES.glob("fit_*.json")}
paper = json.loads((RES / "paper_results.json").read_text())
MAIN = ["martin", "log", "cornell"]
NAME = {"martin": "Martin", "log": "Log", "cornell": "Cornell"}


def tex_state(lab: str) -> str:
    r"""'13P2' -> $1\,^3P_2$ ; '1P1'' -> $1^{+'}$ (mixed)."""
    prime = "'" if lab.endswith("'") else ""
    lab = lab.rstrip("'")
    if len(lab) == 3:
        return rf"$1^{{+{prime}}}$" if lab[1] == "P" else rf"${lab[0]}{lab[1]}_{{{lab[2]}}}{prime}$"
    n, s, L, J = lab
    return rf"${n}\,^{s}{L}_{J}$"


# ---------------------------------------------------------------- Table: data set
rows = []
for pair, table, name in (("cc", pdg.CHARMONIUM, r"$c\bar c$"), ("bb", pdg.BOTTOMONIUM, r"$b\bar b$")):
    for lab, (m, s) in table.items():
        rows.append(f"{name} & {tex_state(lab)} & ${m:.1f} \\pm {s:.1f}$ \\\\")
(OUT / "tables/dataset.tex").write_text("\n".join(rows) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Table: fitted parameters
PNAME = {"B": "B", "alpha": r"\alpha", "m_c": "m_c", "m_b": "m_b", "alpha_s": r"\alpha_s",
         "C": "C", "kappa": r"\kappa", "a": "a", "Lambda": r"\Lambda"}


def fmt(n, x, e):
    if e < 0.1:
        dig = 3
        return f"${PNAME[n]} = {x:.{dig}f}({round(e * 10**dig):.0f})$"
    return f"${PNAME[n]} = {x:.2f} \\pm {e:.2f}$"
lines = []
for fam in MAIN + ["martin_run_geometric", "log_run_geometric"]:
    d = fits[fam]
    params = ", ".join(fmt(n, x, e) for n, x, e in zip(d["names"], d["x"], d["errors"]))
    label = {**NAME, "martin_run_geometric": r"Martin, running $\alpha_s$", "log_run_geometric": r"Log, running $\alpha_s$"}[fam]
    lines.append(f"{label} & {d['chi2']:.1f}/{d['dof']} & {params} \\\\")
(OUT / "tables/params.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Table: Bc spectrum
spec = paper["spectrum"]
order = ["11S0", "13S1", "13P0", "1P1", "1P1'", "13P2", "21S0", "23S1",
         "13D1", "1D2", "1D2'", "13D3", "23P0", "2P1", "2P1'", "23P2", "31S0", "33S1"]
EQ94 = {"11S0": 6264, "13S1": 6337, "13P0": 6700, "1P1": 6730, "1P1'": 6736, "13P2": 6747,
        "21S0": 6856, "23S1": 6899, "13D1": 7012, "1D2": 7009, "1D2'": 7012, "13D3": 7005,
        "23P0": 7108, "2P1": 7135, "2P1'": 7142, "23P2": 7153, "31S0": 7244, "33S1": 7280}
lines = []
for lab in order:
    s = spec[lab]
    cells = [f"${s[f]:.0f}({s['param_err'][f]:.0f})$" for f in MAIN]
    pdg = pdg.BC.get(lab)
    lines.append(f"{tex_state(lab)} & " + " & ".join(cells)
                 + f" & ${s['mean']:.0f} \\pm {s['spread']:.0f}$ & {EQ94.get(lab, '--')} & "
                 + (f"${pdg[0]:.1f} \\pm {pdg[1]}$" if pdg else "--") + r" \\")
(OUT / "tables/bc_spectrum.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Table: hyperfine
hyp = paper["hyperfine"]
lines = []
for fam, label in [("martin", "Martin"), ("log", "Log"), ("cornell", "Cornell"),
                   ("martin_run_geometric", r"Martin, running ($\sqrt{m_cm_b}$)"),
                   ("martin_run_reduced", r"Martin, running ($2\mu$)"),
                   ("log_run_geometric", r"Log, running ($\sqrt{m_cm_b}$)"),
                   ("log_run_reduced", r"Log, running ($2\mu$)")]:
    h = hyp[fam]
    lines.append(f"{label} & {h['d1S']:.1f} & {h['d2S']:.1f} & {h['diff']:.1f} & {h['bcstar']:.0f} \\\\")
(OUT / "tables/hyperfine.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Table: scan
lines = []
for s in paper["scan"]:
    lines.append(f"{s['alpha']:g} & {s['chi2']:.0f}/{s['dof']} & {s['alpha_s']:.3f} & {s['11S0']:.0f} & {s['13S1']:.0f} & "
                 f"{s['21S0']:.0f} & {s['P1']:.0f} & {s['D1']:.0f} & {s['33S1']:.0f} & "
                 f"{s['13S1']-s['11S0']-(s['23S1']-s['21S0']):.1f} \\\\")
(OUT / "tables/scan.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Table: radiative
rad = sorted(paper["radiative"], key=lambda r: -r["width_keV"])
keep = [r for r in rad if r["width_keV"] >= 10 and r["type"] == "E1"] + [r for r in rad if r["type"] == "M1"]
lines = [f"{tex_state(r['initial'])} $\\to$ {tex_state(r['final'])}$\\,\\gamma$ & {r['type']} & {r['k_MeV']:.0f} & {r['width_keV']:.1f} \\\\" for r in keep]
(OUT / "tables/radiative.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Table: uncertainty budget
bud = json.loads((RES / "uncertainty_budget.json").read_text())
b = bud["budget"]["martin_run_geometric"]
srcs = [k for k in b if not k.startswith("_")]
lines = []
for lab in ["11S0", "13S1", "21S0", "23S1", "13P0", "1P1", "1P1'", "13P2", "13D1", "13D3", "33S1", "d1S", "d2S"]:
    name = tex_state(lab) if lab[0].isdigit() else {"d1S": r"$\Delta_{1S}$", "d2S": r"$\Delta_{2S}$"}[lab]
    lines.append(f"{name} & {b['_central'][lab]:.0f} & " + " & ".join(f"{b[k][lab]:.0f}" for k in srcs) + f" & {b['_total'][lab]:.0f} \\\\")
(OUT / "tables/budget.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")
(OUT / "tables/budget_header.tex").write_text(" & ".join({"floor": "floor", "data": "data set", "scale": r"$\alpha_s$ scale", "model": "model", "fit": "fit"}[k] for k in srcs) + "\n")

# ---------------------------------------------------------------- Table: self-consistent (excluded)
lines = []
NAMES = {"martin_sc": "Martin, fixed", "log_sc": "Log, fixed", "martin_run_geometric_sc": "Martin, running",
         "log_run_geometric_sc": "Log, running", "martin_run_geometric_sc_s0.8": r"Martin, running, $\sigma=0.8$",
         "martin_run_geometric_sc_s1.2": r"Martin, running, $\sigma=1.2$"}
for k, v in bud["selfconsistent"].items():
    lines.append(f"{NAMES[k]} & {v['chi2']:.0f}/{v['dof']} & {v['m_c']:.2f} & {v['cc_d1S']:.0f} & {v['bb_d1S']:.0f} & {v['d1S']:.0f} & {v['d2S']:.0f} \\\\")
(OUT / "tables/selfconsistent.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Table: experimental predictions
ex = json.loads((RES / "experimental_predictions.json").read_text())
lines = []
for lab in ["13P0", "1P1", "1P1'", "13P2", "21S0", "23S1", "13D1", "1D2", "1D2'", "13D3", "23P0", "2P1", "2P1'", "23P2", "31S0", "33S1"]:
    st = ex["states"][lab]
    tot = b["_total"].get(lab, None)
    err = f"{tot:.0f}" if tot else f"{np.hypot(st['fit_err'], st['model_err']):.0f}"
    ln = "; ".join(f"{tex_state(r['final'])} ({r['k_MeV']:.0f}, {r['width_keV']:.0f}, {100*r['BR_rad']:.0f}\\%)" for r in st["lines"][:3] if r["width_keV"] >= 1) or "--"
    lines.append(f"{tex_state(lab)} & ${st['mass']:.0f}\\pm{err}$ & {st['v2_c']:.2f} & {ln} \\\\")
(OUT / "tables/predictions.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Table: decay constants + v2
lines = []
for pair, name, lat, ref in (("cc", r"$\eta_c$", "394.7(2.4)", r"\cite{Davies2010}"), ("bb", r"$\eta_b$", "667(6)", r"\cite{McNeile2012}"), ("cb", r"$B_c$", "427(6)", r"\cite{McNeile2012}")):
    lines.append(f"{name} & {ex['decay_constants'][pair]:.0f} & {lat} {ref} \\\\")
(OUT / "tables/decay_constants.tex").write_text("\n".join(lines) + "\n\\bottomrule\n")

# ---------------------------------------------------------------- Figures
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 7.5,
                     "font.family": "serif", "mathtext.fontset": "cm", "axes.linewidth": 0.6,
                     "xtick.direction": "in", "ytick.direction": "in", "pdf.fonttype": 42})

# Fig 1: residuals
fig, ax = plt.subplots(figsize=(7.0, 3.2))
labels = [(r["pair"], r["label"]) for r in fits["martin"]["residuals"]]
x = np.arange(len(labels))
for i, fam in enumerate(MAIN + ["martin_run_geometric"]):
    res = {(r["pair"], r["label"]): 1000 * (r["model"] - r["exp"]) for r in fits[fam]["residuals"]}
    ax.plot(x + 0.18 * (i - 1.5), [res[k] for k in labels], "osd^"[i], ms=4, mfc="none",
            label={"martin_run_geometric": r"Martin, running $\alpha_s$"}.get(fam, NAME.get(fam)))
ax.axhspan(-10, 10, color="0.85", zorder=0)
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels([f"{p} {l}" for p, l in labels], rotation=90, fontsize=6.5)
ax.set_ylabel("model $-$ PDG [MeV]"); ax.legend(ncol=4, loc="upper left")
fig.tight_layout(); fig.savefig(OUT / "figures/fig1_residuals.pdf"); plt.close(fig)

# Fig 2: chi2 vs alpha and hyperfine difference vs alpha
scan = paper["scan"]
a = [s["alpha"] for s in scan]
fig, ax1 = plt.subplots(figsize=(3.4, 2.8))
ax1.plot(a, [s["chi2"] / s["dof"] for s in scan], "o-", color="C0")
ax1.set_xlabel(r"confining power $\alpha$"); ax1.set_ylabel(r"$\chi^2/\mathrm{dof}$", color="C0")
ax1.set_yscale("log")
ax2 = ax1.twinx()
ax2.plot(a, [s["13S1"] - s["11S0"] - (s["23S1"] - s["21S0"]) for s in scan], "s--", color="C3")
ax2.axhspan(29.1 - 1.7, 29.1 + 1.7, color="C3", alpha=0.2)
ax2.set_ylabel(r"$\Delta_{1S}-\Delta_{2S}$ [MeV]", color="C3")
fig.tight_layout(); fig.savefig(OUT / "figures/fig2_chi2_alpha.pdf"); plt.close(fig)

# Fig 3: band plot (re-render in pdf)
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
for key, name, c in (("11S0", "$B_c(1S)$", "C0"), ("13S1", "$B_c^*(1S)$", "C1"), ("21S0", "$B_c(2S)$", "C2"), ("23S1", "$B_c^*(2S)$", "C3")):
    axes[0].errorbar(a, [s[key] for s in scan], yerr=[s["err"][key] for s in scan], fmt="o-", ms=3, capsize=2, color=c, label=name)
axes[0].axhline(pdg.BC["11S0"][0], color="C0", ls=":", lw=1)
axes[0].axhline(pdg.BC["21S0"][0], color="C2", ls=":", lw=1)
axes[0].axhspan(6339.0 - 1.7, 6339.0 + 1.7, color="C1", alpha=0.35)
axes[0].set_xlabel(r"$\alpha$"); axes[0].set_ylabel("mass [MeV]"); axes[0].legend(ncol=2)
for key, name in (("P1", "1P centroid"), ("D1", "1D centroid"), ("33S1", r"$3\,^3S_1$")):
    axes[1].plot(a, [s[key] for s in scan], "o-", ms=3, label=name)
axes[1].set_xlabel(r"$\alpha$"); axes[1].legend()
fig.tight_layout(); fig.savefig(OUT / "figures/fig3_band.pdf"); plt.close(fig)

# Fig 4: hyperfine test
fig, ax = plt.subplots(figsize=(3.4, 2.8))
fams = ["martin", "log", "cornell", "martin_run_geometric", "martin_run_reduced", "log_run_geometric", "log_run_reduced"]
short = ["Martin", "Log", "Cornell", r"Martin$^{\rm run}_{\sqrt{m_cm_b}}$", r"Martin$^{\rm run}_{2\mu}$", r"Log$^{\rm run}_{\sqrt{m_cm_b}}$", r"Log$^{\rm run}_{2\mu}$"]
ax.plot(range(len(fams)), [hyp[f]["diff"] for f in fams], "s", mfc="none", color="C0", label=r"$\Delta_{1S}-\Delta_{2S}$")
ax.axhspan(29.1 - 1.7, 29.1 + 1.7, color="C1", alpha=0.3, label="CMS")
ax.axhline(30.9, color="C3", lw=1, label="LHCb")
ax.set_xticks(range(len(fams))); ax.set_xticklabels(short, rotation=60, fontsize=6.5)
ax.set_ylabel("[MeV]"); ax.legend(loc="upper right")
fig.tight_layout(); fig.savefig(OUT / "figures/fig4_hyperfine.pdf"); plt.close(fig)

# Fig 5: absolute 1S hyperfine vs ATLAS and lattice
fig, ax = plt.subplots(figsize=(3.4, 2.8))
ax.plot(range(len(fams)), [hyp[f]["d1S"] for f in fams], "o", mfc="none", color="C0", label=r"$\Delta_{1S}$")
ax.axhspan(64.5 - 1.7, 64.5 + 1.7, color="C1", alpha=0.35, label="ATLAS 2026")
ax.axhspan(54 - 3, 54 + 3, color="C2", alpha=0.3, label="HPQCD 2012")
ax.set_xticks(range(len(fams))); ax.set_xticklabels(short, rotation=60, fontsize=6.5)
ax.set_ylabel("[MeV]"); ax.legend(loc="upper right")
fig.tight_layout(); fig.savefig(OUT / "figures/fig5_hyperfine_abs.pdf"); plt.close(fig)
print("assets written to", OUT)
