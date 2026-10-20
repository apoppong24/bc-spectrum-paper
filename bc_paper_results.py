#!/usr/bin/env python3
"""Assemble the paper-level results from the refits:

  1. B_c spectrum from every fitted family with parametric (fit) errors
     and the model spread.
  2. The hyperfine test: Δ(1S) − Δ(2S) vs CMS (29.1 ± 1.7) and LHCb (30.9),
     and B_c*(1S) vs lattice (6330 ± 9).
  3. Exponent scan: B_c levels as a function of the confining power α.
  4. Radiative E1/M1 widths for B_c (Martin refit).

Writes results/paper_results.md, results/paper_results.json and
results/fig_band_exponent.png, results/fig_hyperfine_test.png.
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np

from quarkonium import pdg
from quarkonium.fit import FAMILIES, register_power_family
from quarkonium.radiative import radiative_table
from quarkonium.spectrum import spectrum
from quarkonium.spin_dependent import full_spectrum
from quarkonium.uncertainty import predict_with_errors

OUT = Path(__file__).parent / "results"
BC_LABELS = ["11S0", "13S1", "13P0", "1P1", "1P1'", "13P2", "21S0", "23S1",
             "13D1", "1D2", "1D2'", "13D3", "23P0", "2P1", "2P1'", "23P2", "31S0", "33S1"]
CMS_DIFF, CMS_ERR = 29.1, 1.7          # MeV, CMS 1902.00571
LHCB_DIFF = 6872.1 - 6841.2             # MeV, LHCb 1904.00081 (peak difference)
LATTICE_BCSTAR, LATTICE_ERR = 6330.0, 9.0   # HPQCD 0909.4462 (mass)
LATTICE_D1S, LATTICE_D1S_ERR = 54.0, 3.0     # HPQCD Dowdall et al. 2012 (hyperfine)


def load_all():
    fits = {}
    for p in sorted(glob.glob(str(OUT / "fit_*.json"))):
        d = json.loads(Path(p).read_text())
        if d["family"].startswith("power_"):
            register_power_family(d["alpha"])
        fits[d["family"]] = d
    return fits


def main():
    fits = load_all()
    print("families:", list(fits))
    res = {}
    for fam, d in fits.items():
        print("predicting", fam)
        res[fam] = predict_with_errors(fam, d["x"], d["cov"], "cb", BC_LABELS)

    def alpha_of(fam):
        return {"log": 0.0, "martin": 0.1}.get(fam, float(fam.split("_")[1]) if fam.startswith("power_") else None)

    power_fams = sorted([f for f in fits if alpha_of(f) is not None], key=alpha_of)
    main_fams = ["martin", "log", "cornell"]

    md = ["# B_c from charmonium + bottomonium: paper results\n"]

    # ---- 1. spectrum table -------------------------------------------------
    md.append("## 1. B_c spectrum (MeV) with parametric errors; spread = half-range over the three main families\n")
    md.append("| State | Martin | Log | Cornell | mean ± spread | PDG |")
    md.append("|---|---|---|---|---|---|")
    table = {}
    for lab in BC_LABELS:
        vals = [1000 * res[f][lab][0] for f in main_fams]
        errs = [1000 * res[f][lab][1] for f in main_fams]
        spread = (max(vals) - min(vals)) / 2
        pdg = pdg.BC.get(lab)
        table[lab] = dict(zip(main_fams, vals)) | {"mean": float(np.mean(vals)), "spread": spread,
                                                   "param_err": dict(zip(main_fams, errs))}
        md.append(f"| {lab} | " + " | ".join(f"{v:.0f} ± {e:.0f}" for v, e in zip(vals, errs))
                  + f" | {np.mean(vals):.0f} ± {spread:.0f} | " + (f"{pdg[0]:.1f}" if pdg else "-") + " |")
    md.append("")

    # ---- 2. hyperfine test ---------------------------------------------------
    md.append("## 2. Hyperfine test\n")
    md.append("| Family | Δ(1S) | Δ(2S) | Δ(1S)−Δ(2S) | B_c*(1S) |")
    md.append("|---|---|---|---|---|")
    hyp = {}
    for f in main_fams + [p for p in power_fams if p not in main_fams]:
        d1 = 1000 * (res[f]["13S1"][0] - res[f]["11S0"][0])
        d2 = 1000 * (res[f]["23S1"][0] - res[f]["21S0"][0])
        e1 = 1000 * np.hypot(res[f]["13S1"][1], res[f]["11S0"][1])   # conservative (ignores correlation)
        hyp[f] = {"d1S": d1, "d2S": d2, "diff": d1 - d2, "bcstar": 1000 * res[f]["13S1"][0]}
        md.append(f"| {f} | {d1:.1f} | {d2:.1f} | {d1-d2:.1f} | {1000*res[f]['13S1'][0]:.0f} ± {1000*res[f]['13S1'][1]:.0f} |")
    md.append("| *running α_s (one-loop, Λ fitted to cc̄+bb̄)* | | | | |")
    for f in [k for k in fits if "_run_" in k]:
        d1 = 1000 * (res[f]["13S1"][0] - res[f]["11S0"][0])
        d2 = 1000 * (res[f]["23S1"][0] - res[f]["21S0"][0])
        hyp[f] = {"d1S": d1, "d2S": d2, "diff": d1 - d2, "bcstar": 1000 * res[f]["13S1"][0]}
        md.append(f"| {f} | {d1:.1f} | {d2:.1f} | {d1-d2:.1f} | {1000*res[f]['13S1'][0]:.0f} ± {1000*res[f]['13S1'][1]:.0f} |")
    md.append(f"| **CMS** | | | **{CMS_DIFF} ± {CMS_ERR}** | |")
    md.append(f"| **LHCb** | | | **{LHCB_DIFF:.1f}** | |")
    md.append(f"| **ATLAS 2026** | **64.5 ± 1.7** | **35.4 ± 2.4** (with CMS) | | **6339.0 ± 1.7** |")
    md.append(f"| **Lattice (HPQCD 2012 / 2010)** | **{LATTICE_D1S:.0f} ± {LATTICE_D1S_ERR:.0f}** | | | **{LATTICE_BCSTAR:.0f} ± {LATTICE_ERR:.0f}** |")
    md.append("")

    # ---- 3. exponent scan ----------------------------------------------------
    md.append("## 3. Dependence on the confining power α (V = V₀ + B(r^α − 1), pinned at r = 1 GeV⁻¹; α→0 is the log potential, α = 0.1 is Martin)\n")
    md.append("| α | χ²/dof | α_s | B_c(1S) | B_c*(1S) | B_c(2S) | 1P centroid | 1D centroid | 3S |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    scan = []
    for f in power_fams:
        a = alpha_of(f)
        d = fits[f]
        alpha_s = d["x"][d["names"].index("alpha_s")]
        r = res[f]
        p_cen = (r["13P0"][0] + 3 * r["1P1"][0] + 3 * r["1P1'"][0] + 5 * r["13P2"][0]) / 12
        d_cen = (3 * r["13D1"][0] + 5 * r["1D2"][0] + 5 * r["1D2'"][0] + 7 * r["13D3"][0]) / 20
        scan.append(dict(alpha=a, chi2=d["chi2"], dof=d["dof"], alpha_s=alpha_s,
                         **{k: 1000 * r[k][0] for k in BC_LABELS}, P1=1000 * p_cen, D1=1000 * d_cen,
                         err={k: 1000 * r[k][1] for k in BC_LABELS}))
        md.append(f"| {a:g} | {d['chi2']:.0f}/{d['dof']} | {alpha_s:.3f} | {1000*r['11S0'][0]:.0f} | {1000*r['13S1'][0]:.0f} | "
                  f"{1000*r['21S0'][0]:.0f} | {1000*p_cen:.0f} | {1000*d_cen:.0f} | {1000*r['33S1'][0]:.0f} |")
    md.append("")

    # ---- 4. radiative widths (Martin refit) ----------------------------------
    md.append("## 4. B_c radiative widths, Martin refit (keV)\n")
    P = FAMILIES["martin"]["make"](*fits["martin"]["x"])
    S = spectrum(P, "c", "b", levels=[("1S", 1, 0), ("2S", 2, 0), ("1P", 1, 1), ("1D", 1, 2), ("2P", 2, 1)],
                 r_max={"1S": 14, "2S": 14, "1P": 14, "1D": 18, "2P": 18})
    lv = full_spectrum(S)
    rad = radiative_table(S, lv, "c", "b")
    md.append("| Transition | Type | k (MeV) | Γ (keV) |")
    md.append("|---|---|---|---|")
    for r in sorted(rad, key=lambda r: -r["width_keV"]):
        if r["width_keV"] > 0.05:
            md.append(f"| {r['initial']} → {r['final']} γ | {r['type']} | {r['k_MeV']:.0f} | {r['width_keV']:.1f} |")
    md.append("")

    text = "\n".join(md)
    (OUT / "paper_results.md").write_text(text)
    (OUT / "paper_results.json").write_text(json.dumps({"spectrum": table, "hyperfine": hyp, "scan": scan, "radiative": rad}, indent=2))
    print(text)

    # ---- figures ---------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if scan:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        a = [s["alpha"] for s in scan]
        for key, name in (("11S0", "B_c(1S)"), ("13S1", "B_c*(1S)"), ("21S0", "B_c(2S)"), ("23S1", "B_c*(2S)")):
            axes[0].errorbar(a, [s[key] for s in scan], yerr=[s["err"][key] for s in scan], fmt="o-", capsize=3, label=name)
        axes[0].axhline(pdg.BC["11S0"][0], color="C0", ls="--", lw=1)
        axes[0].axhline(pdg.BC["21S0"][0], color="C2", ls="--", lw=1)
        axes[0].axhspan(LATTICE_BCSTAR - LATTICE_ERR, LATTICE_BCSTAR + LATTICE_ERR, color="C1", alpha=0.2, label="lattice B_c*")
        axes[0].set_xlabel("confining power α  (0 = log)"); axes[0].set_ylabel("mass [MeV]"); axes[0].legend(fontsize=8)
        axes[0].set_title("S-wave B_c vs potential shape (dashed: PDG)")
        for key, name in (("P1", "1P centroid"), ("D1", "1D centroid"), ("33S1", "3³S₁")):
            axes[1].plot(a, [s[key] for s in scan], "o-", label=name)
        axes[1].set_xlabel("confining power α"); axes[1].set_ylabel("mass [MeV]"); axes[1].legend(fontsize=8)
        axes[1].set_title("Higher B_c levels vs potential shape")
        fig.tight_layout(); fig.savefig(OUT / "fig_band_exponent.png", dpi=120); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    names = list(hyp)
    ax.errorbar(range(len(names)), [hyp[n]["diff"] for n in names], fmt="s", ms=7, mfc="none", label="model: Δ(1S)−Δ(2S)")
    ax.axhspan(CMS_DIFF - CMS_ERR, CMS_DIFF + CMS_ERR, color="C1", alpha=0.3, label="CMS")
    ax.axhline(LHCB_DIFF, color="C3", lw=1, label="LHCb")
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=30, fontsize=8)
    ax.set_ylabel("Δ_hyp(1S) − Δ_hyp(2S) [MeV]"); ax.set_title("Hyperfine test: predicted from cc̄ + bb̄ only")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(OUT / "fig_hyperfine_test.png", dpi=120); plt.close(fig)
    print("figures written")


if __name__ == "__main__":
    main()
