#!/usr/bin/env python3
"""Refit one potential family to PDG 2026 and save results/fit_<family>.json.

    python refit_pdg.py martin      # also: log, cornell
    python refit_pdg.py --report    # combine the three fits, predict B_c
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).parent / "results"


def run_fit(family, sigma_th=10.0, exclude=(), tag=""):
    from quarkonium.fit import fit
    r = fit(family, sigma_th=sigma_th, verbose=True, exclude=exclude)
    print(r.summary())
    OUT.mkdir(exist_ok=True)
    (OUT / f"fit_{family}{tag}.json").write_text(json.dumps({
        "family": family, "tag": tag, "sigma_th": sigma_th, "exclude": list(exclude),
        "names": r.names, "x": r.x.tolist(), "errors": r.errors.tolist(),
        "cov": None if r.cov is None else r.cov.tolist(),
        "chi2": r.chi2, "dof": r.dof, "seconds": r.seconds,
        "residuals": [{"pair": p, "label": l, "exp": e, "model": m, "pull": q} for p, l, e, m, q in r.residuals],
    }, indent=2))


def report():
    from quarkonium import pdg
    from quarkonium.fit import FAMILIES, predict
    from quarkonium.potentials import CORNELL, LOGARITHMIC, MARTIN

    fits = {f: json.loads((OUT / f"fit_{f}.json").read_text()) for f in ("martin", "log", "cornell")}
    old = {"martin": MARTIN, "log": LOGARITHMIC, "cornell": CORNELL}
    md = ["# Refit of potential models to PDG 2026\n",
          "Joint least-squares fit of each potential family (parameters, m_c, m_b, α_s) to "
          f"{len(fits['martin']['residuals'])} measured charmonium and bottomonium levels "
          "(J-resolved, below open-flavour threshold), σ² = σ_exp² + (10 MeV)².\n"]

    md.append("## Fitted parameters (Eichten & Quigg 1994 values in brackets)\n")
    md.append("| Family | χ²/dof | " + " | ".join(f"p{i+1}" for i in range(5)) + " |")
    md.append("|---|---|" + "---|" * 5)
    for f, d in fits.items():
        P0 = old[f]
        old_vals = {"martin": [6.898, 0.1, P0.m_c, P0.m_b, P0.alpha_s],
                    "log": [0.733, P0.m_c, P0.m_b, P0.alpha_s],
                    "cornell": [0.52, 2.34, P0.m_c, P0.m_b, P0.alpha_s]}[f]
        cells = [f"{n} = {x:.4f} ± {e:.4f} [{o:.3g}]" for n, x, e, o in zip(d["names"], d["x"], d["errors"], old_vals)]
        cells += [""] * (5 - len(cells))
        md.append(f"| {f} | {d['chi2']:.1f}/{d['dof']} | " + " | ".join(cells) + " |")
    md.append("")
    md.append("Parameter errors are from the covariance of the least-squares fit at the 10 MeV "
              "theory floor. For the Martin family B, α, m_c and m_b are strongly correlated "
              "(B·r^α ≈ B(1 + α ln r) for small α, so B carries a constant that is degenerate with "
              "the quark masses); the individual errors are large but the predicted *masses* are "
              "not. The combination m_c + m_b that enters B_c is well constrained by the cc̄ and bb̄ fits.\n")

    md.append("## Fit residuals (MeV)\n")
    md.append("| State | PDG 2026 | Martin | Log | Cornell |")
    md.append("|---|---|---|---|---|")
    labels = [(r["pair"], r["label"], r["exp"]) for r in fits["martin"]["residuals"]]
    res = {f: {(r["pair"], r["label"]): r["model"] for r in d["residuals"]} for f, d in fits.items()}
    for pair, lab, exp in labels:
        md.append(f"| {pair} {lab} | {1000*exp:.1f} | " + " | ".join(
            f"{1000*(res[f][(pair, lab)]-exp):+.1f}" for f in fits) + " |")
    md.append("")

    # B_c predictions with refitted parameters
    bc_labels = ["11S0", "13S1", "13P0", "1P1", "1P1'", "13P2", "21S0", "23S1",
                 "13D1", "1D2", "1D2'", "13D3", "23P0", "2P1", "2P1'", "23P2", "31S0", "33S1"]
    preds = {}
    for f, d in fits.items():
        P = FAMILIES[f]["make"](*d["x"])
        preds[f] = predict(P, "cb", bc_labels, h=1e-3)
    md.append("## B_c spectrum with refitted potentials (MeV)\n")
    md.append("Spread = half the range across the three families (model dependence).\n")
    md.append("| State | Martin | Log | Cornell | Mean ± spread | PDG 2026 |")
    md.append("|---|---|---|---|---|---|")
    bc_out = {}
    for lab in bc_labels:
        v = np.array([preds[f][lab] for f in fits]) * 1000
        pdg = pdg.BC.get(lab)
        bc_out[lab] = {"martin": v[0], "log": v[1], "cornell": v[2], "mean": v.mean(), "spread": (v.max() - v.min()) / 2}
        md.append(f"| {lab} | {v[0]:.1f} | {v[1]:.1f} | {v[2]:.1f} | {v.mean():.1f} ± {(v.max()-v.min())/2:.1f} | "
                  + (f"{pdg[0]:.1f} ± {pdg[1]}" if pdg else "-") + " |")
    text = "\n".join(md)
    (OUT / "refit_pdg.md").write_text(text)
    (OUT / "bc_predictions_refit.json").write_text(json.dumps(bc_out, indent=2))
    print(text)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("family", nargs="?")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--floor", type=float, default=10.0, help="theory floor in MeV")
    ap.add_argument("--exclude", nargs="*", default=[], help="drop levels whose label contains these (e.g. D 3P)")
    a = ap.parse_args()
    if a.report:
        report()
    else:
        tag = ""
        if a.floor != 10.0:
            tag += f"__floor{a.floor:g}"
        if a.exclude:
            tag += "__no" + "".join(a.exclude)
        run_fit(a.family, sigma_th=a.floor, exclude=tuple(a.exclude), tag=tag)
