#!/usr/bin/env python3
"""Systematic uncertainty budget for the B_c predictions.

Sources (each evaluated as a shift of the prediction relative to the
baseline fit of the same family):
  floor     theory floor 5 and 20 MeV instead of 10          (half the range)
  data      drop the charmonium D waves / drop the bb̄ 3P    (max |shift|)
  scale     α_s scale √(m_c m_b) vs 2μ for running families (|difference|)
  model     half-range over the three potential families (fixed-α_s baseline)
  fit       propagated parameter covariance (baseline)

Writes results/uncertainty_budget.md / .json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from quarkonium.fit import FAMILIES, predict, register_power_family
from quarkonium.uncertainty import predict_with_errors

OUT = Path(__file__).parent / "results"
LAB = ["11S0", "13S1", "21S0", "23S1", "13P0", "1P1", "1P1'", "13P2", "13D1", "13D3", "33S1"]


def load(name):
    p = OUT / f"fit_{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def pred(name):
    d = load(name)
    if d is None:
        return None
    fam = d["family"]
    if fam.startswith("power_"):
        register_power_family(d["alpha"])
    P = FAMILIES[fam]["make"](*d["x"])
    return {l: 1000 * v for l, v in predict(P, "cb", LAB, h=1e-3).items()}


def main():
    rows = {}
    base = {f: pred(f) for f in ("martin", "log", "cornell", "martin_run_geometric", "log_run_geometric",
                                 "martin_run_geometric_sc", "martin_sc", "log_sc", "log_run_geometric_sc")}
    base = {k: v for k, v in base.items() if v is not None}
    hyp = lambda p: {"d1S": p["13S1"] - p["11S0"], "d2S": p["23S1"] - p["21S0"]}

    def shift_table(family, variants, kind):
        """kind='halfrange' -> half the range over base+variants; 'max' -> max |variant-base|."""
        b = base[family]
        vs = [pred(v) for v in variants]
        vs = [v for v in vs if v is not None]
        out = {}
        for l in LAB:
            vals = [b[l]] + [v[l] for v in vs]
            out[l] = (max(vals) - min(vals)) / 2 if kind == "halfrange" else max(abs(v[l] - b[l]) for v in vs)
        for k in ("d1S", "d2S"):
            vals = [hyp(b)[k]] + [hyp(v)[k] for v in vs]
            out[k] = (max(vals) - min(vals)) / 2 if kind == "halfrange" else max(abs(hyp(v)[k] - hyp(b)[k]) for v in vs)
        return out

    budget = {}
    for fam in ("martin", "martin_run_geometric"):
        bud = {}
        bud["floor"] = shift_table(fam, [f"{fam}__floor5", f"{fam}__floor20"], "halfrange")
        bud["data"] = shift_table(fam, [f"{fam}__noD", f"{fam}__no3P"], "max")
        if "run" in fam:
            bud["scale"] = shift_table(fam, ["martin_run_reduced"], "max")
        # NB: the self-consistent smeared hyperfine variants are NOT a budget
        # entry: with a mass-independent smearing they fail to describe the
        # bb̄ hyperfine splitting (see the table below) and are reported as an
        # excluded alternative rather than as a systematic.
        # model spread over the three fixed-alpha families (or the two running ones)
        fams = ("martin", "log", "cornell") if "run" not in fam else ("martin_run_geometric", "log_run_geometric")
        mod = {}
        for l in LAB:
            vals = [base[f][l] for f in fams]
            mod[l] = (max(vals) - min(vals)) / 2
        for k in ("d1S", "d2S"):
            vals = [hyp(base[f])[k] for f in fams]
            mod[k] = (max(vals) - min(vals)) / 2
        bud["model"] = mod
        d = load(fam)
        pe = predict_with_errors(fam, d["x"], d["cov"], "cb", LAB, h=1e-3)
        fit_err = {l: 1000 * pe[l][1] for l in LAB}
        fit_err["d1S"] = 1000 * np.hypot(pe["13S1"][1], pe["11S0"][1])
        fit_err["d2S"] = 1000 * np.hypot(pe["23S1"][1], pe["21S0"][1])
        bud["fit"] = fit_err
        budget[fam] = {k: v for k, v in bud.items() if v is not None}

    keys = LAB + ["d1S", "d2S"]
    md = ["# Uncertainty budget for B_c predictions (MeV)\n"]
    for fam, bud in budget.items():
        md.append(f"## {fam}\n")
        srcs = list(bud)
        md.append("| State | " + " | ".join(srcs) + " | total (quad.) | central |")
        md.append("|---|" + "---|" * (len(srcs) + 2))
        for l in keys:
            vals = [bud[s][l] for s in srcs]
            tot = float(np.sqrt(sum(v * v for v in vals)))
            cen = base[fam][l] if l in LAB else hyp(base[fam])[l]
            md.append(f"| {l} | " + " | ".join(f"{v:.1f}" for v in vals) + f" | **{tot:.1f}** | {cen:.1f} |")
        md.append("")
        # keep in budget dict
        for l in keys:
            bud.setdefault("_total", {})[l] = float(np.sqrt(sum(bud[s][l] ** 2 for s in srcs)))
            bud.setdefault("_central", {})[l] = base[fam][l] if l in LAB else hyp(base[fam])[l]
    # self-consistent hyperfine results for the record
    md.append("## Self-consistent (smeared) hyperfine fits — excluded alternative\n")
    md.append("Gaussian-smeared contact term (width σ, mass-independent) solved self-consistently, Λ or α_s refitted. "
              "These fits reproduce the cc̄ 1S splitting but not the bb̄ one, and pull m_c to unphysical values; "
              "they are therefore not used as a systematic.\n")
    md.append("| Family | χ²/dof | Λ or α_s | m_c | cc̄ Δ(1S) fit | bb̄ Δ(1S) fit | B_c Δ(1S) | Δ(2S) | Δ(1S)−Δ(2S) | B_c*(1S) |")
    md.append("|---|---|---|---|---|---|---|---|---|---|")
    sc = {}
    for name in ("martin_sc", "log_sc", "martin_run_geometric_sc", "log_run_geometric_sc",
                 "martin_run_geometric_sc_s0.8", "martin_run_geometric_sc_s1.2"):
        d = load(name)
        if d is None:
            continue
        p = pred(name); h = hyp(p)
        rr = {(x["pair"], x["label"]): x["model"] for x in d["residuals"]}
        cc1 = 1000 * (rr[("cc", "13S1")] - rr[("cc", "11S0")]); bb1 = 1000 * (rr[("bb", "13S1")] - rr[("bb", "11S0")])
        mc = d["x"][d["names"].index("m_c")]
        sc[name] = {"chi2": d["chi2"], "dof": d["dof"], "coupling": d["x"][-1], "m_c": mc, "cc_d1S": cc1, "bb_d1S": bb1, **h, "bcstar": p["13S1"]}
        md.append(f"| {name} | {d['chi2']:.1f}/{d['dof']} | {d['x'][-1]:.3f} | {mc:.2f} | {cc1:.0f} | {bb1:.0f} | {h['d1S']:.1f} | {h['d2S']:.1f} | {h['d1S']-h['d2S']:.1f} | {p['13S1']:.0f} |")
    md.append("| experiment | | | | 112.8 | 61.7 | 64.5 ± 1.8 | 35.4 ± 2.5 | 29.1 ± 1.7 | 6339.0 ± 1.8 |")
    # perturbative fits' cc/bb hyperfine for comparison
    md.append("\nFor reference, the perturbative fits give (cc̄ Δ(1S), bb̄ Δ(1S)): " + "; ".join(
        f"{f}: ({1000*(rr2[('cc','13S1')]-rr2[('cc','11S0')]):.0f}, {1000*(rr2[('bb','13S1')]-rr2[('bb','11S0')]):.0f})"
        for f in ("martin", "log", "martin_run_geometric", "log_run_geometric")
        for rr2 in [{(x["pair"], x["label"]): x["model"] for x in load(f)["residuals"]}]) + " MeV.")
    text = "\n".join(md)
    (OUT / "uncertainty_budget.md").write_text(text)
    (OUT / "uncertainty_budget.json").write_text(json.dumps({"budget": budget, "selfconsistent": sc}, indent=2))
    print(text)


if __name__ == "__main__":
    main()
