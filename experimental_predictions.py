#!/usr/bin/env python3
"""Predictions aimed at experiment, from the Martin running-α_s fit
(the family that reproduces the ATLAS/CMS/LHCb hyperfine data):

  1. For every excited B_c level through 2P/1D: mass ± (fit error, model
     spread), all E1/M1 photon lines with photon energy, partial width and
     branching fraction (radiative only), and the dominant decay chain.
  2. Cascade signatures for the 1P and 1D states: total photon energy and
     the final S-wave state reached.
  3. <v²> for every state (how relativistic each prediction is).
  4. Decay constants f_ηc, f_ηb, f_Bc vs lattice.
  5. Falsifiability thresholds: predictions that would exclude the whole
     cc̄+bb̄-determined potential family.

Writes results/experimental_predictions.md and .json.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from quarkonium import pdg2024
from quarkonium.fit import FAMILIES, predict
from quarkonium.radiative import radiative_table
from quarkonium.solver import solve
from quarkonium.spin_dependent import fine_structure, full_spectrum
from quarkonium.uncertainty import predict_with_errors

OUT = Path(__file__).parent / "results"
LEVELS = [("1S", 1, 0), ("2S", 2, 0), ("3S", 3, 0), ("1P", 1, 1), ("2P", 2, 1), ("1D", 1, 2)]
R_MAX = {"1S": 14, "2S": 14, "3S": 18, "1P": 14, "2P": 18, "1D": 18}
LAB = ["11S0", "13S1", "13P0", "1P1", "1P1'", "13P2", "21S0", "23S1",
       "13D1", "1D2", "1D2'", "13D3", "23P0", "2P1", "2P1'", "23P2", "31S0", "33S1"]
LATTICE_F = {"cc": ("f_ηc", *pdg2024.LATTICE_F["cc"]), "bb": ("f_ηb", *pdg2024.LATTICE_F["bb"]), "cb": ("f_Bc", *pdg2024.LATTICE_F["cb"])}   # HPQCD, MeV


def tex(lab):
    return lab


def main():
    fits = {f: json.loads((OUT / f"fit_{f}.json").read_text()) for f in ("martin_run_geometric", "log_run_geometric", "martin", "log", "cornell")}
    main_fam = "martin_run_geometric"
    P = FAMILIES[main_fam]["make"](*fits[main_fam]["x"])

    # ---- masses with parametric errors (main family) and model spread (all five families)
    central = predict_with_errors(main_fam, fits[main_fam]["x"], fits[main_fam]["cov"], "cb", LAB, h=1e-3)
    others = {f: predict(FAMILIES[f]["make"](*fits[f]["x"]), "cb", LAB, h=1e-3) for f in fits}
    spread = {l: (max(others[f][l] for f in fits) - min(others[f][l] for f in fits)) / 2 for l in LAB}

    # ---- radiative transitions from the solved states of the main family
    S = {lab: solve(P, P.m_c, P.m_b, n, L, r_max=R_MAX[lab], h=1e-3) for lab, n, L in LEVELS}
    levels = full_spectrum(S)
    rad = radiative_table(S, levels, "c", "b")
    by_init = {}
    for r in rad:
        by_init.setdefault(r["initial"], []).append(r)
    for k in by_init:
        tot = sum(r["width_keV"] for r in by_init[k])
        for r in by_init[k]:
            r["BR_rad"] = r["width_keV"] / tot if tot > 0 else 0.0
        by_init[k].sort(key=lambda r: -r["width_keV"])

    md = ["# Predictions for experiment — B_c excited states\n",
          f"Central values: Martin potential with one-loop running α_s refitted to PDG 2026 cc̄+bb̄ "
          f"(Λ = {fits[main_fam]['x'][-1]:.3f} GeV). Errors: (fit) propagated parameter covariance; "
          f"(model) half-range over the five fitted families. Widths are radiative only; hadronic "
          f"transitions (ππ) and the B_c(2S)/B_c*(2S) → B_c(*)ππ channels are not included.\n"]

    # ---- 1. per-state table
    md.append("## 1. Masses and photon lines\n")
    md.append("| State | J^P | Mass [MeV] | ± fit | ± model | ⟨v_c²⟩ | ⟨v_b²⟩ | Photon lines: final state, E_γ [MeV], Γ_E1 [keV], BR_rad |")
    md.append("|---|---|---|---|---|---|---|---|")
    JP = {"S0": "0⁻", "S1": "1⁻", "P0": "0⁺", "P1": "1⁺", "P2": "2⁺", "D1": "1⁻", "D2": "2⁻", "D3": "3⁻"}
    state_rows = {}
    for lab in LAB:
        core = lab.rstrip("'")
        jp = JP[core[-2:]] if len(core) == 4 else JP[core[-2:]]
        n = int(core[0]); L = "SPD".index(core[-2])
        st = S[f"{n}{'SPD'[L]}"]
        v = st.v2()
        lines = by_init.get(lab, [])
        line_txt = "; ".join(f"{r['final']} {r['k_MeV']:.0f} {r['width_keV']:.0f} {100*r['BR_rad']:.0f}%" for r in lines[:4]) or "—"
        m = 1000 * central[lab][0]
        state_rows[lab] = dict(mass=m, fit_err=1000 * central[lab][1], model_err=1000 * spread[lab],
                               v2_c=v[0], v2_b=v[1], lines=lines[:6])
        md.append(f"| {lab} | {jp} | {m:.0f} | {1000*central[lab][1]:.0f} | {1000*spread[lab]:.0f} | {v[0]:.2f} | {v[1]:.3f} | {line_txt} |")
    md.append("")

    # ---- 2. cascades
    md.append("## 2. Cascade signatures (dominant radiative chains)\n")
    md.append("| Initial | Chain | Photon energies [MeV] | Σ E_γ [MeV] | Final |")
    md.append("|---|---|---|---|---|")
    casc = []
    for lab in ["13P2", "1P1'", "1P1", "13P0", "13D1", "1D2", "1D2'", "13D3"]:
        chain, es, cur = [lab], [], lab
        for _ in range(3):
            nxt = by_init.get(cur, [])
            if not nxt:
                break
            best = nxt[0]
            es.append(best["k_MeV"]); cur = best["final"]; chain.append(cur)
        casc.append(dict(initial=lab, chain=chain, E=es))
        md.append(f"| {lab} | {' → '.join(chain)} | {', '.join(f'{e:.0f}' for e in es)} | {sum(es):.0f} | {chain[-1]} |")
    md.append("")

    # ---- 3. decay constants
    md.append("## 3. Decay constants (pseudoscalar, van Royen–Weisskopf with one-loop correction)\n")
    md.append("| Meson | f [MeV] | lattice (HPQCD) |")
    md.append("|---|---|---|")
    fdec = {}
    for pair, (q1, q2) in (("cc", ("c", "c")), ("bb", ("b", "b")), ("cb", ("c", "b"))):
        st = solve(P, P.mass(q1), P.mass(q2), 1, 0, r_max=12, h=1e-3)
        f = 1000 * st.decay_constant()
        name, lat, laterr = LATTICE_F[pair]
        fdec[pair] = f
        md.append(f"| {name} | {f:.0f} | {lat:.1f} ± {laterr:.1f} |")
    md.append("")

    # ---- 4. hyperfine summary for the record
    d1 = 1000 * (central["13S1"][0] - central["11S0"][0]); d2 = 1000 * (central["23S1"][0] - central["21S0"][0])
    md.append("## 4. Hyperfine sector\n")
    md.append(f"Δ(1S) = {d1:.1f} MeV (ATLAS {pdg2024.ATLAS_HYP_1S[0]} ± {pdg2024.ATLAS_HYP_1S[1]}); "
              f"Δ(2S) = {d2:.1f} MeV (ATLAS+CMS 35.4 ± 2.4); Δ(1S)−Δ(2S) = {d1-d2:.1f} MeV "
              f"(CMS {pdg2024.CMS_HYP_DIFF[0]} ± {pdg2024.CMS_HYP_DIFF[1]}, LHCb {pdg2024.LHCB_HYP_DIFF}).\n")

    # ---- 5. falsifiability
    p1 = [1000 * others[f]["13P0"] for f in fits]; p2 = [1000 * others[f]["13P2"] for f in fits]
    cen = [(1000 * (others[f]["13P0"] + 3 * others[f]["1P1"] + 3 * others[f]["1P1'"] + 5 * others[f]["13P2"]) / 12) for f in fits]
    s3 = [1000 * others[f]["33S1"] for f in fits]
    md.append("## 5. What would falsify the cc̄+bb̄-determined description\n")
    md.append(f"* A 1P centroid outside {min(cen):.0f}–{max(cen):.0f} MeV, or a ³P₂–³P₀ splitting outside "
              f"{min(b-a for a,b in zip(p1,p2)):.0f}–{max(b-a for a,b in zip(p1,p2)):.0f} MeV: no potential in the scanned family, "
              f"fitted to cc̄+bb̄, produces it.")
    md.append(f"* A B_c*(2S)–B_c(2S) splitting below {d2-5:.0f} MeV or above {d2+5:.0f} MeV (currently 35.4 ± 2.4 from ATLAS+CMS).")
    md.append(f"* A 3³S₁ below {min(s3):.0f} MeV or above {max(s3):.0f} MeV.")
    md.append("* A 1D multiplet whose centroid lies more than 30 MeV from 7035 MeV, or with the ³D₁ above the ³D₃.")
    md.append("* Any B_c excited state observed below 6690 MeV.")
    md.append("")

    text = "\n".join(md)
    (OUT / "experimental_predictions.md").write_text(text)
    (OUT / "experimental_predictions.json").write_text(json.dumps(
        {"states": state_rows, "cascades": casc, "decay_constants": fdec, "hyperfine": {"d1S": d1, "d2S": d2}}, indent=2, default=float))
    print(text)


if __name__ == "__main__":
    main()
