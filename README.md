# B_c spectrum from charmonium and bottomonium — paper repository

[![ORCID](https://img.shields.io/badge/ORCID-0009--0009--4678--0034-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0009-4678-0034)

**P. A. Oppong, "How well do charmonium and bottomonium determine the B_c
spectrum? A potential-model test with the B_c(2S) hyperfine difference"**
(manuscript in `paper/`, draft PDF built from `paper/main.tex`).

This repository is the publication version of the code; it continues the
history of https://github.com/apoppong24/bc-quarkonium (the thesis refactor)
and adds the self-consistent hyperfine test, the uncertainty budget, ⟨v²⟩,
decay constants and the experimental-prediction tables. The original MSc
thesis and notebooks are archived at
https://github.com/apoppong24/bc-meson-msc-thesis.

Reproduce everything (≈1 h on one core, or ≈2 min using the shipped fits):

```bash
pip install -r requirements.txt
python bc_paper_results.py            # spectrum, hyperfine test, α-scan, E1/M1
python uncertainty_budget.py          # systematic budget + self-consistent test
python experimental_predictions.py    # Table 10 and cascade signatures
python paper/make_paper_assets.py     # LaTeX tables + figures
cd paper && latexmk -pdf main.tex
```

---

# quarkonium — Numerov solver for heavy-quark mesons and the B_c spectrum

Refactored code for the MSc project **"B_c Meson Masses"** (Peter Agyemang
Oppong, University of Glasgow): the eight original Jupyter notebooks are
replaced by one small, tested Python package plus a script that regenerates
every table in the thesis. The thesis itself and the untouched original
notebooks are archived at
https://github.com/apoppong24/bc-meson-msc-thesis.

## What changed

| Before (notebooks) | After (this package) |
|---|---|
| Numerov loop copy-pasted ~20× with different masses/potentials | One `numerov()` + one `solve()` |
| Eigenvalues found by hand from 3 trial energies and a plot | Automatic node-counting bisection (the method described in thesis §3.3), converged to 1e-10 GeV |
| S states only; P-state code absent | Any `(n, l)` — the centrifugal term is built in |
| `input()` calls, can't run unattended | Pure functions; `python reproduce_thesis.py` regenerates Tables 2–15 in ~30 s |
| Results copied into tables by hand | `results/tables.md` + `results/results.json` written automatically |
| No tests | 31 regression tests against the thesis numbers and physics invariants |
| S-state hyperfine only | **Fine structure for any L**: spin–orbit, tensor and hyperfine shifts (Eichten & Quigg 1994 formalism) with ³L_L–¹L_L mixing for B_c |

## Install & run

```bash
pip install numpy matplotlib pytest      # numpy is the only hard dependency
python reproduce_thesis.py --plots       # tables -> results/, Figs. 3, 8, 9, 10-12 -> results/*.png
python fine_structure_tables.py --potential martin   # J-resolved cc, bb, Bc spectra through 4S/3P/2D
python refit_pdg.py martin      # refit one family to PDG 2026 (~4 min; also: log, cornell)
python refit_pdg.py --report    # combine fits, residual table, Bc predictions with model spread
python scan_exponent.py 0.35        # refit V = V0 + B(r^a - 1) at fixed confining power a (~4 min each)
python refit_pdg.py martin_run_geometric   # running-alpha_s variants: {martin,log}_run_{reduced,geometric}
python bc_paper_results.py          # paper tables: spectrum + errors, hyperfine test, a-scan, E1/M1 widths
python refit_pdg.py martin --floor 5      # budget variants: --floor 5|20, --exclude D|3P; families *_sc = smeared hyperfine
python uncertainty_budget.py        # systematic budget; self-consistent-hyperfine comparison
python experimental_predictions.py  # per-state masses, photon lines, cascades, decay constants, falsifiability
pytest                                   # ~25 s
```

## Usage

```python
from quarkonium import MARTIN, LOGARITHMIC, CORNELL, solve, spectrum
from quarkonium import hyperfine_splitting, singlet_triplet

# one state: bottomonium 1S in the Martin (power-law) potential
s = solve(MARTIN, MARTIN.m_b, MARTIN.m_b, n=1, l=0)
s.E, s.mass, s.R_at_origin_sq()     # -0.8742, 9.4738, 4.423

# the five thesis levels for B_c
bc = spectrum(MARTIN, "c", "b")      # {'1S': RadialState, '2S': ..., '1P': ..., '2P': ...}
bc["1P"].E                           # -0.21518

# pseudoscalar / vector B_c masses (thesis Table 15)
singlet_triplet(bc["1S"].mass, hyperfine_splitting(bc["1S"]))   # (6.2768, 6.3474)

# fine structure: every spin state built on a radial level (any L)
from quarkonium.spin_dependent import fine_structure
for f in fine_structure(bc["1P"]):
    print(f.label, f.J, round(f.mass, 4))     # 13P0, 1P1, 1P1', 13P2 (J=1 states mixed)

# wavefunction for plotting
import matplotlib.pyplot as plt
plt.plot(s.r, s.R)
```

Define your own potential by constructing a `Potential(name, V, m_c, m_b, alpha_s)`.

## Layout

```
quarkonium/
  potentials.py   Cornell, Logarithmic, Martin + Eichten–Quigg parameters   (thesis §2.5)
  numerov.py      the three-term recursion, eq. (35)                        (thesis §3.2)
  solver.py       V_eff, boundary conditions, bisection, normalisation,
                  |R(0)|²                                                   (thesis §3.1, 3.3)
  observables.py  hyperfine splitting (41), spin averages (38,39),
                  singlet/triplet (43), B_c formula (42), level splittings  (thesis §4)
  spectrum.py     solve the 1S 2S 3S 1P 2P set for a flavour pair
  plots.py        Figs. 3 (potentials), 8 (shooting), 9 (wavefunctions), 10-12 (splittings)
  spin_dependent.py  first-order spin-dependent shifts, E&Q 1994 eqs. (2.15)-(2.25)
  pdg.py      PDG 2026 cc̄ / bb̄ / B_c masses used in the refit
  fit.py          joint least-squares refit of (potential params, m_c, m_b, α_s); fixed-exponent families
  uncertainty.py  propagate fit covariance to predicted masses
  radiative.py    E1/M1 widths (E&Q 1994 §III), Wigner 6-j statistical factors
reproduce_thesis.py   regenerate Tables 2–15 and Figs. 3, 8, 9, 10–12
fine_structure_tables.py  J-resolved spectra (1S-4S, 1P-3P, 1D-2D) vs PDG and E&Q Table C
refit_pdg.py      refit to PDG 2026 and predict B_c; writes results/refit_pdg.md
scan_exponent.py      refit at fixed confining power α
bc_paper_results.py   assemble paper-level results -> results/paper_results.md, fig_band_exponent.png, fig_hyperfine_test.png
tests/                regression tests vs. thesis values
```

## Notes found while refactoring

* **Table 4, Cornell cc̄ 1S** is printed as 0.2657 GeV; the notebook and the
  mass in Table 5 (3.9375 = 2×1.84 + 0.2575) both correspond to **0.2575**.
* **Eq. (42)** with the *model* masses `M_bb = 2m_b + E_bb` reduces
  algebraically to `m_b + m_c + E_Bc` (eq. 40), so it does not cancel the
  quark-mass dependence as stated. It only does so if the *experimental*
  spin-averaged masses are inserted; `bc_mass()` supports both.
* The hand-tuned eigenvalues are typically within 0.1–0.5 MeV of the converged
  ones; the largest gaps are log bb̄ 2S/3S (≈2 MeV) and log cb̄ 3S (≈6 MeV),
  which were searched on a short `r_max`.
* With a log potential the level spacing is exactly independent of reduced
  mass, so the Logarithmic columns of Tables 7, 8 and 15 are identical and
  Table 13's log column is a constant (0.06088 GeV). The converged results
  show this exactly; the thesis values show it approximately.
* The error bars in the regenerated Figs. 10–12 are *defined here* as each
  model's deviation from the mean of the three potential models for that
  splitting (a model-dependence estimate). The thesis figures' bars are
  described only as "the deviation that occurred in the potential models",
  so this is a reconstruction, not a reproduction.
* The `alpha_s` used for the hyperfine splitting differs per potential in
  the notebooks (Cornell 0.31, log 0.37, Martin 0.43); these are kept as
  attributes of each `Potential` and can be overridden in
  `hyperfine_splitting(state, alpha_s=...)`.

## Extension beyond the thesis: fine structure

`quarkonium/spin_dependent.py` adds the first-order spin-dependent shifts of
Eichten & Quigg (1994), eqs. (2.15)–(2.25):

    Δ = Σ T_k,   T̃1 = −⟨V'/r⟩ + 2T̃2,  T̃2 = (4α_s/3)⟨r⁻³⟩,  T̃3 = (8α_s/9)|R(0)|²,  T̃4 = (α_s/3)⟨r⁻³⟩

evaluated on the solved wavefunctions. For unequal masses (B_c) the ³L_L and
¹L_L states mix; the 2×2 matrix is built in the LS basis and diagonalised
(trace and determinant agree with E&Q's jj matrices, eq. 2.21).

Checks built into the tests: Landé centroid rule, ¹P₁ at the centroid for
equal masses, S-wave reduction to the thesis hyperfine formula, trace of the
B_c mixing matrix, and B_c 1P levels within 20–40 MeV of E&Q Table C.

## Refit to PDG 2026

`refit_pdg.py` refits each family (potential parameters, m_c, m_b, α_s)
jointly to 27 J-resolved charmonium and bottomonium levels from PDG 2026
(`quarkonium/pdg.py`), with a 10 MeV theory floor added in quadrature.
Results are in `results/refit_pdg.md` / `fit_<family>.json`, and
`quarkonium.fit.load_fit("martin")` returns the refitted `Potential`.

Headline: with refitted parameters the three families predict
B_c(1S) = 6276 ± 14 MeV (PDG 6274.5), B_c*(1S) = 6356 ± 17 MeV,
B_c(2S) = 6884 ± 6 MeV (PDG 6871.2), where ± is the model spread.

## Paper-level results (`results/paper_results.md`)

* B_c spectrum through 3S/2P/1D from cc̄+bb̄ data only, with parametric
  errors (covariance propagation) and model spread.
* Hyperfine test: Δ_hyp(1S) − Δ_hyp(2S) predicted 28–40 MeV by the three
  families vs CMS 29.1 ± 1.7 / LHCb 30.9 MeV. With a *fixed* α_s the absolute
  B_c* hyperfine (76–82 MeV) overshoots lattice (~55 MeV); with a one-loop
  *running* α_s(μ) (Λ fitted to cc̄+bb̄, `*_run_*` families; μ = 2μ_red or
  √(m_c m_b)) the Martin family gives Δ(1S) = 63–67 MeV vs **ATLAS 2026:
  64.5 ± 1.7 MeV** (arXiv:2605.16228), Δ(2S) = 35–38 vs 35.4 ± 2.4 (ATLAS+CMS),
  and the fit χ² improves (85→70, 74→59).
* Confining-power scan: χ² prefers α ≲ 0.1 (near-log); the hyperfine
  difference falls from 37 MeV (α=0.05) to 0 (α=1) and independently prefers
  α ≈ 0.1–0.2.
* E1/M1 widths for B_c (1P→1S 60–110 keV, 1D→1P 20–95 keV).
  Mixed J=L states are treated as their dominant S=1 component in E1 rates;
  a mixing-angle-weighted version is a straightforward refinement.

Known limitations: first-order perturbation theory with a fixed α_s leaves
10–40 MeV residuals (worst for cc̄ D-waves, which sit at the DD̄ threshold);
no relativistic corrections; no coupled-channel effects.
