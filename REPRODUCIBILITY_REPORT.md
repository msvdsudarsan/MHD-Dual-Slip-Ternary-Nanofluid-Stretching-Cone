# Reproducibility report: `ternary_cone_master.py`

This note documents an independent, fresh re-run of the code deposited in
this repository, confirming that it reproduces every numerical value
reported in the manuscript and summarizing the checks used to establish
solution accuracy.

## Environment

- Python 3.12, SciPy 1.17, NumPy, Matplotlib
- Single core, 2.1 GHz Intel Xeon-class processor
- Full run time (all tables and figures): 9.5 s

## What was re-run

```
python ternary_cone_master.py
```

was executed end-to-end from a clean checkout, and every printed
numerical value was compared against the corresponding number quoted in
the manuscript text and tables.

## Result

| Table | Status |
|---|---|
| Effective property ratios (A1-A5, A7) | Match |
| Baseline solve (f_inf, Qc, Cf, Nu, Ns(0), Be(0), decay rates) | Match |
| Admissibility map (Table 2, S x Q) | Match |
| Domain-truncation study (Table 3, binary & ternary, eta_inf = 5..30) | Match |
| Mesh-independence study (Table 4, N = 40..1280) | Match |
| Clear-fluid validation (Table 5) | Match |
| TiO2 incremental sweep (Table 6) | Match |
| sigma_TiO2 sensitivity (Table 7) | Match |
| Hamilton-Crosser shape-factor sweep (Table 8) | Match |
| Full single-parameter sweep (Tables 9-10) | Match, all 58 rows |

Baseline state (Eq. (39) of the manuscript):

```
f_inf   = 1.300178
Qc      = 0.190374
Cf      = 0.232017
Nu      = 1.463872
Ns(0)   = 0.668808
Be(0)   = 0.947692
```

## How solution accuracy is established

This code is the single computational record behind the manuscript, so
accuracy is not established by comparison against a second solver but by
three internal checks, each reproduced by this script:

1. **Mesh/tolerance refinement.** The collocation mesh is refined from
   N = 40 to N = 1280 nodes with the solver tolerance tightened from
   1e-6 to 1e-10 (`mesh_table`). The wall skin-friction coefficient is
   converged to 8 decimal places by N = 40; the Nusselt number is
   converged to 8 decimal places by N = 160 and does not change further.
2. **Clear-fluid limiting case.** Setting all volume fractions, M, Rd,
   Q, L1, L2, and S to zero reduces the model to the classical
   clear-fluid cone-flow problem (`validation_table`); the wall
   gradients vary monotonically with the stretching ratio as expected.
3. **Analytical far-field comparison.** The domain-truncation study
   (`trunc_table`, eta_inf = 5..30) is compared against the closed-form
   decay rate of Proposition 1 in the manuscript; the empirical and
   analytical decay rates agree to within the discretization error
   expected from a finite set of eta_inf values.

## How to reproduce

```
pip install -r requirements.txt
python ternary_cone_master.py
```

This prints the full console report (property ratios, baseline solve,
and all ten tables) and writes `velocity_profile.pdf`,
`temperature_profile.pdf`, `entropy_Ns.pdf`, `bejan_number.pdf`, and
`convergence.pdf` to the working directory.

## V17 update: fresh re-run and mixing-order permutation check

For the V17 revision, `ternary_cone_master.py` was re-run end-to-end
again from a clean checkout, on a Python 3.12-class environment with
NumPy 2.4.4, SciPy 1.17.1, and Matplotlib 3.10.8 (see `requirements.txt`
for the exact pinned versions). All values reproduced exactly:

```
f_inf   = 1.300178
Qc      = 0.190374
Cf      = 0.232017
Nu      = 1.463872
Ns(0)   = 0.668808
Be(0)   = 0.947692
```

In addition, `mixing_order_check.py` was added and run to directly
re-solve the coupled boundary-value problem under all six permutations
of the ternary mixing order (previously the manuscript only estimated
this effect from the separate shape-factor sensitivity study). Results
are saved to `data/mixing_order_permutations.csv`:

| Order | A3 | A5 | Cf | Nu | Ns(0) | Be(0) |
|---|---|---|---|---|---|---|
| Al2O3-Cu-TiO2 (adopted) | 2.370370 | 2.890824 | 0.232017 | 1.463872 | 0.668808 | 0.947692 |
| Al2O3-TiO2-Cu | 2.370370 | 3.030939 | 0.232017 | 1.470795 | 0.644295 | 0.945702 |
| Cu-Al2O3-TiO2 | 2.370370 | 2.860959 | 0.232017 | 1.462270 | 0.674246 | 0.948114 |
| Cu-TiO2-Al2O3 | 2.370370 | 2.956137 | 0.232017 | 1.467218 | 0.657181 | 0.946767 |
| TiO2-Al2O3-Cu | 2.370370 | 3.116850 | 0.232017 | 1.474586 | 0.630034 | 0.944473 |
| TiO2-Cu-Al2O3 | 2.370370 | 3.068660 | 0.232017 | 1.472500 | 0.637965 | 0.945164 |

`Cf` is identical to six decimal places across all six orderings, as
expected since `A5` does not enter the momentum equation. `Nu` spreads
by only 0.84% (confirming, rather than merely estimating, the
manuscript's "well under 1%" claim). `Ns(0)` is more sensitive, with a
7.02% spread across orderings — this is reported in the manuscript as a
genuine, moderate sensitivity rather than folded into the same
"negligible" bracket as `Cf` and `Nu`.

To reproduce this check:

```
python mixing_order_check.py
```
