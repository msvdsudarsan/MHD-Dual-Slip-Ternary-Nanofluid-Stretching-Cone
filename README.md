# Ternary hybrid nanofluid stretching-cone model — Python code

Python code accompanying:

> S. V. D. S. Madhyannapu, K. Subbarao, "MHD Dual-Slip Flow and Entropy
> Generation in a Cu–Al2O3–TiO2/Water Ternary Hybrid Nanofluid over a
> Stretching Cone," submitted to the *International Journal of Heat and
> Mass Transfer*.

## Contents

- `ternary_cone_master.py` — single-file driver. Computes the three-stage
  sequential thermophysical property ratios (density, viscosity,
  electrical conductivity, heat capacity, thermal conductivity), solves
  the similarity-reduced momentum/energy boundary-value problem with
  `scipy.integrate.solve_bvp`, evaluates the far-field solvability
  criterion (Proposition 1 of the paper), and reproduces every table in
  the manuscript (admissibility map, domain-truncation study,
  mesh-independence study, clear-fluid validation, TiO2 sweep, sigma_TiO2
  sensitivity, shape-factor sensitivity, and the full single-parameter
  sweep of Tables 9-10).
- `mixing_order_check.py` — re-solves the boundary-value problem under
  all six permutations of the Al2O3/Cu/TiO2 mixing order (added in V17;
  see the manuscript's mixing-order discussion in Section 2). Writes
  `data/mixing_order_permutations.csv`.
- `REPRODUCIBILITY_REPORT.md` — a fresh, independently re-run log of the
  script's output, confirming that every number reported in the
  manuscript is reproduced by this code, together with the
  mesh-refinement and clear-fluid-limit checks that establish solution
  accuracy.
- `requirements.txt` — pinned package versions for archival
  reproducibility (added in V17; V16 listed unpinned package names only).

## Requirements

```
pip install -r requirements.txt
```

## Usage

```
python ternary_cone_master.py
python mixing_order_check.py
```

The first command prints the full numerical report to the console and
saves five figures (`velocity_profile.pdf`, `temperature_profile.pdf`,
`entropy_Ns.pdf`, `bejan_number.pdf`, `convergence.pdf`) to the working
directory. On a single core of a 2.1 GHz desktop/cloud-class processor
(Python 3.12, SciPy 1.17) the full run, covering every table and figure
in the paper, completes in under 10 seconds. The second command
re-solves the flow under all six mixing-order permutations and writes
`data/mixing_order_permutations.csv`.

## Citation / permanent archive

This repository accompanies the manuscript's Data Availability
statement. A permanent DOI (via Zenodo or an equivalent archive) has
not yet been minted for this repository as of V17 — see the
manuscript's Data Availability section, which will be updated with the
repository URL and DOI once the repository is made public, following
the same GitHub-Zenodo release workflow used for the authors' companion
rotating-disk-nanofluid paper.

## License / attribution

Released alongside the manuscript for reproducibility. See the paper's
Data Availability statement.
