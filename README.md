# Ternary hybrid nanofluid stretching-cone model — Python code

Python code accompanying:

> S. V. D. S. Madhyannapu, K. Subbarao, P. K. Cintaginjala, A. Savitha Mary,
> G. Arunagiri, V. Sravana Kumar, "Nonlinear MHD Dual-Slip Flow and Entropy
> Generation in a Cu-Al2O3-TiO2/Water Ternary Hybrid Nanofluid over a
> Stretching Cone," submitted to the *International Journal of Heat and
> Mass Transfer*.

## Contents

- `ternary_cone_master.py` — single-file driver. Computes the three-stage
  sequential thermophysical property ratios (density, viscosity,
  electrical conductivity, heat capacity, thermal conductivity), solves
  the similarity-reduced momentum/energy boundary-value problem with
  `scipy.integrate.solve_bvp`, evaluates the far-field solvability
  criterion, and reproduces every table in the manuscript (admissibility
  map, domain-truncation study, mesh-independence study, clear-fluid
  validation, TiO2 sweep, sigma_TiO2 sensitivity, shape-factor
  sensitivity, and the full single-parameter parametric sweep).
- `mixing_order_check.py` — re-solves the boundary-value problem under
  all six permutations of the Al2O3/Cu/TiO2 mixing order. Writes
  `data/mixing_order_permutations.csv`.
- `make_combined_figures.py` — combines the velocity and temperature
  profiles into one 2-panel figure, and the local entropy generation and
  Bejan number profiles into another, using the same `solve_full()`
  baseline/binary computation as `make_figures()` in
  `ternary_cone_master.py`.
- `REPRODUCIBILITY_REPORT.md` — a log of the script's output, confirming
  that every number reported in the manuscript is reproduced by this
  code, together with the mesh-refinement and clear-fluid-limit checks
  that establish solution accuracy.
- `requirements.txt` — pinned package versions for archival
  reproducibility.

## Requirements

```
pip install -r requirements.txt
```

## Usage

```
python ternary_cone_master.py
python mixing_order_check.py
python make_combined_figures.py
```

The first command prints the full numerical report to the console and
saves five figures (`velocity_profile.pdf`, `temperature_profile.pdf`,
`entropy_Ns.pdf`, `bejan_number.pdf`, `convergence.pdf`) to the working
directory. The second command re-solves the flow under all six
mixing-order permutations and writes `data/mixing_order_permutations.csv`.
The third command produces `velocity_temperature.pdf` and
`entropy_bejan.pdf`, the two combined figures embedded in the manuscript.

## Citation / permanent archive

This repository accompanies the manuscript's Data Availability statement
and is published at:
https://github.com/msvdsudarsan/MHD-Dual-Slip-Ternary-Nanofluid-Stretching-Cone

Release v1.0.0 of this repository is archived on Zenodo with the
permanent identifier:
https://doi.org/10.5281/zenodo.22048410

## License / attribution

Released alongside the manuscript for reproducibility. See the paper's
Data Availability statement.
