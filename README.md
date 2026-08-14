# Cluster Hydro Simulation

This repository evolves a spherical stellar, binary, and dark-matter system with a conduction/interaction, hydrostatic-relaxation, and grid-realignment scheme. The active executable runs the isolated and tidally truncated Ardi--Baumgardt (AB) comparison cases.

## Build

The Ubuntu build requires GCC with C++20 support, Eigen, LAPACKE, LAPACK, and OpenBLAS. Boost.PFR 1.84 is pinned as a repository submodule.

```bash
git submodule update --init --recursive
make -j6
```

## Run and test

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main
make test
make test-ab
```

The AB runs write packed little-endian `Real64` arrays under `output/baseline_AB/` and `output/with_tidal_AB/`. `make test` checks local interaction-energy and remapping-mass conservation. `make test-ab` also performs a complete AB run and checks the output profiles.

## Plot

`plot_packed.nb` is the original interactive Mathematica analysis. A headless equivalent for the AB section is available when a Mathematica front-end or license is unavailable:

```bash
python3 script/plot_ab.py --input output --output output/plots
```

The script writes density, velocity-dispersion, and Lagrangian-radius plots as PDF and PNG files.

## AB interaction correction

For fluids `f` and `f2`, the dynamical exchange term contains

```text
(m_f U_f - m_f2 U_f2) / m_s.
```

The original band matrix included `m_f / m_s` on the diagonal but omitted `m_f2 / m_s` from the off-diagonal entry. With the AB ratio `m_d / m_s = 0.1` and initially equal specific energies, that implementation heated the dark matter without the compensating stellar cooling, injecting energy and making the stellar profile expand. The corrected matrix conserves interaction energy and produces stellar contraction with dark-matter expansion.

Realignment now also preserves each component's total mass, and its common outer boundary follows the reference algorithm by using the largest component radius.
