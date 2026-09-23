# Moving-fluid tidal capture and Statler comparison

`MovingThreeFluidSim` supports `BINARY_FORMATION_OFF` and
`BINARY_FORMATION_POWER_LAW`. The latter uses the same externally initialized
dimensionless coefficient as the existing Statler reproduction:
`S = capture_coefficient * rho_s^2 / U_s^0.6`. No three-body formation,
encounter-driven ejection, or binary disruption is added.

Mass, radial momentum, and random-plus-bulk energy transfer enter the existing
banded solve. Birth velocity is `v_s` and birth random energy is `U_s/2`.
Mass and momentum are conserved locally; capture removes resolved random
energy, tracked by `last_capture_energy` and `cumulative_capture_energy`.
Bulk velocity mixing is handled by conservative recovery. The per-donor
capture frequency is frozen and the donor sink implicit, so this is
first-order equivalent, not bitwise equivalent, to the hydrostatic solver's
explicit formation step. There is no reaction substep or retry loop.

The half-bandwidth remains (13,14). Bounded iterative refinement reuses the
band factorization when needed; the original residual threshold is unchanged.
The state layout, boundary fluxes, and non-formation hydrostatic code are unchanged.

The full local matrix/RHS derivation is in
[moving-binary-formation.tex](moving-binary-formation.tex), with executable
[Mathematica checks](../script/mathematica/moving_capture.wls).

## Reproduction

From the repository root, use fresh output directories:

```sh
make -j2 main-strict check
make check-moving-symbolic
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
./main-strict moving-statler output/moving_statler_500 500 10000 2000000
./main-strict moving-statler output/moving_statler_1000 1000 10000 2000000
python3 script/plot_moving_statler.py path/to/Statler_Ostriker_Cohn.pdf output/moving_statler_500
python3 script/plot_moving_statler.py path/to/Statler_Ostriker_Cohn.pdf output/moving_statler_1000
```

The runner calls `evolve()` and saves evolution, physical initialization,
radial-grid, and observer parameters separately. It uses the paper's
`N_star=3e5`, `M_star=0.7 Msun`, `R_star=0.57 Rsun`, Plummer length `1.13 pc`,
and reference half-mass relaxation time `225 Myr`. The moving grid extends
from its first nonzero face at `1e-6` to `1e4`, unlike the standard initializer's
`1e-2` to `1e3`. An analytic Plummer initializer supplies the pressure and
density on that grid; the standard Statler initializer supplies coefficients.
Tiny binary and inert DM seeds retain the solver's positive-state contract.

Output contains the existing Statler binary-array diagnostic schema, plus
`budget.csv` with formed count, capture-energy loss, total mass including
outflow, and the gas-energy/work residual. An error preserves accepted-history
diagnostics but returns failure; it is not a successful completed run.

Plots are PDFs under `plots/`: one overlay for each of Figs. 11--15 and
`statler_11_15_overlays.pdf`. They reuse the public scan registration, log
axes, limits, and aspect ratios without fitted shifts. `summary.json` reports
the actual endpoint, density peak, binary count, and conservation diagnostics.

## What is and is not matched

This is the direct-heating-only comparison. Initial physical scales, capture
law, and observable units are matched. It is not an implementation of the
paper's orbit-averaged encounter operators: fluid conduction and local c4
heating remain approximations. In particular, the reference suppresses
encounters that would cause ejection; the fluid c4 closure does not resolve
those individual encounters. The fixed coefficient Coulomb logarithms also
do not implement the paper's evolving core-number prescription.

The moving solver still has an open vacuum/thermal outer boundary. Outflow
must be measured, not interpreted as encounter ejection. The energy residual
tests gas energy against discrete gravity work and sources, not conservation
of total gravitational binding energy. The estimated binary age assumes
age-neutral loss; it is not used in Figs. 11--15. Profile legends report the
actual saved epochs. The moving runner saves near `18,50,350,2000 t_rh`,
approximate square locations read from Fig. 11a; the plotting script selects
only profiles within 5% of those targets rather than substituting a distant
endpoint. The paper does not tabulate exact times, so these scan-read epochs
remain approximate and must be refined before a quantitative profile fit.

No Heggie reproduction is claimed by this runner; it requires a separately
specified reference model and settings.
