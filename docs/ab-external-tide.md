# AB comparison with an explicit spherical tide

## Implemented force

Set `sim.param.tidal_q` to a positive value to enable the outward acceleration
`f = q*r`; zero (the default) disables it. This is independent of
`tidal_cutoff`, which controls the older explicit density sink. No sink is
enabled in the new AB comparison.

For a dimensionful coefficient Q and fixed HydroSim reference density rho0,

    Phi_ext = -Q*r^2/2
    dP_f/dr = -rho_f*(G*M_total(<r)/r^2 - Q*r)
    tidal_q = Q/(4*pi*G*rho0)

The dimensionless coefficient is not Q times the square of HydroSim's
relaxation-time unit. In code units the effective gravity numerator is
`C_i = M_total[i] - tidal_q*R[f][i]^3`. This is not a change to the actual
enclosed mass. All components experience the same external acceleration.

The existing central, interior and one-sided outer hydrostatic equations
use this numerator. Initialization and pressure reconstruction after
realignment use exactly the same modified equations. The relaxation
Jacobian additionally includes `dC_i/dr_i = -3*tidal_q*r_i^2`, including
the outer row. Shell mass and specific entropy are preserved by the
existing finite-volume updates. No extra irreversible heating is added
to conduction by a static external potential.

`externalPotentialEnergy()` returns the code-unit external energy, integrating
piecewise constant shell densities:

    E_ext = -tidal_q/10 * sum_f,i rho_f,i*(r_i^5-r_(i-1)^5)

There is no additional self-gravity factor of one half. For an isolated
energy budget, add this term to random kinetic plus self-gravitational
energy. Existing conduction-boundary and realignment errors are not repaired
by adding the force.

## What can be inferred from AB?

E. Ardi and H. Baumgardt, *Depletion of dark matter within globular clusters*,
J. Phys.: Conf. Ser. **1503** (2020) 012023,
[DOI](https://doi.org/10.1088/1742-6596/1503/1/012023), Sections 2.1--2.2,
specify NBODY4, mutual gravity plus an external potential, and a circular
orbit at 8.5 kpc around a Milky-Way-like galaxy. They do **not** specify the
galaxy potential, circular speed, NBODY4 tidal switches or escape prescription.
The exact coefficient used in their run cannot be recovered from the paper.

We therefore choose a stated surrogate: a spherical flat-rotation galaxy
with circular speed 220 km/s. In the rotating circular-orbit approximation,
the escape-axis stretching coefficient is

    Q = 4*Omega^2 - kappa^2 = 2*(v_c/R_g)^2.

It is approximately `1.4013e-3 Myr^-2`. Applying this stretching in every
radial direction is not the full NBODY4 tide: it omits vertical compression,
the distinction between radial and tangential directions, Coriolis forces,
and delayed escape. It is also not the angular average of that field.

## Initializer and units

The runner replaces the old AB runner in `main.cpp` and calls the shared
`evolve()`. It uses canonical `initPlummer` and `initCoeffs`, retaining the
old example's equal stellar/DM mass, `rho_s(0)=rho_d(0)=0.5`, Plummer scale
ratios of one, and `m_d/m_s=0.1`. There is a negligible binary trace and no
formation or binary heating. The old hand-adjusted Yiming coefficients are
not used for the new physically scaled experiment.

`ABInitParam` owns physical initialization assumptions; `ABObserverParam`
owns output conversions and the six mass fractions. Both records are saved
separately, along with every evolution parameter including `tidal_q`.

The reference untruncated model has 50,000 particles, equal component masses,
and particle masses of 1 and 0.1 solar masses. Thus

    M = 2*N/(1/m_s + 1/m_d) = 9090.9091 solar masses
    Rh = 3 pc
    a = Rh*sqrt(2^(2/3)-1) = 2.2993 pc
    rho0 = 3*M/(4*pi*a^3), M0 = 3*M
    q = Q*a^3/(3*G*M) = 1.38839708e-4.

The 3 pc half-mass radius is an illustrative physical scale mentioned in
AB, not a specified input for Fig. 1. AB also use a Kroupa stellar mass
spectrum from 1 to 100 solar masses; our one stellar fluid does not reproduce
that spectrum. Matching the particle count does not fix those differences.

For the continuum Plummer profile, the initial force-balance radius satisfies

    rJ/a = sqrt((1/(3*q))^(2/3)-1) = 13.3528.

The initial computational edge is 0.45 of that radius, `Rmax=6.00878`.
Both the isolated and tidal runs use this same edge. `initPlummer` now accepts
an optional outer radius; its old default of 1000 is unchanged. The finite
domain excludes about 4% of the reference Plummer mass and slightly reduces
its represented half-mass radius. Particle masses and the physical conversion
units remain fixed to the reference model, not renormalized to this truncated
or subsequently evolving mass. This initialization truncation is not an
evolving escape prescription.

With `lnLambda = ln(0.8*(M/m_s)/(1+m_d/m_s))`, the initializer computes

    t0 = sqrt(M/(3*m_s))*a^(3/2)/(sqrt(G*m_s)*lnLambda) = 325.296 Myr
    t_fric = 0.035*Rh^(3/2)*sqrt(M/G)/m_s = 258.535 Myr
    t_plot = t_code*t0/t_fric.

The friction time is AB Eq. (2). For the reference Plummer N-body convention
`G=M=-4E=1`, `a_Nbody=3*pi/16`, so the plotted radius is
`r_Nbody=(3*pi/16)*r_code`. These fixed reference conversions are not fits
to the curves. The finite, tidally modified initialization is not claimed
to have exactly the isolated reference binding energy. Radii enclose fractions
1%, 5%, 10%, 20%, 50% and 70% of each component's **current represented mass**.

## Boundary and validity limits

The force-balance radius is not an absorbing boundary. With the existing
one-sided pressure condition, the last finite-density shell has positive
pressure only when `M_total/r^3 > tidal_q`. A positive-pressure hydrostatic
Plummer tail cannot extend to infinity in a uniformly outward `q*r` field.
Initialization or evolution outside the admissible domain therefore throws
an explicit error; it does not clamp the force, erase matter or silently retry.
Only already accepted observations are saved if a run fails.

The existing per-fluid projection freezes enclosed masses by shell index
during each linear solve. The tidal Jacobian is correct for that existing
approximation, not a new fully coupled common-radius gravity Jacobian.
Likewise, the existing fixed outer-temperature conduction equation and
nonconservative density realignment are unchanged. A quantitative tidal
escape calculation requires a separate escape/boundary model.

## Reproduce and plot

Use the strict build for numerical validation:

```sh
make -j3 main-strict check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main-strict ab isolated output/ab/isolated 1000000 4 150
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main-strict ab tidal output/ab/tidal 1000000 4 150
python3 script/plot_ab_comparison.py Ardi_Baumgardt.pdf output/ab/isolated output/ab/tidal --output output/ab/ab_fig1_overlay.pdf
```

The last three runner arguments are maximum steps, final time in friction
units, and zones; they default to 1,000,000, 4 and 150. The initial step is
`1e-4` code time and the shared U-only controller is unchanged. A final step
may overshoot the time limit. The reference scan is supplied by the user;
generated data and PDFs are not committed. The plot uses frame-registered
reference scans with their aspect ratios preserved, not digitized curves.

## Results and checks

All five runs reached four friction times with finite, ordered states.
Representative endpoint results (radii in reference N-body units) are:

| Run | Zones | Stellar mass drift | DM mass drift | Stellar r50 | DM r50 | DM r70 |
|---|---:|---:|---:|---:|---:|---:|
| Isolated | 150 | -0.307% | -11.46% | 0.4692 | 1.4297 | 1.8846 |
| Tidal | 150 | -0.402% | -11.14% | 0.4751 | 1.4499 | 1.9203 |
| Isolated | 300 | -0.084% | -5.71% | 0.4643 | 1.4201 | 1.8677 |
| Tidal | 300 | -0.131% | -5.45% | 0.4699 | 1.4411 | 1.9042 |
| Tidal | 600 | -0.061% | -2.68% | 0.4666 | 1.4341 | 1.8915 |

The DM radii expand and the inner stellar radii contract, but the strong
late-time outer expansion in AB is not reproduced. Here the tidal field
increases the 150-zone final DM r70 by about 1.9% relative to the isolated
control. This is a spherical, finite-domain, no-escape comparison, not a
claim to reproduce AB's unspecified NBODY4 setup. In particular, the listed
mass drift is numerical, **not physical stripping**. Its reduction with grid
refinement exposes the existing remapping limitation; the 600-zone result
still has appreciable DM drift and is not a convergence-certified solution.

`check_tidal_force` independently checks central/interior/outer Jacobian
entries against finite differences (relative discrepancy about `2.5e-9`),
negligible relaxation/realignment change to an already tidal equilibrium,
shell mass and entropy invariants, the external-energy factor, invalid q and
supercritical-boundary rejection, and bounded shared-driver runs. Existing
isolated regressions are retained. Python tests verify the AB conversions,
mass-fraction shapes, plot bounds, preserved scan aspect ratios and PDF output.
