# Mode-2 binary-formation validation

This note records a conservation assessment of the checked-in mode-2 example
at code commit `02bd417`. We ran `binary_formation()` from `src/main.cpp`
without changing its initialization, parameters, operator ordering, stopping
condition, or production compiler flags.

## Run outcome

The production executable first runs the no-formation baseline and then mode
2. Both cases complete normally:

```text
baseline: t=5.73374, steps=19642
mode 2:   t=5.751016726960865, steps=42567
```

Mode 2 stops when the maximum central density reaches
`1.0003593735937466e12`, just above `StopDensity=1e12`. It does not stop due
to a NaN, a nonpositive field, shell crossing, or the one-million-step limit.
The proposed next timestep is `2.3146e-12`.

Across the trajectory, the minimum density is `1.1360e-25`, the minimum `U`
is `3.7184e-5`, and all component grids remain ordered. The maximum relative
change accepted from one conduction solve is `5.9944e-3`. This exceeds the
nominal `1e-3` target because the current controller uses a completed step to
choose the next timestep rather than rejecting and retrying it.

## Conserved quantities

We measure component mass with the code's shell convention,

```text
M_f = sum_i Rho[f][i] (R[f][i]^3 - R[f][i-1]^3)/3,
```

where the inner radius of zone zero is zero. The resolved energy is

```text
E_resolved = sum_f integral Rho_f U_f r^2 dr
             - integral M(r) dM(r)/r.
```

The gravitational term is integrated exactly for the piecewise-constant
density represented by the numerical shells. It is not an accounting of
binary binding energy or energy exchanged with the prescribed outer thermal
boundary.

To locate conservation errors, we evaluate mass and energy immediately before
and after each split operator:

```text
conduction -> mode-2 formation -> relaxation -> realignment.
```

The resulting ledgers close to `6.1e-20` in mass and `1.3e-20` in energy,
relative to their initial scales.

## Mass ledger

The total mass changes as

```text
initial total mass       0.3333326863622548
final total mass         0.2751260660918499
relative change         -0.1746201997338716
```

| Operator | Mass change / initial mass |
|---|---:|
| Conduction | `0` |
| Mode-2 formation | `+1.1634e-13` |
| Relaxation | `-9.3635e-16` |
| Realignment | `-0.1746201997340` |

Mode-2 formation transfers `7.8136882e-7` from singles to binaries. Its
maximum per-call total-mass residual is `2.342e-17` of the initial mass, and
the maximum pointwise residual in total density is `2.220e-16`. Formation is
therefore conservative to roundoff.

Realignment removes `0.05820662` of total mass, approximately 74,493 times the
mass transferred by formation. It also removes `7.3789e-9`, or about 0.94%,
of the binary mass created during the trajectory. The final binary mass is
`7.7402331e-7`.

## Energy ledger

The initial and final resolved energies are

```text
                              initial                final
random energy          +0.0164627327038830   +0.0147443010305511
gravitational energy   -0.0326995094245412   -0.0289601379892112
total energy           -0.0162367767206582   -0.0142158369586601
```

The total energy increases by `0.002020939762`, or `0.124467` times the
initial absolute total energy. The system therefore becomes less bound.

| Operator | Energy change / absolute initial energy |
|---|---:|
| Conduction, exchange, and binary heating | `+6.56648e-4` |
| Mode-2 formation | `-5.01739e-13` |
| Relaxation | `+4.07272e-4` |
| Realignment | `+0.123402888` |

Formation's maximum per-call random-energy residual is `5.248e-17` of the
initial energy scale, and its maximum resolved-total-energy residual is
`1.195e-16`. This agrees with the pointwise update: total density is unchanged,
so gravitational energy is unchanged, and the updated binary `U` preserves
`Rho_s U_s + Rho_b U_b`.

Binary heating injects resolved energy without evolving a binary-binding-
energy reservoir. The fixed outer-`U` row can also exchange heat with the
system. These effects prevent the resolved energy above from being a complete
physical conservation law, but they are small in this run compared with the
realignment contribution.

## Source stiffness

The largest formation updates are

```text
maximum global single-mass depletion per call   7.2250e-10
maximum local single-density depletion           1.5955e-4
maximum local binary-density growth               7.3535e-4
```

The source remains positive in this trajectory. It is nevertheless absent
from timestep selection. A robust controller should bound local donor
depletion and recipient growth, reject a step that violates those bounds, and
retry with a smaller timestep.

## Interpretation

Mode-2 formation itself passes the local mass and resolved-energy checks. The
complete trajectory does not: logarithmic interpolation in `realign()`
dominates both the mass loss and resolved-energy change. A conservative
finite-volume remap and a grid-convergence test are required before the
mode-2 trajectory can be interpreted quantitatively.
