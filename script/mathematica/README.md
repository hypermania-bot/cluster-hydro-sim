# MovingThreeFluidSim derivation

This derives the experimental C++ moving solver built on upstream
`hypermania/cluster-hydro-sim` revision `49dd74e`. The updated alternating
acoustic flux removes checkerboarding in the tested contraction runs;
coefficient checks alone are not nonlinear stability or convergence evidence.

Open `moving_three_fluid_compiler.nb` and evaluate its first cell. Keep the
notebook and `moving_three_fluid_compiler.wl` in the same directory. The
notebook uses the `MovingThreeFluid` context explicitly; it does not require
variables defined by a previous interactive session.

For a command-line check, from this directory:

```sh
wolframscript -code 'Get["moving_three_fluid_compiler.wl"]; If[TrueQ[MovingThreeFluid`allChecksPassed],Exit[0],Exit[1]]'
```

The source regenerates:

- `moving_three_fluid_coefficients.wl`: full symbolic coefficient/RHS association.
- `moving_three_fluid_entries.txt`: every nonzero scalar coefficient and each RHS.
- `moving_three_fluid_checks.wl`: verification results and numerical derivative errors.
- `moving_three_fluid_storage_coefficients.wl`: all blocks and RHS in optimized C++ order.
- `moving_three_fluid_storage_entries.txt`: every nonzero entry and RHS in that order.

It creates the companion notebook only if missing, preserving later notebook
annotations. The TeX derivation is `docs/moving_three_fluid_scheme.tex`.

## Reading the tables

Cell offsets are -1, 0, +1. Within each cell, species order is s, b, d,
then fields rho, v, u, M. The 12 row order is continuity, momentum, energy,
enclosed-mass constraint, repeated for each species.

Select `rowsByRegion["bulk"]` (in the notebook context), or the keys
`centre`, `outer_fan`, `outer_supersonic`, `outer_vacuum`.
Each result includes `Aminus`, `Adiag`, `Aplus`, `IncrementRHS` and
`AbsoluteRHS`. Unlisted scalar entries are zero. Species can select
different outer regimes independently.

`IncrementRHS` solves A deltaX = b. `AbsoluteRHS` solves for the primitive
predictor X*. The note's local conservative recovery is still required.

The implemented scalar half-bandwidths are 13 lower and 14 upper; the
logical, unpermuted derivation has 19 lower and 13 upper. The complexity claim
is O(N) per step, not O(N) for an entire resolution-refined simulation.

## Implemented choices

- Fixed common geometry; radial fluid motion and a diagnostic moving tidal
  radius, not an unrepresented moving-mesh degree of freedom.
- Regular central parity and the first-cell mass/Robin relation.
- Gas-to-vacuum material boundary and an independently specified absorbing
  Robin thermal boundary, with extrapolation length B in the proposed baseline.
- One linearly implicit solve plus local conservative primitive recovery.
- Current monatomic EOS and c1/c2/c4 thermal closure. Binary mass conversion,
  mean-relative-velocity drag and populated/empty-cell transitions require
  their own subsequent designs/tests; they are not implemented here.

All 37 algebraic checks passed on 22 September 2026. The largest high-precision
finite-difference Jacobian discrepancy in the supplied test states was
9.07e-16. The tests validate algebra and coefficients, not production-run
stability.

## C++ audit and implementation changes

From the repository root, run `make check-moving-symbolic`. This regenerates
the symbolic tables, exports the C++ fixture, compares the two assemblies,
and audits the bandwidth. The independent comparison
includes all matrix entries, RHS entries, centre and all three outer
branches in one five-cell fixture. Matrix/RHS discrepancies are below
5e-16 in the tested build.

The compiler includes the frozen reference acceleration in both momentum
and its mechanical work, and alternating left-velocity/right-pressure
acoustic traces with their pressure work. Raw primitives replace the old
reference offsets. Advective diffusion and velocity viscosity retain their
energy fluxes. Nonlinear finite-difference tests include nonzero reference
acceleration and viscosity. New checks verify flux consistency, absence of
the centred checkerboard null space and the backward-Euler damping polynomial.
The generic symbols
`speed[f,side,0]` and `viscosity[f,side]` are frozen face coefficients;
their numerical construction is given in the TeX note and C++ audit.

## Within-cell bandwidth audit

The original derivation tables retain their original ordering. A simultaneous
row/column permutation to
`[u_s,v_s,rho_s,u_b,v_b,rho_b,u_d,v_d,rho_d,M_s,M_b,M_d]`
reduces `(kl,ku)` from `(19,13)` to `(13,14)`. Equation order becomes
energy, momentum, continuity for each species, followed by the three mass
constraints. Apply the row permutation to the RHS as well.

This is optimal among identical within-cell row/column permutations:
bidirectional primitive coupling in each neighbour block forces
each half-bandwidth to be at least 13. Three primitive variables occupy
positions spanning at least two slots, and every primitive pair couples
in the same direction in both neighbour blocks. Thus one half-bandwidth
must be at least 14. The proposed ordering attains both bounds.
It minimizes `kl+ku` (27), `max(kl,ku)` (14), and LAPACK general-band
factorization storage `2*kl+ku+1` (41, previously 52).
This does not establish optimality for independent row/column permutations
or cross-cell reorderings. Runtime speedup has not been benchmarked.

Reproduce the structural checks and a 7776-candidate block-order search:

```sh
wolframscript -file export_moving_band_pattern.wls
python3 check_moving_bandwidth.py
```

The audit uses the union of the symbolic centre, bulk and all outer regimes;
it does not exploit zeros specific to equilibrium or disabled physics.
