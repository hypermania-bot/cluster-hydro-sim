# Cluster Hydro Simulation

`cluster-hydro-sim` is a research code for the secular evolution of a
spherically symmetric, self-gravitating system with three components:

- single stars (`s`),
- stellar binaries (`b`), and
- dark matter (`d`).

The code evolves dimensionless density, enclosed mass, pressure, and specific
kinetic energy profiles on radial Lagrangian grids. Its numerical structure
extends the two-fluid conduction method of Zhong and Shapiro to three fluids
and includes experimental binary-formation, binary-heating, and tidal-removal
terms.

The current executable is a research prototype rather than a general-purpose
simulation package. In particular, experiment selection and parameters are
set in `src/main.cpp`, and there is not yet an automated test suite or command
line interface.

## Numerical method

For each fluid `f`, the principal arrays are

```text
R[f][i]     radial shell boundary
Rho[f][i]   shell-averaged mass density
Menc[f][i]  enclosed mass
P[f][i]     pressure
U[f][i]     specific kinetic energy
```

The code uses

```text
P = (2/3) Rho U,        U = (3/2) P/Rho.
```

One evolution step consists of:

1. A semi-implicit conduction and inter-fluid energy-exchange solve. The
   resulting block-banded system is solved with `LAPACKE_dgbsv`.
2. Optional binary formation and tidal mass removal.
3. Two linearized hydrostatic-relaxation solves per component, using
   `LAPACKE_dgtsv` while preserving shell mass and specific entropy within
   each relaxation solve.
4. Logarithmic-grid realignment, followed by reconstruction of enclosed mass
   and hydrostatic pressure.

The center uses a regularity boundary condition, `U[f][0] = U[f][1]`, and the
outer value of `U` is held fixed during the conduction solve. The timestep
starts at `1e-3` and is adjusted using a target maximum relative change in
`U` of `1e-3`. Evolution ends when a component's central density exceeds
`1e12` or the configured step limit is reached.

The implementation and notation follow:

- Yi-Ming Zhong and Stuart L. Shapiro, [*Dynamical Evolutions in Globular
  Clusters and Dwarf Galaxies: Conduction Fluid
  Simulations*](https://arxiv.org/abs/2505.18251).

## Repository layout

```text
.
|-- Makefile             GCC/C++20 release build
|-- external/            pinned Eigen and Boost.PFR submodules
|-- plot_packed.nb       Mathematica import and plotting notebook
`-- src/
    |-- main.cpp         experiment definitions and active entry point
    |-- three_fluid.*    simulation state and numerical operators
    |-- io.hpp           binary array I/O
    |-- param.hpp        parameter serialization for Mathematica
    `-- utility.hpp      output-directory and timing helpers
```

## Requirements

The tested Ubuntu build uses:

- GCC with C++20 support,
- GNU Make,
- Eigen,
- Boost.PFR,
- LAPACKE and LAPACK, and
- OpenBLAS.

Eigen and Boost.PFR are pinned as Git submodules. On Ubuntu, the remaining
development packages can be installed with

```bash
sudo apt-get install build-essential liblapacke-dev liblapack-dev libopenblas-dev
```

## Clone and build

Clone the repository together with its submodules:

```bash
git clone --recurse-submodules https://github.com/hypermania/cluster-hydro-sim.git
cd cluster-hydro-sim
make -j6
```

If the repository was cloned without submodules, initialize them before
building:

```bash
git submodule update --init --recursive
make -j6
```

The compiler can be selected through Make's standard `CXX` variable, for
example `make CXX=g++ -j6`. The default build enables `-O3`, native CPU
instructions, OpenMP, `-ffast-math`, and `NDEBUG`, so its executable is not
portable across all CPU architectures and does not retain Eigen's debug
assertions.

## Run

Run from the repository root so the relative output paths resolve correctly:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./main
```

The active entry point runs two one-fluid tidal comparison cases:

| Output directory | Tidal removal | Tidal radius | Removal factor |
| --- | ---: | ---: | ---: |
| `output/one_fluid_without_tidal/` | off | -- | -- |
| `output/one_fluid_with_tidal/` | on | `10` | `50` |

Both cases use 150 zones, `m_d = m_s`, and initial central-density ratios
`rho_s:rho_b:rho_d = 1:1e-10:1e-10`. The dark-matter conduction coefficient
is set equal to the single-star coefficient so that the negligible binary and
dark-matter populations leave an effectively one-fluid calculation. The
tidal case applies the explicit density sink outside `r = 10`; the no-tidal
case additionally records the central density and specific energy at every
step. Binary conduction, binary heating, and binary formation remain
disabled. Other experiment functions are present in `src/main.cpp` but are
not selected by a runtime option.

On an Ubuntu 22 system with GCC 11, the current default build completed both
cases in approximately 15 seconds. The no-tidal case stopped at `t = 5.56517`
after 15,795 steps, and the tidal case stopped at `t = 5.53821` after 40,135
steps. These values are reference observations, not regression-test
tolerances, and may vary with the compiler and target architecture.

`make clean` removes the executable and object files. It deliberately leaves
simulation output in place.

## Output format

A run directory contains parameter metadata and the packed, headerless binary
arrays produced by its attached observers:

| File | Contents | Shape |
| --- | --- | --- |
| `t.dat` | profile snapshot times | `n_snapshot` |
| `R.dat` | common radial grid | `n_snapshot x N` |
| `Rho_s.dat`, `Rho_b.dat`, `Rho_d.dat` | component densities | `n_snapshot x N` |
| `U_s.dat`, `U_b.dat`, `U_d.dat` | component specific energies | `n_snapshot x N` |
| `radii_t.dat` | Lagrangian-radius sample times | `n_radius_sample` |
| `radii_fraction.dat` | enclosed-mass fractions | `n_fraction` |
| `radii_s.dat`, `radii_b.dat`, `radii_d.dat` | component Lagrangian radii | `n_radius_sample x n_fraction` |
| `central_t.dat` | central-value sample times | `n_state` |
| `central_Rho_s.dat`, `central_Rho_b.dat`, `central_Rho_d.dat` | component central densities | `n_state` |
| `central_U_s.dat`, `central_U_b.dat`, `central_U_d.dat` | component central specific energies | `n_state` |

Numeric arrays are native binary `double` values; on the tested x86-64 system
they are little-endian IEEE-754 `float64`. Snapshot times are approximate:
the observer records the first evolved state whose time reaches each requested
value. If a run ends before a requested time, that snapshot is absent. The
Lagrangian-radius observer, by contrast, records every evolution step and can
therefore consume substantial memory in a long run.

The `central_*` files are currently written only by the no-tidal run. Its
central-value observer is called at the same cadence as the Lagrangian-radius
observer, including the initial and terminal states, so `central_t.dat` and
`radii_t.dat` have the same length and timestamps for that run.

The parameter metadata consists of `param.dat`, `paramNames.txt`,
`paramTypes.txt`, and `paramOffsets.txt`. `param.dat` is a raw C++
`ThreeFluidParam` object, including implementation-dependent layout and
padding; use the accompanying type and offset files when importing it.

For example, NumPy can load one run as follows:

```python
from pathlib import Path
import numpy as np

run = Path("output/baseline_AB")
N = int(np.fromfile(run / "param.dat", dtype="<i8", count=1)[0])
t = np.fromfile(run / "t.dat", dtype="<f8")
r = np.fromfile(run / "R.dat", dtype="<f8").reshape(len(t), N)
rho_s = np.fromfile(run / "Rho_s.dat", dtype="<f8").reshape(len(t), N)
u_s = np.fromfile(run / "U_s.dat", dtype="<f8").reshape(len(t), N)
```

`plot_packed.nb` imports the parameter metadata and packed arrays and contains
sections for the one-fluid split, single-fluid tidal comparison, and AB
comparison. It requires an activated Mathematica installation.

## Current validation limits

Before using the results for scientific inference, add invariant and
convergence checks appropriate to the experiment. Important current limits
include:

- the log-density interpolation used during realignment is not
  mass-conservative;
- there are no automated unit, regression, or grid-convergence tests;
- binary formation still contains an unresolved dimensionless-conversion
  note and is disabled in the active setup;
- the tidal prescription is a phenomenological explicit density sink rather
  than an external gravitational potential; and
- output writes do not currently report short writes or file-open failures.

The repository does not currently declare a software license.
