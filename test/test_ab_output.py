#!/usr/bin/env python3
"""Regression checks for a completed AB run."""

from pathlib import Path

import numpy as np


def read(path: Path) -> np.ndarray:
    return np.fromfile(path, dtype=np.float64)


run = Path("output/baseline_AB")
n = int(np.fromfile(run / "param.dat", dtype=np.int64, count=1)[0])
t = read(run / "t.dat")
r = read(run / "R.dat").reshape(-1, n)
rho_s = read(run / "Rho_s.dat").reshape(-1, n)
rho_d = read(run / "Rho_d.dat").reshape(-1, n)
fractions = read(run / "radii_fraction.dat")
radii_s = read(run / "radii_s.dat").reshape(-1, len(fractions))
radii_d = read(run / "radii_d.dat").reshape(-1, len(fractions))

assert len(t) >= 6
assert np.all(np.diff(t) > 0)
assert np.all(np.isfinite(r)) and np.all(r > 0)
assert np.all(np.diff(r, axis=1) > 0)
assert np.all(np.isfinite(rho_s)) and np.all(rho_s >= 0)
assert np.all(np.isfinite(rho_d)) and np.all(rho_d >= 0)
assert np.all(radii_s[-1] < radii_s[0])
assert np.all(radii_d[-1] > radii_d[0])

for rho in (rho_s, rho_d):
    mass = rho[:, 0] * r[:, 0] ** 3 / 3
    mass += np.sum(rho[:, 1:] * (r[:, 1:] ** 3 - r[:, :-1] ** 3), axis=1) / 3
    assert np.max(np.abs(mass / mass[0] - 1)) < 1e-10

print("AB output regression passed")
