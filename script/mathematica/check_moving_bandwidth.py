"""Check within-cell permutations against the exported symbolic sparsity.

Generate moving_band_pattern.json with the companion WL command in the note.
Rows follow the corresponding variable permutation; no equations change.
"""
import itertools
import json
from pathlib import Path
import numpy as np

patterns = json.loads(Path(__file__).with_name("moving_band_pattern.json").read_text())
edges = sorted({(r, c % 12, c // 12 - 1)
                for matrix in patterns.values()
                for r, row in enumerate(matrix)
                for c, value in enumerate(row) if value})
r, c, offset = np.array(edges).T

def bandwidth(order):
    position = np.argsort(order)
    difference = position[r] - position[c] - 12 * offset
    return int(max(difference)), int(max(-difference))

original = tuple(range(12))
recommended = (2, 1, 0, 6, 5, 4, 10, 9, 8, 3, 7, 11)
assert bandwidth(original) == (19, 13)
assert bandwidth(recommended) == (13, 14)

# A bidirectional primitive pair in each neighbour block forces kl,ku >=13.
# All three primitive pairs have at least one directed coupling in BOTH
# neighbour blocks. Their positions span >=2: at least one half-band >=14.
# Thus kl+ku >=27, max(kl,ku)>=14 and 2*kl+ku+1 >=41 for ANY permutation.
edge_set = set(edges)
for f in range(3):
    rho, vel, energy = (4*f+k for k in range(3))
    for side in (-1, 1):
        assert any((a, b, side) in edge_set and (b, a, side) in edge_set
                   for a, b in itertools.combinations((rho, vel, energy), 2))
        for a, b in itertools.combinations((rho, vel, energy), 2):
            assert (a, b, side) in edge_set or (b, a, side) in edge_set
    for a, b in itertools.combinations((rho, vel, energy), 2):
        assert all((a, b, side) in edge_set for side in (-1, 1)) or all(
            (b, a, side) in edge_set for side in (-1, 1))

# Independent search of 7776 natural block orders. The graph lower bound,
# not this restricted enumeration, establishes global optimality.
best = 100
count = 0
for fluids in itertools.permutations(range(3)):
    for local in itertools.product(list(itertools.permutations(range(3))), repeat=3):
        for masses in itertools.permutations((3, 7, 11)):
            order = tuple(4*f+k for f in fluids for k in local[f]) + masses
            kl, ku = bandwidth(order)
            best = min(best, 2*kl+ku+1)
            count += 1
assert count == 7776 and best == 41
print(f"original={bandwidth(original)} recommended={bandwidth(recommended)}")
print(f"searched={count} minimum_LDAB={best}; graph lower bound=41")
print("All bandwidth checks passed.")
