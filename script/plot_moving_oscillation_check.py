#!/usr/bin/env python3
"""Compare raw cell velocities at common code time; write PDF and JSON.

The second-difference estimator includes smooth curvature. It is not a
grid-independent spectral norm, and v/(sqrt(u)) is NOT a Mach number.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_moving_comparison import load_frames, COMPONENTS


def oscillation_metrics(frame, component):
    r, v = frame[:, component, 0], frame[:, component, 4]
    mask = (r[1:-1] > .1) & (r[1:-1] < 1)
    if not np.any(mask):
        raise ValueError("no cells in diagnostic radial window")
    rough = (2*v[1:-1] - v[:-2] - v[2:])/4
    core = v[(r > .1) & (r < 1)]
    slopes = np.diff(core)
    return dict(median=float(np.median(np.abs(rough[mask]))),
                maximum=float(np.max(np.abs(rough[mask]))),
                sign_changes=int(np.count_nonzero(core[:-1]*core[1:] < 0)),
                slope_turns=int(np.count_nonzero(slopes[:-1]*slopes[1:] < 0)),
                full_sign_changes=int(np.count_nonzero(v[:-1]*v[1:] < 0)),
                full_slope_turns=int(np.count_nonzero(np.diff(v)[:-1]*np.diff(v)[1:] < 0)))


def frame_at(times, frames, time):
    if time < times[0] or time > times[-1]:
        raise ValueError("refusing extrapolation")
    hi = min(np.searchsorted(times, time, side="right"), len(times)-1)
    lo = max(0, hi-1)
    if hi == lo:
        return frames[hi]
    if not np.array_equal(frames[lo][:, :, 0], frames[hi][:, :, 0]):
        raise ValueError("diagnostic requires a fixed grid")
    weight = (time-times[lo])/(times[hi]-times[lo])
    return (1-weight)*frames[lo]+weight*frames[hi]


def compare(old, new, destination):
    destination = Path(destination);destination.mkdir(parents=True, exist_ok=True)
    data = [load_frames(path) for path in (old, new)]
    time = min(x[0][-1] for x in data)
    end = [frame_at(*d, time) for d in data]
    report = dict(common_time=float(time), old_end=float(data[0][0][-1]),
                  new_end=float(data[1][0][-1]), components={})
    fig, axes = plt.subplots(3, 2, figsize=(10, 9), squeeze=False)
    for f in range(3):
        report["components"][COMPONENTS[f]] = {}
        for (times, frames), frame, label, colour in zip(data, end,
                ("Centred (old)", "Alternating (new)"), ("tab:red", "tab:blue")):
            report["components"][COMPONENTS[f]][label] = oscillation_metrics(frame, f)
            r = frame[:, f, 0];mask = (r > .08) & (r < 1.2)
            axes[f, 0].plot(r[mask], frame[mask, f, 4], color=colour, label=label, lw=1)
            amplitude = [oscillation_metrics(x, f)["median"] for x in frames]
            axes[f, 1].semilogy(times, np.maximum(amplitude, 1e-16), color=colour, label=label)
        axes[f, 0].set_title(f"{COMPONENTS[f]}, t={time:.5g}")
        axes[f, 0].set_ylabel(r"$v/(r_0/t_0)$")
        axes[f, 1].set_ylabel("Median second-difference estimator")
        for ax in axes[f]:ax.grid(alpha=.2)
    axes[0, 1].legend();axes[-1, 0].set_xlabel(r"$r/r_0$")
    axes[-1, 1].set_xlabel(r"$t/t_0$")
    fig.tight_layout();fig.savefig(destination/"velocity_oscillations.pdf");plt.close(fig)
    (destination/"velocity_metrics.json").write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("old", type=Path);parser.add_argument("new", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args();compare(args.old, args.new, args.destination)
