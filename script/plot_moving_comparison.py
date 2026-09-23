#!/usr/bin/env python3
"""Matched-code-time, matched-radius comparison; PDF only by default.

Input is the moving-comparison runner's output directory. Each frame stores
float64 [cell, species, (centre radius, rho, U, enclosed mass, radial velocity)].
Neither time nor amplitude is shifted to align collapse.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COMPONENTS = ("Single stars", "Binaries", "DM")

def load_history(directory):
    data = np.genfromtxt(Path(directory)/"history.csv", delimiter=",", names=True)
    data = np.atleast_1d(data)
    if not np.all(np.isfinite(data["time"])) or np.any(np.diff(data["time"]) <= 0):
        raise ValueError("history times must be finite and strictly increasing")
    return data

def load_frames(directory):
    directory = Path(directory)
    table = np.atleast_1d(np.genfromtxt(directory/"snapshots.csv", delimiter=",", names=True))
    frames = []
    for row in table:
        frame = np.fromfile(directory/f"profile_{int(row['snapshot'])}.dat", dtype="<f8")
        frame = frame.reshape(int(row["zones"]), 3, 5)
        if not np.isfinite(frame).all() or np.any(frame[:, :, :3] <= 0):
            raise ValueError("nonfinite/nonpositive profile")
        if np.any(np.diff(frame[:, :, 0], axis=0) <= 0):
            raise ValueError("unordered profile radii")
        frames.append(frame)
    if np.any(np.diff(table["time"]) <= 0):
        raise ValueError("unordered snapshot times")
    return table["time"], frames

def matched_profile(times, frames, time, radii, component, quantity):
    if time < times[0] or time > times[-1]:
        raise ValueError("refusing temporal extrapolation")
    right = min(int(np.searchsorted(times, time, side="right")), len(times)-1)
    left = max(0, right-1)
    if right == 0 or time == times[right]:
        left = right
    values = []
    for k in (left, right):
        data = frames[k][:, component]
        if radii[0] < data[0, 0] or radii[-1] > data[-1, 0]:
            raise ValueError("refusing spatial extrapolation")
        values.append(np.exp(np.interp(np.log(radii), np.log(data[:, 0]), np.log(data[:, quantity]))))
    weight = 0 if left == right else (time-times[left])/(times[right]-times[left])
    return (1-weight)*values[0]+weight*values[1]

def plot_case(directory):
    from plot_output import read_parameters
    directory = Path(directory)
    histories = [load_history(directory/name) for name in ("hydro", "moving")]
    datasets = [load_frames(directory/name) for name in ("hydro", "moving")]
    final = min(x[0][-1] for x in datasets)
    requested = float(read_parameters(directory/"moving")["maxTime"])
    completed = histories[1]["time"][-1] >= requested
    times = np.array([0, 0.1, 0.25*final, 0.5*final, 0.75*final, final])
    times = np.unique(times[(times >= 0) & (times <= final)])
    colours = plt.cm.rainbow(np.linspace(0, 1, len(times)))
    initial = datasets[0][1][0]
    mass_fraction = initial[-1, :, 3]/initial[-1, :, 3].sum()
    active = np.flatnonzero(mass_fraction > 1e-6)
    fig, axes = plt.subplots(2, len(active), squeeze=False, figsize=(5*len(active), 7), sharex=True)
    for col, f in enumerate(active):
        # Common physical radius range, trimmed for readable core/halo profiles.
        low = max(frame[0, f, 0] for _, frames in datasets for frame in frames)
        high = min(30, *(frame[-1, f, 0] for _, frames in datasets for frame in frames))
        radii = np.geomspace(low, high, 250)
        for t, colour in zip(times, colours):
            for (ts, frames), style in zip(datasets, ("-", "--")):
                rho = matched_profile(ts, frames, t, radii, f, 1)
                sigma = np.sqrt((2/3)*matched_profile(ts, frames, t, radii, f, 2))
                axes[0, col].loglog(radii, rho, style, color=colour, lw=1.4)
                axes[1, col].loglog(radii, sigma, style, color=colour, lw=1.4)
        axes[0, col].set_title(COMPONENTS[f])
        axes[1, col].set_xlabel(r"$r/r_0$ (cell-centre radius)")
        axes[0, col].set_ylabel(r"$\rho/\rho_0$")
        axes[1, col].set_ylabel(r"$\sigma_{1D}/\sqrt{u_0}$")
        for ax in axes[:, col]:
            ax.grid(alpha=.2)
    handles = [Line2D([], [], color=c, label=f"t={t:.3g}") for t, c in zip(times, colours)]
    handles += [Line2D([], [], color="black", ls=s, label=l) for s, l in
                (("-", "Hydrostatic"), ("--", "Moving fluid"))]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5,.97), ncol=4, fontsize=8)
    if not completed:fig.suptitle(f"Partial moving run: t={histories[1]['time'][-1]:.4g} of {requested:g}", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, .87));fig.savefig(directory/"density_dispersion_profiles.pdf");plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    names = ("s", "b", "d")
    metrics = {"common_profile_end_time": float(final), "requested_end_time": requested,
               "moving_reached_requested_time": bool(completed), "components": {}}
    for f in active:
        name = names[f];colour = ("tab:red", "tab:purple", "tab:blue")[f]
        for h, style in zip(histories, ("-", "--")):
            axes[0, 0].semilogy(h["time"], h[f"rho_{name}"], style, color=colour,
                                label=f"{COMPONENTS[f]} / {'hydro' if style=='-' else 'moving'}")
            axes[0, 1].plot(h["time"], np.sqrt(2*h[f"u_{name}"]/3), style, color=colour)
        end = min(h["time"][-1] for h in histories)
        t = np.linspace(0, end, 501)
        ratios = []
        for key, ax, power in ((f"rho_{name}", axes[1, 0], 1), (f"u_{name}", axes[1, 1], .5)):
            a, b = [np.interp(t, h["time"], h[key])**power for h in histories]
            ratio = b/a;ratios.append(ratio)
            ax.plot(t, ratio, color=colour)
            ax.axhline(1, color="grey", lw=.7)
        metrics["components"][name] = {
            "central_density_final_ratio": float(ratios[0][-1]),
            "central_dispersion_final_ratio": float(ratios[1][-1]),
            "fraction_time_density_within_10_percent": float(np.mean(abs(ratios[0]-1) < .1)),
            "fraction_time_dispersion_within_10_percent": float(np.mean(abs(ratios[1]-1) < .1)),
        }
    axes[0, 0].set_ylabel(r"Central $\rho/\rho_0$")
    axes[0, 1].set_ylabel(r"Central $\sigma_{1D}/\sqrt{u_0}$")
    for ax in axes[1]:ax.set_ylabel("Moving / hydrostatic");ax.set_xlabel(r"$t/t_0$")
    for ax in axes.flat:ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=7)
    if not completed:fig.suptitle(f"Partial moving run: t={histories[1]['time'][-1]:.4g} of {requested:g}")
    fig.tight_layout();fig.savefig(directory/"central_evolution.pdf");plt.close(fig)

    # Do not conceal failures in dynamically negligible components by only
    # plotting the mass-dominant species. Preserve raw cell-to-cell velocities.
    fig, axes = plt.subplots(3, 3, figsize=(12, 9))
    for f in range(3):
        low = max(frame[0, f, 0] for _, frames in datasets for frame in frames)
        high = min(3, *(frame[-1, f, 0] for _, frames in datasets for frame in frames))
        radii = np.geomspace(low, high, 400)
        for (ts, frames), style, label in zip(datasets, ("-", "--"), ("Hydrostatic", "Moving")):
            axes[0,f].loglog(radii, matched_profile(ts,frames,final,radii,f,1), style, label=label)
            axes[1,f].loglog(radii, np.sqrt(2*matched_profile(ts,frames,final,radii,f,2)/3), style)
        frame = datasets[1][1][-1][:, f]
        mask = (frame[:,0] >= low) & (frame[:,0] <= high)
        axes[2,f].semilogx(frame[mask,0], frame[mask,4], ".-", ms=2, lw=.7)
        axes[0,f].set_title(COMPONENTS[f])
        axes[0,f].set_ylabel(r"$\rho/\rho_0$");axes[1,f].set_ylabel(r"$\sigma/\sqrt{u_0}$")
        axes[2,f].set_ylabel(r"Moving $v t_0/r_0$");axes[2,f].set_xlabel(r"$r/r_0$")
    axes[0,0].legend(fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle(f"All components: profiles at t={final:.5g}; raw velocity at t={datasets[1][0][-1]:.5g}", fontsize=10)
    fig.tight_layout();fig.savefig(directory/"all_components_final.pdf");plt.close(fig)

    h = histories[1]
    metrics["mass_ledger_relative_error"] = {
        f: float(np.max(np.abs(h[f"m_{f}"]+h[f"out_{f}"]-h[f"m_{f}"][0]))/h[f"m_{f}"][0])
        for f in names}
    metrics["escaped_mass_fraction"] = {f: float(h[f"out_{f}"][-1]/h[f"m_{f}"][0]) for f in names}
    metrics["maximum_linear_backward_error"] = float(np.max(h["linear_residual"]))
    if "energy_ledger_error" in h.dtype.names:
        metrics["maximum_energy_work_ledger_error"] = float(np.max(h["energy_ledger_error"]))
    metrics["hydro_end_time"] = float(histories[0]["time"][-1])
    metrics["moving_end_time"] = float(histories[1]["time"][-1])
    (directory/"comparison_metrics.json").write_text(json.dumps(metrics, indent=2)+"\n")
    return metrics

def plot_split_check(single, split):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    metrics = {}
    for solver, colour in (("hydro", "tab:red"), ("moving", "tab:blue")):
        histories = [load_history(Path(p)/solver) for p in (single, split)]
        time = np.linspace(0, min(h["time"][-1] for h in histories), 501)
        totals, dispersions = [], []
        for h, style, label in zip(histories, ("-", "--"), ("single", "split")):
            rho = sum(h[f"rho_{f}"] for f in ("s", "b", "d"))
            sigma = np.sqrt((2/3)*sum(h[f"rho_{f}"]*h[f"u_{f}"] for f in ("s", "b", "d"))/rho)
            totals.append(np.interp(time, h["time"], rho))
            dispersions.append(np.interp(time, h["time"], sigma))
            axes[0,0].semilogy(time, totals[-1], style, color=colour, label=f"{solver}: {label}")
            axes[0,1].plot(time, dispersions[-1], style, color=colour)
        drho = np.abs(totals[1]/totals[0]-1);dsigma = np.abs(dispersions[1]/dispersions[0]-1)
        axes[1,0].semilogy(time, np.maximum(drho, 1e-16), color=colour)
        axes[1,1].semilogy(time, np.maximum(dsigma, 1e-16), color=colour)
        metrics[solver] = {"maximum_density_relative_difference": float(max(drho)),
                           "maximum_dispersion_relative_difference": float(max(dsigma))}
    axes[0,0].set_ylabel(r"Central total $\rho/\rho_0$")
    axes[0,1].set_ylabel(r"Mass-weighted $\sigma/\sqrt{u_0}$")
    for ax in axes[1]:ax.set_ylabel("Absolute relative difference");ax.set_xlabel(r"$t/t_0$")
    for ax in axes.flat:ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8);fig.tight_layout()
    destination=Path(single).parent
    fig.savefig(destination/"single_split_consistency.pdf");plt.close(fig)
    (destination/"single_split_consistency.json").write_text(json.dumps(metrics, indent=2)+"\n")
    return metrics

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="*", type=Path)
    parser.add_argument("--split-check", nargs=2, type=Path, metavar=("SINGLE", "SPLIT"))
    args = parser.parse_args()
    for directory in args.directories:
        print(directory, json.dumps(plot_case(directory)))
    if args.split_check:print(json.dumps(plot_split_check(*args.split_check)))
    if not args.directories and not args.split_check:parser.error("provide comparison directories or --split-check")
