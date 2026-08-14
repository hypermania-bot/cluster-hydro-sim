#!/usr/bin/env python3
"""Plot packed AB outputs using the layout and labels from plot_packed.nb."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_real64(path: Path) -> np.ndarray:
    return np.fromfile(path, dtype=np.float64)


def read_run(path: Path) -> dict[str, np.ndarray | int]:
    n = int(np.fromfile(path / "param.dat", dtype=np.int64, count=1)[0])
    fractions = read_real64(path / "radii_fraction.dat")
    return {
        "n": n,
        "t": read_real64(path / "t.dat"),
        "r": read_real64(path / "R.dat").reshape(-1, n),
        "rho_s": read_real64(path / "Rho_s.dat").reshape(-1, n),
        "rho_d": read_real64(path / "Rho_d.dat").reshape(-1, n),
        "u_s": read_real64(path / "U_s.dat").reshape(-1, n),
        "u_d": read_real64(path / "U_d.dat").reshape(-1, n),
        "fractions": fractions,
        "radii_t": read_real64(path / "radii_t.dat"),
        "radii_s": read_real64(path / "radii_s.dat").reshape(-1, len(fractions)),
        "radii_d": read_real64(path / "radii_d.dat").reshape(-1, len(fractions)),
    }


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Latin Modern Roman", "Computer Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "font.size": 13,
            "axes.labelsize": 17,
            "axes.linewidth": 1.5,
            "legend.fontsize": 9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.bbox": "tight",
        }
    )


def plot_profiles(data: dict[str, np.ndarray | int], output: Path) -> None:
    t = data["t"]
    r = data["r"]
    rho_s = data["rho_s"]
    rho_d = data["rho_d"]
    u_s = data["u_s"]
    u_d = data["u_d"]
    assert isinstance(t, np.ndarray)
    assert isinstance(r, np.ndarray)
    assert isinstance(rho_s, np.ndarray)
    assert isinstance(rho_d, np.ndarray)
    assert isinstance(u_s, np.ndarray)
    assert isinstance(u_d, np.ndarray)

    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(t)))
    fig, ax = plt.subplots(figsize=(6, 6))
    for i, (time, color) in enumerate(zip(t, colors)):
        valid_s = rho_s[i] > 0
        valid_d = rho_d[i] > 0
        ax.loglog(r[i, valid_s], rho_s[i, valid_s], color=color, lw=1.8,
                  label=rf"stars, $\hat t={time:.3g}$")
        ax.loglog(r[i, valid_d], rho_d[i, valid_d], color=color, lw=1.5, ls="--",
                  label=rf"DM, $\hat t={time:.3g}$")
    ax.set_xlim(1e-3, 1e3)
    ax.set_ylim(1e-15, 1e3)
    ax.set_xlabel(r"$\hat r$")
    ax.set_ylabel(r"$\hat\rho$")
    ax.text(0.04, 0.05, "AB parameters, isolated", transform=ax.transAxes, fontsize=14)
    ax.legend(loc="upper right", ncol=2, frameon=False)
    fig.savefig(output / "AB_rho.pdf")
    fig.savefig(output / "AB_rho.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    for i, (time, color) in enumerate(zip(t, colors)):
        sigma_s = np.sqrt(2.0 * u_s[i] / 3.0)
        sigma_d = np.sqrt(2.0 * u_d[i] / 3.0)
        ax.loglog(r[i], sigma_s, color=color, lw=1.8,
                  label=rf"stars, $\hat t={time:.3g}$")
        ax.loglog(r[i], sigma_d, color=color, lw=1.5, ls="--",
                  label=rf"DM, $\hat t={time:.3g}$")
    ax.set_xlim(1e-3, 1e3)
    ax.set_ylim(0.03, 1.5)
    ax.set_xlabel(r"$\hat r$")
    ax.set_ylabel(r"$\hat\sigma$")
    ax.legend(loc="upper right", ncol=2, frameon=False)
    fig.savefig(output / "AB_sigma.pdf")
    fig.savefig(output / "AB_sigma.png", dpi=180)
    plt.close(fig)


def plot_radii(
    isolated: dict[str, np.ndarray | int],
    tidal: dict[str, np.ndarray | int],
    output: Path,
) -> None:
    fractions = isolated["fractions"]
    assert isinstance(fractions, np.ndarray)
    colors = plt.cm.plasma(np.linspace(0.08, 0.9, len(fractions)))
    fig, ax = plt.subplots(figsize=(7.2, 6))
    for i, (fraction, color) in enumerate(zip(fractions, colors)):
        ax.semilogy(isolated["radii_t"], isolated["radii_s"][:, i], color=color,
                    lw=1.8, label=rf"isolated ({fraction:g})")
        ax.semilogy(tidal["radii_t"], tidal["radii_s"][:, i], color=color,
                    lw=1.3, ls="--", label=rf"tidal ({fraction:g})")
    ax.set_xlabel(r"$\hat t$")
    ax.set_ylabel(r"$\hat r$")
    ax.set_title("Stellar Lagrangian radii (AB)")
    ax.legend(loc="best", ncol=2, frameon=False)
    fig.savefig(output / "AB_radii.pdf")
    fig.savefig(output / "AB_radii.png", dpi=180)
    plt.close(fig)


def plot_before_after(
    before: dict[str, np.ndarray | int],
    after: dict[str, np.ndarray | int],
    output: Path,
) -> None:
    fractions = after["fractions"]
    assert isinstance(fractions, np.ndarray)
    selected = [0, 3, 5]
    colors = plt.cm.tab10(np.arange(len(selected)))
    fig, ax = plt.subplots(figsize=(7.2, 6))
    for color, i in zip(colors, selected):
        fraction = fractions[i]
        ax.semilogy(before["radii_t"], before["radii_s"][:, i], color=color,
                    lw=1.4, ls="--", label=rf"before, $f={fraction:g}$")
        ax.semilogy(after["radii_t"], after["radii_s"][:, i], color=color,
                    lw=2.0, label=rf"fixed, $f={fraction:g}$")
    ax.set_xlim(0, min(3.2, float(before["radii_t"][-1])))
    ax.set_xlabel(r"$\hat t$")
    ax.set_ylabel(r"$\hat r$")
    ax.set_title("AB stellar profile: interaction bug before and after")
    ax.legend(loc="best", ncol=2, frameon=False)
    fig.savefig(output / "AB_radii_before_after.pdf")
    fig.savefig(output / "AB_radii_before_after.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("output"))
    parser.add_argument("--output", type=Path, default=Path("output/plots"))
    parser.add_argument("--before", type=Path)
    args = parser.parse_args()

    configure_style()
    args.output.mkdir(parents=True, exist_ok=True)
    isolated = read_run(args.input / "baseline_AB")
    tidal = read_run(args.input / "with_tidal_AB")
    plot_profiles(isolated, args.output)
    plot_radii(isolated, tidal, args.output)
    if args.before is not None:
        plot_before_after(read_run(args.before / "baseline_AB"), isolated, args.output)
    print(f"Wrote AB plots to {args.output}")


if __name__ == "__main__":
    main()
