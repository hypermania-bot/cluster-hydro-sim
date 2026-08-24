#!/usr/bin/env python3
"""Plot the packed binary output produced by cluster-hydro-sim.

The visual conventions follow ``plot_packed.nb``: dimensionless hatted axis
labels, log-log density and velocity-dispersion profiles, log-radius
Lagrangian evolution, component/time legends, and Latin Modern typography.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COMPONENTS = ("s", "b", "d")
COMPONENT_LABELS = {
    "s": "star",
    "b": "binary",
    "d": "DM",
}
COMPONENT_STYLES = {
    "s": {"color": "#aa0000", "linestyle": "--"},
    "b": {"color": "#006400", "linestyle": "-."},
    "d": {"color": "#666666", "linestyle": "-"},
}
PARAM_DTYPES = {
    "Integer64": np.dtype("<i8"),
    "Real64": np.dtype("<f8"),
}


@dataclass(frozen=True)
class ProfileData:
    time: np.ndarray
    radius: np.ndarray
    rho: dict[str, np.ndarray]
    energy: dict[str, np.ndarray]


@dataclass(frozen=True)
class RadiiData:
    time: np.ndarray
    fraction: np.ndarray
    radius: dict[str, np.ndarray]


@dataclass(frozen=True)
class CentralData:
    time: np.ndarray
    rho: dict[str, np.ndarray]
    energy: dict[str, np.ndarray]


@dataclass(frozen=True)
class RunData:
    directory: Path
    parameters: dict[str, object]
    profiles: ProfileData | None
    radii: RadiiData | None
    central: CentralData | None


class OutputFormatError(ValueError):
    """Raised when packed output is absent, incomplete, or malformed."""


def _read_lines(path: Path) -> list[str]:
    if not path.is_file():
        raise OutputFormatError(f"missing required metadata file: {path}")
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def read_parameters(directory: Path) -> dict[str, object]:
    """Read ``param.dat`` using the accompanying names, types, and offsets."""
    names = _read_lines(directory / "paramNames.txt")
    types = _read_lines(directory / "paramTypes.txt")
    offset_text = _read_lines(directory / "paramOffsets.txt")
    param_path = directory / "param.dat"
    if not param_path.is_file():
        raise OutputFormatError(f"missing required metadata file: {param_path}")

    if not (len(names) == len(types) == len(offset_text)):
        raise OutputFormatError(
            "parameter metadata lengths differ: "
            f"{len(names)} names, {len(types)} types, {len(offset_text)} offsets"
        )
    try:
        offsets = [int(value) for value in offset_text]
    except ValueError as exc:
        raise OutputFormatError("parameter offsets must be decimal byte offsets") from exc
    if offsets != sorted(offsets) or (offsets and offsets[0] < 0):
        raise OutputFormatError("parameter offsets must be nonnegative and ordered")

    raw = param_path.read_bytes()
    parameters: dict[str, object] = {}
    for index, (name, type_name, start) in enumerate(zip(names, types, offsets)):
        try:
            dtype = PARAM_DTYPES[type_name]
        except KeyError as exc:
            raise OutputFormatError(
                f"unsupported parameter type {type_name!r} for {name!r}"
            ) from exc
        stop = offsets[index + 1] if index + 1 < len(offsets) else len(raw)
        if start >= stop or stop > len(raw):
            raise OutputFormatError(
                f"invalid byte range [{start}, {stop}) for parameter {name!r}"
            )
        byte_count = stop - start
        count = byte_count // dtype.itemsize
        if count == 0:
            raise OutputFormatError(f"parameter {name!r} is shorter than {dtype.itemsize} bytes")
        values = np.frombuffer(raw, dtype=dtype, count=count, offset=start).copy()
        parameters[name] = values.item() if values.size == 1 else values
    return parameters


def _read_float64(path: Path) -> np.ndarray:
    if not path.is_file():
        raise OutputFormatError(f"missing required data file: {path}")
    if path.stat().st_size % 8:
        raise OutputFormatError(f"{path} size is not a multiple of eight bytes")
    return np.fromfile(path, dtype="<f8")


def _optional_family(directory: Path, filenames: list[str]) -> bool:
    present = [name for name in filenames if (directory / name).is_file()]
    if present and len(present) != len(filenames):
        missing = sorted(set(filenames) - set(present))
        raise OutputFormatError(
            f"incomplete data family in {directory}; missing: {', '.join(missing)}"
        )
    return bool(present)


def _reshape_exact(values: np.ndarray, shape: tuple[int, ...], path: Path) -> np.ndarray:
    expected = int(np.prod(shape))
    if values.size != expected:
        raise OutputFormatError(
            f"{path} has {values.size} values, expected {expected} for shape {shape}"
        )
    return values.reshape(shape)


def _warn_nonfinite(name: str, values: np.ndarray) -> None:
    count = int(values.size - np.count_nonzero(np.isfinite(values)))
    if count:
        warnings.warn(f"{name} contains {count} non-finite values; they will be omitted")


def load_output(directory: str | Path) -> RunData:
    """Load and validate every complete observer family in an output directory."""
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise OutputFormatError(f"output directory does not exist: {directory}")
    parameters = read_parameters(directory)

    profiles: ProfileData | None = None
    profile_files = ["t.dat", "R.dat"] + [
        f"{field}_{component}.dat"
        for field in ("Rho", "U")
        for component in COMPONENTS
    ]
    if _optional_family(directory, profile_files):
        if "N" not in parameters:
            raise OutputFormatError("parameter metadata does not contain N")
        n_grid = int(parameters["N"])
        if n_grid <= 0:
            raise OutputFormatError(f"N must be positive, got {n_grid}")
        time = _read_float64(directory / "t.dat")
        if time.size == 0:
            raise OutputFormatError("t.dat contains no snapshots")
        shape = (time.size, n_grid)
        radius = _reshape_exact(_read_float64(directory / "R.dat"), shape, directory / "R.dat")
        rho = {
            component: _reshape_exact(
                _read_float64(directory / f"Rho_{component}.dat"),
                shape,
                directory / f"Rho_{component}.dat",
            )
            for component in COMPONENTS
        }
        energy = {
            component: _reshape_exact(
                _read_float64(directory / f"U_{component}.dat"),
                shape,
                directory / f"U_{component}.dat",
            )
            for component in COMPONENTS
        }
        _warn_nonfinite("R", radius)
        for component in COMPONENTS:
            _warn_nonfinite(f"Rho_{component}", rho[component])
            _warn_nonfinite(f"U_{component}", energy[component])
        profiles = ProfileData(time=time, radius=radius, rho=rho, energy=energy)

    radii: RadiiData | None = None
    radii_files = ["radii_t.dat", "radii_fraction.dat"] + [
        f"radii_{component}.dat" for component in COMPONENTS
    ]
    if _optional_family(directory, radii_files):
        time = _read_float64(directory / "radii_t.dat")
        fraction = _read_float64(directory / "radii_fraction.dat")
        if time.size == 0 or fraction.size == 0:
            raise OutputFormatError("Lagrangian-radius times and fractions must be nonempty")
        shape = (time.size, fraction.size)
        radius = {
            component: _reshape_exact(
                _read_float64(directory / f"radii_{component}.dat"),
                shape,
                directory / f"radii_{component}.dat",
            )
            for component in COMPONENTS
        }
        for component in COMPONENTS:
            _warn_nonfinite(f"radii_{component}", radius[component])
        radii = RadiiData(time=time, fraction=fraction, radius=radius)

    central: CentralData | None = None
    central_files = ["central_t.dat"] + [
        f"central_{field}_{component}.dat"
        for field in ("Rho", "U")
        for component in COMPONENTS
    ]
    if _optional_family(directory, central_files):
        time = _read_float64(directory / "central_t.dat")
        if time.size == 0:
            raise OutputFormatError("central_t.dat contains no states")
        rho = {
            component: _reshape_exact(
                _read_float64(directory / f"central_Rho_{component}.dat"),
                (time.size,),
                directory / f"central_Rho_{component}.dat",
            )
            for component in COMPONENTS
        }
        energy = {
            component: _reshape_exact(
                _read_float64(directory / f"central_U_{component}.dat"),
                (time.size,),
                directory / f"central_U_{component}.dat",
            )
            for component in COMPONENTS
        }
        for component in COMPONENTS:
            _warn_nonfinite(f"central_Rho_{component}", rho[component])
            _warn_nonfinite(f"central_U_{component}", energy[component])
        central = CentralData(time=time, rho=rho, energy=energy)

    if profiles is None and radii is None and central is None:
        raise OutputFormatError(f"no observer data found in {directory}")
    return RunData(
        directory=directory,
        parameters=parameters,
        profiles=profiles,
        radii=radii,
        central=central,
    )


def configure_style() -> None:
    """Apply the typography and framed-axis style used by plot_packed.nb."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Latin Modern Roman", "CMU Serif", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "font.size": 14,
            "axes.labelsize": 18,
            "axes.titlesize": 18,
            "axes.edgecolor": "black",
            "axes.linewidth": 1.0,
            "legend.fontsize": 10,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.facecolor": "white",
        }
    )


def _snapshot_indices(count: int, maximum: int) -> np.ndarray:
    if maximum <= 0:
        raise ValueError("maximum number of profile snapshots must be positive")
    if count <= maximum:
        return np.arange(count, dtype=int)
    return np.unique(np.linspace(0, count - 1, maximum, dtype=int))


def _plot_positive(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    *,
    require_positive_x: bool = True,
    **kwargs: object,
) -> bool:
    mask = np.isfinite(x) & np.isfinite(y) & (y > 0.0)
    if require_positive_x:
        mask &= x > 0.0
    if np.count_nonzero(mask) < 2:
        return False
    ax.plot(x[mask], y[mask], **kwargs)
    return True


def _profile_figure(
    profiles: ProfileData,
    quantity: str,
    title: str,
    maximum_snapshots: int,
) -> plt.Figure:
    indices = _snapshot_indices(profiles.time.size, maximum_snapshots)
    fig, ax = plt.subplots(figsize=(8.0, 8.0))
    for rank, index in enumerate(indices):
        opacity = 0.48 + 0.52 * (rank + 1) / len(indices)
        for component in COMPONENTS:
            values = (
                profiles.rho[component][index]
                if quantity == "density"
                else np.sqrt((2.0 / 3.0) * profiles.energy[component][index])
            )
            style = COMPONENT_STYLES[component]
            _plot_positive(
                ax,
                profiles.radius[index],
                values,
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=1.7,
                alpha=opacity,
                label=(
                    f"{COMPONENT_LABELS[component]} "
                    rf"$\hat{{t}}={profiles.time[index]:.4f}$"
                ),
            )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\hat{r}$")
    ax.set_ylabel(r"$\hat{\rho}$" if quantity == "density" else r"$\hat{\sigma}$")
    ax.set_title(title)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="best", frameon=False, ncol=2)
    fig.tight_layout()
    return fig


def _radii_figure(radii: RadiiData, title: str) -> plt.Figure:
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 12.0), sharex=True)
    line_styles = ("-", "--", "-.", ":")
    for ax, component in zip(axes, COMPONENTS):
        for index, fraction in enumerate(radii.fraction):
            _plot_positive(
                ax,
                radii.time,
                radii.radius[component][:, index],
                color=COMPONENT_STYLES[component]["color"],
                linestyle=line_styles[index % len(line_styles)],
                linewidth=1.5 + 0.15 * (index // len(line_styles)),
                alpha=0.55 + 0.45 * (index + 1) / radii.fraction.size,
                label=rf"$M(r)/M(\infty)={fraction:g}$",
                require_positive_x=False,
            )
        ax.set_yscale("log")
        ax.set_ylabel(r"$\hat{r}$")
        ax.set_title(COMPONENT_LABELS[component])
        ax.legend(loc="best", frameon=False, ncol=2)
    axes[-1].set_xlabel(r"$\hat{t}$")
    fig.suptitle(f"{title}: Lagrangian radii", fontsize=18)
    fig.tight_layout()
    return fig


def _central_figure(central: CentralData, quantity: str, title: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8.0, 8.0))
    for component in COMPONENTS:
        values = (
            central.rho[component]
            if quantity == "density"
            else np.sqrt((2.0 / 3.0) * central.energy[component])
        )
        style = COMPONENT_STYLES[component]
        _plot_positive(
            ax,
            central.time,
            values,
            color=style["color"],
            linestyle=style["linestyle"],
            linewidth=1.8,
            label=COMPONENT_LABELS[component],
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\hat{t}$")
    ax.set_ylabel(
        r"$\hat{\rho}(0)$" if quantity == "density" else r"$\hat{\sigma}(0)$"
    )
    ax.set_title(title)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="best", frameon=False)
    fig.tight_layout()
    return fig


def _save_figure(
    figure: plt.Figure,
    plot_directory: Path,
    stem: str,
    formats: tuple[str, ...],
    dpi: int,
) -> list[Path]:
    paths: list[Path] = []
    for output_format in formats:
        path = plot_directory / f"{stem}.{output_format}"
        save_options: dict[str, object] = {"bbox_inches": "tight"}
        if output_format == "png":
            save_options["dpi"] = dpi
        figure.savefig(path, **save_options)
        paths.append(path)
    plt.close(figure)
    return paths


def generate_plots(
    run: RunData,
    plot_directory: str | Path | None = None,
    formats: tuple[str, ...] = ("png", "pdf"),
    title: str | None = None,
    maximum_snapshots: int = 8,
    dpi: int = 160,
) -> list[Path]:
    """Generate every plot supported by the observer data present in ``run``."""
    if not formats or any(value not in {"png", "pdf", "svg"} for value in formats):
        raise ValueError("formats must be a nonempty selection of png, pdf, or svg")
    if dpi <= 0:
        raise ValueError("dpi must be positive")
    plot_directory = (
        Path(plot_directory).expanduser().resolve()
        if plot_directory is not None
        else run.directory / "plots"
    )
    plot_directory.mkdir(parents=True, exist_ok=True)
    title = title or run.directory.name.replace("_", " ")
    configure_style()

    paths: list[Path] = []
    if run.profiles is not None:
        paths.extend(
            _save_figure(
                _profile_figure(run.profiles, "density", title, maximum_snapshots),
                plot_directory,
                "density_profiles",
                formats,
                dpi,
            )
        )
        paths.extend(
            _save_figure(
                _profile_figure(run.profiles, "sigma", title, maximum_snapshots),
                plot_directory,
                "velocity_dispersion_profiles",
                formats,
                dpi,
            )
        )
    if run.radii is not None:
        paths.extend(
            _save_figure(
                _radii_figure(run.radii, title),
                plot_directory,
                "lagrangian_radii",
                formats,
                dpi,
            )
        )
    if run.central is not None:
        paths.extend(
            _save_figure(
                _central_figure(run.central, "density", f"{title}: central density"),
                plot_directory,
                "central_density",
                formats,
                dpi,
            )
        )
        paths.extend(
            _save_figure(
                _central_figure(
                    run.central,
                    "sigma",
                    f"{title}: central velocity dispersion",
                ),
                plot_directory,
                "central_velocity_dispersion",
                formats,
                dpi,
            )
        )
    return paths


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot one cluster-hydro-sim packed output directory."
    )
    parser.add_argument("output_directory", type=Path)
    parser.add_argument(
        "--plot-directory",
        type=Path,
        help="destination directory (default: OUTPUT_DIRECTORY/plots)",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "pdf", "svg"),
        default=("png", "pdf"),
        help="one or more output formats (default: png pdf)",
    )
    parser.add_argument("--title", help="figure title (default: output directory name)")
    parser.add_argument(
        "--max-profile-snapshots",
        type=int,
        default=8,
        help="maximum profile snapshots to draw, sampled across the run (default: 8)",
    )
    parser.add_argument("--dpi", type=int, default=160, help="PNG resolution (default: 160)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_argument_parser().parse_args(argv)
    try:
        run = load_output(args.output_directory)
        paths = generate_plots(
            run,
            plot_directory=args.plot_directory,
            formats=tuple(dict.fromkeys(args.formats)),
            title=args.title,
            maximum_snapshots=args.max_profile_snapshots,
            dpi=args.dpi,
        )
    except (OutputFormatError, OSError, ValueError) as exc:
        print(f"plot_output.py: error: {exc}", file=sys.stderr)
        return 2
    print(f"Loaded {run.directory}")
    for path in paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
