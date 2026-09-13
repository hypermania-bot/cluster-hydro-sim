#!/usr/bin/env python3
"""Overlay saved HydroSim Lagrangian radii on Ardi & Baumgardt (2020), Fig. 1."""
from pathlib import Path
import argparse

import fitz
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from plot_output import load_output, read_parameters


# Plot frames on PDF page 3 (595 x 842 points), excluding axis labels.
FRAMES = [(111.3, 414.7, 283.0, 537.3), (332.0, 414.7, 504.0, 537.3)]
FRACTIONS = np.array([.01, .05, .1, .2, .5, .7])


def load_ab(directory):
    directory = Path(directory)
    run = load_output(directory)
    units = read_parameters(directory / "observer")
    if run.radii is None or not np.allclose(run.radii.fraction, FRACTIONS):
        raise ValueError("AB comparison requires the six published mass fractions")
    time_scale = units["time_unit_myr"] / units["friction_time_myr"]
    length_scale = units["radius_to_nbody"]
    if not np.isfinite(time_scale*length_scale) or min(time_scale,length_scale)<=0:
        raise ValueError("invalid AB observer units")
    return run.radii.time*time_scale, {f: run.radii.radius[f]*length_scale for f in ("s", "d")}, run.parameters


def make_comparison(reference, isolated, tidal):
    fig, axes = plt.subplots(1, 2, figsize=(11.69, 5.8))
    fig.subplots_adjust(left=.075, right=.98, top=.85, bottom=.27, wspace=.24)
    for axis, directory, frame, title in zip(axes, (isolated, tidal), FRAMES,
                                           ("Isolated", "Spherical tidal surrogate"), strict=True):
        time, radii, params = load_ab(directory)
        page = reference[2]
        rectangle = fitz.Rect(*(x * (page.rect.width/595 if i%2==0 else page.rect.height/842)
                               for i,x in enumerate(frame)))
        pix = page.get_pixmap(matrix=fitz.Matrix(3,3), clip=rectangle, alpha=False)
        scan = np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,pix.n)[...,:3]
        axis.imshow(scan, extent=(0,1,0,1), transform=axis.transAxes,
                    origin="upper", aspect="auto", alpha=.35, zorder=0)
        axis.set_box_aspect(scan.shape[0]/scan.shape[1])
        for component,color,label in [("s","#d62728","Stars"),("d","#1f77b4","DM")]:
            for i in range(len(FRACTIONS)):
                axis.plot(time, radii[component][:,i], color=color, linestyle="--",
                          linewidth=1.0, label=label if i==0 else None)
        axis.set(xlim=(0,4),ylim=(.1,4),yscale="log",xlabel=r"$t/t_{\rm fric}$",
                 ylabel="Lagrangian radius (N-body units)",title=f"{title}: q={params['tidal_q']:.3g}")
        axis.legend(loc="upper left",fontsize=8,framealpha=.9)
    fig.suptitle("Ardi & Baumgardt Fig. 1: Lagrangian radii",fontsize=14)
    fig.text(.075,.15,"Dashed: HydroSim. Faded scan: AB. Fractions: 1, 5, 10, 20, 50, 70% of each component's current mass.",fontsize=9)
    fig.text(.075,.09,"AB does not specify the galactic field strength or the physical scale of Fig. 1.\n"
             "HydroSim uses a stated flat-rotation surrogate, a finite hydrostatic domain, and one stellar mass; no escape law.",fontsize=8)
    fig.text(.075,.025,"Reference: E. Ardi & H. Baumgardt, J. Phys.: Conf. Ser. 1503 (2020) 012023, Fig. 1. DOI: 10.1088/1742-6596/1503/1/012023. CC BY 3.0.",fontsize=7)
    return fig


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_pdf",type=Path)
    parser.add_argument("isolated_directory",type=Path)
    parser.add_argument("tidal_directory",type=Path)
    parser.add_argument("--output",type=Path,default=Path("output/ab_tide/ab_fig1_overlay.pdf"))
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with fitz.open(args.reference_pdf) as document:
        fig=make_comparison(document,args.isolated_directory,args.tidal_directory)
        fig.savefig(args.output);plt.close(fig)
    print(args.output)


if __name__=="__main__":
    main()
