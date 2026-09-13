"""Synthetic AB overlay checks; no reference paper or simulation run required."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "script"))
import plot_ab_comparison as plot


def parameters(path, fields):
    path.mkdir(parents=True, exist_ok=True)
    (path / "paramNames.txt").write_text("\n".join(fields)+"\n")
    (path / "paramTypes.txt").write_text("Real64\n"*len(fields))
    (path / "paramOffsets.txt").write_text("\n".join(str(8*i) for i in range(len(fields)))+"\n")
    np.array(list(fields.values()), dtype="<f8").tofile(path / "param.dat")


class ABPlotTest(unittest.TestCase):
    def test_units_frames_and_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)
            parameters(path, {"N": 3, "tidal_q": .001})
            parameters(path/"observer", {"time_unit_myr": 20, "friction_time_myr": 10, "radius_to_nbody": .5})
            np.array([0,1,2],dtype="<f8").tofile(path/"radii_t.dat")
            plot.FRACTIONS.astype("<f8").tofile(path/"radii_fraction.dat")
            for f in ("s","b","d"):
                np.ones((3,6),dtype="<f8").tofile(path/f"radii_{f}.dat")
            time,radii,params=plot.load_ab(path)
            np.testing.assert_array_equal(time,[0,2,4])
            np.testing.assert_array_equal(radii["s"],np.full((3,6),.5))
            self.assertEqual(params["tidal_q"],.001)
            with plot.fitz.open() as document:
                for _ in range(3): document.new_page(width=595,height=842)
                fig=plot.make_comparison(document,path,path)
                try:
                    fig.canvas.draw()
                    for axis in fig.axes:
                        self.assertEqual(len(axis.lines),12)
                        self.assertEqual(axis.get_yscale(),"log")
                        self.assertEqual(axis.get_xlim(),(0,4))
                        self.assertEqual(axis.get_ylim(),(.1,4))
                        h,w=axis.images[0].get_array().shape[:2]
                        box=axis.get_window_extent()
                        self.assertAlmostEqual(box.height/box.width,h/w,places=12)
                    fig.savefig(path/"overlay.pdf")
                finally: plot.plt.close(fig)
            with plot.fitz.open(path/"overlay.pdf") as document:
                self.assertEqual(len(document),1)
            parameters(path/"observer", {"time_unit_myr": -20, "friction_time_myr": 10, "radius_to_nbody": .5})
            with self.assertRaises(ValueError): plot.load_ab(path)


if __name__=="__main__": unittest.main()
