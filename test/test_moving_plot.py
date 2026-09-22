import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec = importlib.util.spec_from_file_location("moving_plot", Path(__file__).parents[1]/"script/plot_moving_comparison.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class MovingPlots(unittest.TestCase):
    def test_matched_physical_time_radius(self):
        a = np.ones((3, 3, 5));a[:, :, 0] = np.array([1, 2, 4])[:, None]
        a[:, :, 1] = np.array([1, .25, .0625])[:, None]
        b = a.copy();b[:, :, 1] *= 3
        result = module.matched_profile(np.array([0, 2]), [a, b], 1, np.array([np.sqrt(2)]), 0, 1)
        np.testing.assert_allclose(result, [1])
    def test_no_extrapolation(self):
        a = np.ones((2, 3, 5));a[1, :, 0] = 2
        for time, radii in ((3, [1]), (1, [.5]), (1, [3])):
            with self.assertRaises(ValueError):
                module.matched_profile(np.array([0, 2]), [a, a], time, np.array(radii), 0, 1)

if __name__ == "__main__":unittest.main()
