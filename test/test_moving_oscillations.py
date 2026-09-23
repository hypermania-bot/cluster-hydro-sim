from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).parents[1]/"script"))
from plot_moving_oscillation_check import oscillation_metrics, frame_at


class OscillationDiagnostics(unittest.TestCase):
    def test_checkerboard_and_smooth_profile(self):
        x = np.ones((20, 3, 5));x[:, :, 0] = np.linspace(.11, .99, 20)[:, None]
        x[:, 0, 4] = np.where(np.arange(20) % 2, 2., -2.)
        self.assertEqual(oscillation_metrics(x, 0)["median"], 2)
        self.assertEqual(oscillation_metrics(x, 0)["sign_changes"], 19)
        x[:, 0, 4] = np.arange(20)
        self.assertEqual(oscillation_metrics(x, 0)["maximum"], 0)

    def test_signed_velocity_interpolation(self):
        x = np.ones((3, 3, 5));y = x.copy();x[:, :, 4] = -2;y[:, :, 4] = 2
        np.testing.assert_array_equal(frame_at(np.array([0., 1.]), [x, y], .5)[:, :, 4], 0)
        with self.assertRaises(ValueError):frame_at(np.array([0., 1.]), [x, y], 2)


if __name__ == "__main__":unittest.main()
