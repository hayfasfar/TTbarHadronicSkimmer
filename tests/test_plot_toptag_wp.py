import unittest

import numpy as np

import plot_toptag_wp


class TopTagPlotClassificationTest(unittest.TestCase):
    def test_data_is_not_wp_background(self):
        meta = {
            "QCD_PT-600to800": {"sample": "QCD", "is_mc": True},
            "TTbar": {"sample": "TTbar", "is_mc": True},
            "Data_C": {"sample": "Data", "is_mc": False},
        }

        sig, bkg = plot_toptag_wp.classify_datasets(meta)

        self.assertEqual(sig, ["TTbar"])
        self.assertEqual(bkg, ["QCD_PT-600to800"])
        self.assertEqual(plot_toptag_wp.classify_data_datasets(meta), ["Data_C"])

    def test_mc_shape_scale_matches_data_integral(self):
        scale = plot_toptag_wp.mc_shape_scale_to_data(
            data_total=100.0,
            mc_total=80.0,
        )

        self.assertAlmostEqual(scale, 1.25)
        self.assertEqual(
            plot_toptag_wp.mc_shape_scale_to_data(10.0, 0.0),
            0.0,
        )

    def test_data_mc_ratio_skips_empty_mc_bins(self):
        ratio, ratio_err = plot_toptag_wp.data_mc_ratio(
            data_counts=np.array([10.0, 4.0, 2.0]),
            data_variance=np.array([10.0, 4.0, 2.0]),
            mc_counts=np.array([5.0, 0.0, 4.0]),
        )

        self.assertAlmostEqual(ratio[0], 2.0)
        self.assertAlmostEqual(ratio_err[0], np.sqrt(10.0) / 5.0)
        self.assertTrue(np.isnan(ratio[1]))
        self.assertTrue(np.isnan(ratio_err[1]))
        self.assertAlmostEqual(ratio[2], 0.5)


if __name__ == "__main__":
    unittest.main()
