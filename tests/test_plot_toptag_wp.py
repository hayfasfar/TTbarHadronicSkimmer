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

    def test_qcd_scale_matches_data_minus_ttbar(self):
        scale = plot_toptag_wp.qcd_scale_to_data_minus_ttbar(
            data_total=100.0,
            ttbar_total=15.0,
            qcd_total=50.0,
        )

        self.assertAlmostEqual(scale, 1.7)
        self.assertEqual(
            plot_toptag_wp.qcd_scale_to_data_minus_ttbar(10.0, 15.0, 50.0),
            0.0,
        )
        self.assertEqual(
            plot_toptag_wp.qcd_scale_to_data_minus_ttbar(10.0, 1.0, 0.0),
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
