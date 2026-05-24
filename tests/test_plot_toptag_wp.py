import unittest

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


if __name__ == "__main__":
    unittest.main()
