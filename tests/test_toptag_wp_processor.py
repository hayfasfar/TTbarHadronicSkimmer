import os
import sys
import unittest

import awkward as ak

sys.path.append(os.path.join(os.getcwd(), "python"))
import toptag_wp_processor  # noqa: E402


class TopTagWPProcessorWeightFilterTest(unittest.TestCase):
    def test_qcd_genweight_mask_rejects_large_outlier(self):
        mask = toptag_wp_processor._qcd_genweight_mask(
            ak.Array([1.0, 1.0, 1.0, 1.0, 1.0, 1000.0])
        )

        self.assertEqual(ak.to_list(mask), [True, True, True, True, True, False])

    def test_qcd_genweight_mask_keeps_flat_weights(self):
        mask = toptag_wp_processor._qcd_genweight_mask(ak.Array([1.0, 1.0, 1.0]))

        self.assertEqual(ak.to_list(mask), [True, True, True])

    def test_jet_pt_hist_uses_fine_control_bins(self):
        hpt = toptag_wp_processor._make_jet_pt_hist()

        self.assertIn("pt", [axis.name for axis in hpt.axes])
        self.assertGreater(hpt.axes["pt"].size, len(toptag_wp_processor.PT_BIN_EDGES) - 1)


if __name__ == "__main__":
    unittest.main()
