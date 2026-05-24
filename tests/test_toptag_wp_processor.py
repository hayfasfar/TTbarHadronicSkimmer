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


if __name__ == "__main__":
    unittest.main()
