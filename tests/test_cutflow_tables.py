import os
import sys
import unittest

sys.path.append(os.path.join(os.getcwd(), "python"))

from cutflow import build_cutflow_rows, format_latex_table, format_markdown_table  # noqa: E402


class CutflowTableFormattingTest(unittest.TestCase):
    def test_build_rows_uses_scaled_weights_and_efficiencies(self):
        output = {
            "cutflow_table_steps": [
                {"key": "input_events", "label": "Input events", "group": "Bookkeeping"},
                {"key": "trigger", "label": "Trigger", "group": "Preselection"},
                {"key": "tag_2tag", "label": "Two top-tagged jets", "group": "Tagging"},
            ],
            "cutflow_unweighted": {
                "input_events": 100,
                "trigger": 80,
                "tag_2tag": 20,
            },
            "cutflow_weighted": {
                "input_events": 50.0,
                "trigger": 40.0,
                "tag_2tag": 10.0,
            },
            "cutflow_weighted_scaled": {
                "input_events": 500.0,
                "trigger": 400.0,
                "tag_2tag": 100.0,
            },
        }

        rows = build_cutflow_rows(output)

        self.assertEqual([row["key"] for row in rows], ["input_events", "trigger", "tag_2tag"])
        self.assertEqual(rows[1]["weighted"], 400.0)
        self.assertAlmostEqual(rows[1]["eff_prev_percent"], 80.0)
        self.assertAlmostEqual(rows[2]["eff_total_percent"], 20.0)

    def test_markdown_and_latex_render_core_columns(self):
        rows = [
            {
                "group": "Bookkeeping",
                "step": "Input events",
                "events": 100.0,
                "weighted": 500.0,
                "eff_prev_percent": None,
                "eff_total_percent": 100.0,
                "weight_label": "Yield",
            }
        ]

        markdown = format_markdown_table(rows, title="Cutflow")
        latex = format_latex_table(rows, caption="Cutflow")

        self.assertIn("| Group | Step | Events | Yield |", markdown)
        self.assertIn("Input events", markdown)
        self.assertIn(r"\begin{table}", latex)
        self.assertIn("Input events", latex)

    def test_legacy_cutflow_keys_are_adapted(self):
        output = {
            "cutflow": {
                "all events 1": 100,
                "all events": 90,
                "trigger": 80,
                "after_eventCut": 40,
                "after_ttbarcandCuts": 20,
                "2t0bcen": 5,
            },
            "cutflow_scaled": {
                "all events 1": 1000.0,
                "all events": 900.0,
                "trigger": 800.0,
                "after_eventCut": 400.0,
                "after_ttbarcandCuts": 200.0,
                "2t0bcen": 50.0,
            },
        }

        rows = build_cutflow_rows(output)

        self.assertIn("preselection", [row["key"] for row in rows])
        self.assertIn("ttbarcand", [row["key"] for row in rows])
        self.assertIn("category_2t0bcen", [row["key"] for row in rows])
        self.assertEqual(rows[-1]["weighted"], 50.0)
        self.assertAlmostEqual(rows[-1]["eff_prev_percent"], 25.0)


if __name__ == "__main__":
    unittest.main()
