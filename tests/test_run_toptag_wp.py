import json
import tempfile
import unittest
from pathlib import Path

import run_toptag_wp


class TopTagRunnerFilesetTest(unittest.TestCase):
    def test_local_fileset_uses_rootdir_and_maxfiles(self):
        with tempfile.TemporaryDirectory() as tmp:
            rootdir = Path(tmp)
            sample_dir = rootdir / "2024" / "mc" / "TTto4Q"
            sample_dir.mkdir(parents=True)
            for idx in range(3):
                (sample_dir / f"file{idx}.root").touch()

            fileset = run_toptag_wp.build_local_fileset(rootdir, "2024", maxfiles=2)

        self.assertEqual(list(fileset), ["TTto4Q"])
        self.assertEqual(len(fileset["TTto4Q"]["files"]), 2)
        self.assertEqual(fileset["TTto4Q"]["metadata"]["xsec_pb"], 350.6)

    def test_manifest_fileset_uses_redirector_metadata_and_maxfiles(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            qcd_json = tmpdir / "QCD.json"
            ttbar_json = tmpdir / "TTbar.json"
            qcd_json.write_text(json.dumps({
                "2024": {
                    "QCD_PT-600to800": {
                        "files": ["/store/qcd/a.root", "/store/qcd/b.root", "/store/qcd/c.root"],
                        "metadata": {
                            "sample": "QCD",
                            "subsample": "QCD_PT-600to800",
                            "year": "2024",
                            "is_mc": True,
                            "xsec_pb": 178.7,
                        },
                    }
                }
            }))
            ttbar_json.write_text(json.dumps({
                "2024": {
                    "inclusive": {
                        "files": ["/store/tt/a.root", "/store/tt/b.root"],
                        "metadata": {
                            "sample": "TTbar",
                            "subsample": "inclusive",
                            "year": "2024",
                            "is_mc": True,
                            "xsec_pb": 350.6,
                        },
                    }
                }
            }))

            fileset = run_toptag_wp.build_manifest_fileset(
                qcd_json,
                ttbar_json,
                "2024",
                redirector="root://cmsxrootd.fnal.gov/",
                maxfiles=1,
            )

        self.assertEqual(set(fileset), {"QCD_PT-600to800", "TTbar"})
        self.assertEqual(fileset["QCD_PT-600to800"]["files"], [
            "root://cmsxrootd.fnal.gov//store/qcd/a.root"
        ])
        self.assertEqual(fileset["TTbar"]["metadata"]["xsec_pb"], 350.6)


if __name__ == "__main__":
    unittest.main()
