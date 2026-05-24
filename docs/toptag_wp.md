# Deriving GloParTv3 top-tagging working points

This documents the standalone framework for deriving our own 2024 top-tag
working points (the `_TAGGER_WPS['2024']` values in `ttbarprocessor.py` are
placeholders inherited from the Run-2 `tau3/tau2` convention and do not apply to
the GloParTv3 `TopvsQCD` score the analysis actually uses).

The design and rationale live in `toptag_wp_derivation_plan.md`. This file is
the operational how-to.

## What it does

- **WP threshold** for each target QCD mis-tag rate
  (`very_tight 0.1%`, `tight 0.5%`, `medium 1.0%`, `loose 2.5%`, `very_loose 5.0%`),
  derived **per AK8 pT bin**, inside the `105 < mSD < 210` mass window,
  `pT > 400`, `|eta| < 2.5`.
- **Signal efficiency** at each threshold, from AK8 jets gen-matched to a
  hadronic top (`dR < 0.8`) in TTbar.
- Discriminant is `(TopbWqq + TopbWq) / (TopbWqq + TopbWq + QCD)` — identical to
  `TTbarResProcessor._tscore`.
- The WP threshold is defined purely by the **QCD background**; the signal sample
  only supplies the efficiency reported alongside.
- QCD events with large `genWeight` outliers are rejected with the same
  two-standard-deviation filter used in `TTbarResProcessor`; the runner prints
  the kept/raw event counts and the number rejected per dataset.

No ntuples are produced — histograms only.

## Files

| File | Role |
|---|---|
| `python/toptag_wp_processor.py` | `TopTagWPProcessor` — fills per-jet tagger histograms. Year-parametrized via `TAGGER_CONFIG` (no hard-coded year). |
| `run_toptag_wp.py` | Builds a local or manifest-backed NanoAOD v15 fileset, runs locally/LPC/casa, saves a `.coffea`. |
| `plot_toptag_wp.py` | Reads the `.coffea`, derives WPs, writes JSON + plots. |

## Requirements

- Python env with the analysis stack (coffea ≥ 2026, hist, uproot, mplhep).
  Local laptop examples use `coffea-dask/bin/python`; LPC and coffea.casa
  provide the environment, so use `python` there.
- **NanoAOD v15** inputs — the `globalParT3_*` branches only exist from v15.
- Local input layout (default `--rootdir ~/Projects/rootfiles/ttbar`):
  ```
  <rootdir>/<iov>/mc/<SAMPLE>/*.root
  ```
  Directory names starting with `TT` are treated as signal (gen-matched tops);
  everything else is background (QCD).
- LPC and coffea.casa MC runs read `data/nanoAOD/QCD.json` and
  `data/nanoAOD/TTbar.json` by default. Data-only runs use
  `data/nanoAOD/data.json`. LPC uses
  `root://cmsxrootd.fnal.gov/`; coffea.casa uses `root://xcache/`.

## Run it

1. **Local smoke test on a laptop/Mac** (1 local file per dataset):
   ```bash
   coffea-dask/bin/python run_toptag_wp.py \
       --env local --test --maxfiles 1 \
       --rootdir /Users/aritra/Projects/rootfiles/ttbar \
       --out outputs/toptag_wp_2024_local_smoke.coffea
   ```

2. **Full local run over local ROOT files:**
   ```bash
   coffea-dask/bin/python run_toptag_wp.py \
       --env local --workers 4 \
       --rootdir /Users/aritra/Projects/rootfiles/ttbar \
       --out outputs/toptag_wp_2024_local.coffea
   ```

3. **Full LPC run over the 2024 manifests:**
   ```bash
   python run_toptag_wp.py \
       --env lpc \
       --out outputs/toptag_wp_2024_full.coffea
   ```

4. **coffea.casa run over the 2024 manifests:**
   ```bash
   python run_toptag_wp.py \
       --env casa \
       --out outputs/toptag_wp_2024_full.coffea
   ```

5. **Remote smoke test** (2 files per manifest dataset, 1 chunk per dataset):
   ```bash
   python run_toptag_wp.py \
       --env lpc --test \
       --out outputs/toptag_wp_2024_lpc_smoke.coffea
   ```

6. **Data-only run for the TopvsQCD data/MC comparison:**
   ```bash
   python run_toptag_wp.py \
       --env lpc --sample Data \
       --out outputs/toptag_score_data_2024.coffea
   ```
   Use `--test` here as well for a two-files-per-era smoke test.

7. **Derive WPs and make plots, including Data/MC score shapes if data exists:**
   ```bash
   MPLCONFIGDIR=/tmp/mplconfig python plot_toptag_wp.py \
       outputs/toptag_wp_2024_full.coffea \
       --data-infile outputs/toptag_score_data_2024.coffea \
       --iov 2024 \
       --json data/toptag/toptag_wp_2024.json \
       --plotdir plots/images/toptag_wp/2024 \
       --score-rebin 10
   ```
   Writes `data/toptag/toptag_wp_2024.json` and plots under
   `plots/images/toptag_wp/2024/`, and prints a threshold/signal-eff table.
   `--score-rebin` is display-only for `score_distributions.png`; the WP
   derivation still uses the fine score histogram.

## Outputs

- `data/toptag/toptag_wp_<iov>.json` — pT-binned thresholds, threshold
  uncertainty (MC-stat), achieved mis-tag, and signal efficiency per WP, plus
  provenance.
- Plots in `plots/images/toptag_wp/<iov>/`:
  - `score_distributions.png` — signal vs background `TopvsQCD`, per pT bin.
  - `data_mc_score_distributions.png` — data points over stacked QCD+TTbar
    `TopvsQCD`, per pT bin, when `--data-infile` is supplied or data is present
    in the input accumulator. The full stacked MC is scaled to the data integral
    in each pT panel for a shape comparison, while the nominal QCD/TTbar relative
    composition is preserved. The hatched band shows the MC statistical
    uncertainty, with a `Data/MC` ratio panel below each score distribution.
  - `roc.png` — signal eff vs mis-tag, log-x, target points marked.
  - `wp_threshold_vs_pt.png` — derived thresholds vs pT (the parametrization).
  - `signal_eff_vs_pt.png` — signal eff at each WP vs pT.
  - `mistag_closure_vs_pt.png` — achieved mis-tag (should sit on each target).
  - `mistag_vs_msd.png` — mis-tag vs mSD at the medium WP (decorrelation check).

These values are **stored, not wired into the analysis** — integrating them into
`_TAGGER_WPS` is a deliberate later step.

## Adding another year

Once NanoAOD v15 exists for another IOV:

1. Add an entry to `TAGGER_CONFIG` in `python/toptag_wp_processor.py` (the
   discriminant is the same `globalParT3` formula; just list it under the new
   IOV key, with `required_fields`).
2. Add that year's QCD/TTbar manifests under `data/nanoAOD/`, or drop local
   files under `<rootdir>/<iov>/mc/<SAMPLE>/` for laptop tests.
3. Add the year's luminosity to `LUMI_PB` in `run_toptag_wp.py` and make sure
   manifest metadata contains per-sample `xsec_pb`.
4. Rerun the three commands with `--iov <year>`.

## Known limitations

- The local Mac sample is only a smoke-test subset. Use `--env lpc` or
  `--env casa` for the full 2024 QCD mixture and TTbar statistics.
- **Scale factors are not produced here.** Top-tag data/MC SFs require a
  separate semileptonic-ttbar tag-and-probe measurement; see §6 of
  `toptag_wp_derivation_plan.md`.
