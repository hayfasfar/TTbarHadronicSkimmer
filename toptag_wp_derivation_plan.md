# Plan: Deriving In-House GloParTv3 Top-Tagging Working Points (2024)

## 1. Problem statement

The 2024 top-tag working points currently hard-coded in `ttbarprocessor.py`
(`_TAGGER_WPS`, lines 99–118) are **placeholders**. They were inherited from
the Run-2 / `tau3/tau2` mistag convention and do not correspond to any measured
working point for the tagger we actually use in 2024.

Two things are conflated today:

1. **The official top-tag definition** (the text the user quoted): AK8 PUPPI(v15)
   jets in `105 < mSD < 210 GeV` with a cut on the N-subjettiness ratio
   `tau3/tau2`. The five WPs (Very tight 0.38 → Very loose 0.69) are defined by a
   **target QCD mis-tag efficiency** (0.1%, 0.5%, 1.0%, 2.5%, 5.0%).

2. **What our analysis actually uses**: the mass-decorrelated GloParTv3
   `TopvsQCD` score,
   `(TopbWqq + TopbWq) / (TopbWqq + TopbWq + QCD)`, computed in
   `TTbarResProcessor._tscore` (`ttbarprocessor.py:272`). This is a different
   discriminant from `tau3/tau2`, so the `tau3/tau2` thresholds **do not map onto
   it at all**, and official GloParTv3 2024 top WPs are not yet available to us.

**Goal:** build a small, reproducible framework that derives our own GloParTv3
`TopvsQCD` working points for 2024 using the *same convention* as the official
table — fix a target QCD mis-tag efficiency, find the score threshold that
achieves it (**as a function of jet pT**), and report the corresponding
top-tagging (signal) efficiency. Then emit the validation plots, and lay out how
the matching **data/MC scale factors** are measured.

This document is a **plan only** — no code is written yet.

## 2. Locked-in decisions (from user)

| Item | Decision |
|---|---|
| Discriminant `D` | GloParTv3 `TopvsQCD` exactly as in `_tscore` (`ttbarprocessor.py:272`). |
| Jet pre-selection | `pT > 400 GeV`, `|eta| < 2.5`. |
| Mass window | **Apply** `105 < mSD < 210 GeV` when deriving (same as the analysis signal region). |
| pT dependence | **pT-dependent** WPs (binned thresholds), not a single scalar. |
| Background | **QCD MC only** (2024). |
| Signal | TTbar MC 2024, AK8 matched to a gen hadronic top (`get_hadronic_tops`, `dR < 0.8`). |
| Scale factors | **Yes, needed** — see §6 for the derivation method. |
| Implementation | **New standalone processor**, 2024 only for now, but **year-parametrized** (no hard-coded year) so other years drop in once NanoAOD v15 exists for them. |
| Integration into analysis | **Not now.** Just produce and store the WP values (and SFs) in a file. |

**Definitions:**

| Quantity | Definition |
|---|---|
| Signal jet | AK8 jet in TTbar MC matched within `dR < 0.8` to a gen-level hadronic top, passing pre-selection + mass window. |
| Background jet | AK8 jet in QCD MC passing pre-selection + mass window. |
| Mis-tag eff `ε_bkg(t)` | weighted fraction of background jets with `D > t`. |
| Signal eff `ε_sig(t)` | weighted fraction of signal jets with `D > t`. |
| Working point `t*` | threshold where `ε_bkg(t*)` = target (0.1 / 0.5 / 1.0 / 2.5 / 5.0 %), **per pT bin**. |

## 3. Inputs / samples

- **TTbar MC 2024** (signal jets) — must be produced; no coffea output exists yet.
- **QCD multijet MC 2024**, HT-binned (`QCD-4Jets_HT-*`; see wiki *ttbarhadronic
  Run-3 QCD samples*) — mis-tag denominator. **HT bins must be stitched with
  cross-section × lumi weights**, or the inclusive mis-tag rate is meaningless.
- Both available as filesets under `data/nanoAOD/*.json`.
- **NanoAOD v15 is required** — the `FatJet_globalParT3_*` branches only exist
  from v15 onward. This is also why only 2024 is in scope right now.

Per-jet weight = generator weight × MC normalization. Carry it through; the WP
inversion and all efficiencies are weighted.

## 4. Framework architecture — new standalone processor

A dedicated, **year-parametrized** processor (e.g.
`python/toptag_wp_processor.py`, class `TopTagWPProcessor`) separate from
`TTbarResProcessor`. It does **not** run the full event selection — it dumps
per-AK8-jet records *before* any top-tag cut, since that is exactly the
population we need to scan.

Year handling — no hard-coded `'2024'`. A module-level config keyed by IOV:

```python
TAGGER_CONFIG = {
    '2024': {
        'score': lambda fj: (fj.globalParT3_TopbWqq + fj.globalParT3_TopbWq)
                            / (fj.globalParT3_TopbWqq + fj.globalParT3_TopbWq + fj.globalParT3_QCD),
        'required_fields': ['globalParT3_TopbWqq', 'globalParT3_TopbWq', 'globalParT3_QCD'],
    },
    # future: '2023', '2025', ... add an entry once NanoAOD v15 is available
}
```

The processor reads `self.iov`, looks up `TAGGER_CONFIG[self.iov]`, and raises a
clear error if the IOV is absent or the required fields are missing from
`events.FatJet.fields`. Reuse the `_tscore` formula verbatim so the derived WP
matches what the analysis applies.

**Per-jet output record:**
`pt, eta, msd, globalParT3_TopbWqq, globalParT3_TopbWq, globalParT3_QCD, D,
is_genmatched_top, weight, sample, iov`.

Storing the **raw score components** (not just `D`) means alternative
discriminants can be re-derived later without re-running the processor.

**Selection inside the processor:**
- Pre-selection mask: `pt > 400`, `|eta| < 2.5`, `105 < msd < 210`.
- Signal flag: for TTbar samples, match each AK8 to gen hadronic tops
  (`get_hadronic_tops` + `dR < 0.8`, reusing `python/truthstudy.py` logic).
  For QCD samples, `is_genmatched_top = False` for all.
- Output: a flat table (parquet, or a coffea `column_accumulator` set) of all
  surviving jets. This is small (one row per jet, a handful of floats), so no
  chunked-ROOT machinery is needed — accumulate and write directly.

Analysis of the table (WP inversion + plots + SFs) lives in a standalone
notebook, `notebooks/toptag_wp.ipynb`, not in the processor.

## 5. WP derivation procedure

pT bins (proposed): `400–500, 500–600, 600–800, 800–1200, >1200 GeV`.

Per pT bin:

1. Build the weighted background `D` distribution; compute the tail
   `ε_bkg(t)` from high `D` downward (cumulative).
2. Invert to find `t*` for each target mis-tag (0.1, 0.5, 1.0, 2.5, 5.0 %) by
   interpolation.
3. Evaluate `ε_sig(t*)` on the signal table in the same pT bin.
4. Estimate `t*` uncertainty from MC stats (binomial/bootstrap on the tail
   count) — important for the high-pT bins where QCD stats thin out.

Name the five WPs to mirror the official table: Very tight / Tight / Medium /
Loose / Very loose ↔ 0.1 / 0.5 / 1.0 / 2.5 / 5.0 %.

## 6. Scale-factor derivation (details requested)

A working point only tells you *where to cut*. A **scale factor** corrects the
fact that the top-tag efficiency in **data** differs from the efficiency in
**MC**:

```
SF(pT) = ε_top-tag^data(pT) / ε_top-tag^MC(pT)
```

measured per WP and per AK8 pT bin, with up/down systematic variations. The
analysis then weights each tagged top in MC by `SF` (and the SF uncertainty
becomes a systematic — exactly the `ttag_pt*` systematics already scaffolded in
`python/weights.py`).

**The core difficulty:** to measure `ε_data` you need a *pure, known* sample of
real hadronic tops in data. Our signal region is fully hadronic and
QCD-dominated, so it cannot provide that. The standard CMS solution is a
**semileptonic ttbar tag-and-probe**, which is a *separate selection* from the
all-hadronic analysis:

1. **Tag side (selects a clean ttbar sample):** require one isolated
   high-pT muon (or electron) + MET + a b-tagged AK4 jet. This "tags" the event
   as a leptonically-decaying top and makes the *other* side very likely a
   genuine hadronically-decaying top.
2. **Probe side:** the boosted AK8 jet recoiling against the lepton (large
   `dPhi`, `pT > 400`, in the mass window). This is the probe top.
3. **Measure efficiency:** `ε = N(probe passes TopvsQCD WP) / N(all probes)`,
   in pT bins, separately in data and MC.
4. **Background subtraction / fit:** the probe sample is not 100% real tops
   (W+jets, single-top, QCD, and partially-merged tops contaminate it). So:
   - Subtract non-ttbar backgrounds in data using their MC predictions.
   - Split signal MC by **gen-merging category** via gen-matching:
     *fully-merged* top (all 3 daughter quarks inside the AK8 cone),
     *partially-merged* (e.g. only the W, or 2 of 3 quarks), and *not-merged /
     other*. Each category has a very different tag efficiency.
   - Fit the soft-drop mass distribution (pass vs fail the WP) in data to
     extract the data efficiency of the **merged-top** component; the
     not-merged component typically gets its own SF or is constrained by MC.
5. **Result:** `SF(pT) = ε_data / ε_MC` for the merged-top component, per WP,
   per pT bin, with systematics from the background normalization, the
   mass-fit model, the merged/not-merged split, and JES/JER on the probe jet.

**Mis-tag (background) SF.** Because this analysis estimates QCD from data via
2DAlphabet, the QCD mis-tag rate is taken from data directly — so a separate
*mis-tag* SF is generally **not** needed for the all-hadronic search. The SF
that matters is the **signal top-tag efficiency SF** above (applied to TTbar MC
and to the signal model). Confirm this against AN-22-167's systematics list.

**Scope note:** the SF measurement is a **larger, separate task** than the WP
derivation — it needs a new semileptonic selection, lepton/MET objects, b-tag
WPs, background MC, and a fit. This plan delivers the **WPs first**; the SF
measurement should be its own follow-up processor + fit notebook. The reusable
gen-merging classification (3-quark containment) can be prototyped now since
`get_hadronic_tops` and the daughter-quark logic already exist in
`python/truthstudy.py`.

## 7. Storing the outputs (no analysis integration yet)

Do **not** edit `_TAGGER_WPS` for now. Instead write a versioned file:

- `data/toptag/toptag_wp_2024.json` — pT-binned thresholds:
  ```json
  {
    "iov": "2024",
    "discriminant": "globalParT3_TopvsQCD",
    "mass_window": [105, 210],
    "eta_max": 2.5,
    "provenance": {"sample_qcd": "...", "sample_ttbar": "...", "date": "...", "notebook": "notebooks/toptag_wp.ipynb"},
    "working_points": {
      "very_tight": {"target_mistag": 0.001, "pt_bins": [[400,500], ...], "thresholds": [...], "sig_eff": [...]},
      ...
    }
  }
  ```
- A correctionlib-schema JSON is a good forward-compatible option (pT → threshold
  evaluator), but a plain JSON table is fine for a first pass.
- Save the plots under `plots/images/toptag_wp/2024/`.
- When SFs exist later, store them in a parallel `toptag_sf_2024.json`.

## 8. Plots to produce

1. Score `D` distributions — signal vs background, normalized, log-y.
2. ROC curve — `ε_sig` vs `ε_bkg`, log-x; mark the five target points.
3. WP threshold vs pT — `t*` per target mis-tag (the parametrization itself).
4. Signal eff vs pT at each derived WP.
5. Mis-tag eff vs pT at each derived WP (closure: should sit on the target).
6. Mis-tag eff vs `mSD` at a fixed WP (decorrelation cross-check, even though we
   derive *inside* the window — confirms the window choice is benign).
7. *(Optional)* mis-tag eff in the `(pT, |eta|)` plane.

## 9. How to run

> Exact CLI/notebook wiring to be finalized when the code is written; this is the
> intended workflow.

1. **Produce the per-jet tables** (one run per sample group, IOV = 2024):
   - In the runner (`ttbaranalysis.ipynb` / `.py`), select the new
     `TopTagWPProcessor`, `iov="2024"`, and the QCD and TTbar 2024 filesets from
     `data/nanoAOD/*.json`.
   - Local smoke test: `env="local"`, a single small file, a few thousand events.
   - Full run: `env="lpc"` or `coffea.casa` Dask, as with the main analysis.
   - Outputs: `outputs/toptag_wp_qcd_2024.parquet` and
     `outputs/toptag_wp_ttbar_2024.parquet` (or `.coffea` with the accumulated
     tables).
2. **Derive WPs + plots:** open `notebooks/toptag_wp.ipynb`, point it at the two
   tables, run all cells. It writes `data/toptag/toptag_wp_2024.json` and the
   plots under `plots/images/toptag_wp/2024/`.
3. **Add a new year later:** add an entry to `TAGGER_CONFIG`, ensure the
   NanoAOD v15 fileset exists in `data/nanoAOD/`, rerun steps 1–2 with the new
   `iov`.

A short usage section will be appended to `README.md` (or a
`docs/toptag_wp.md`) once the code lands.

## 10. Out of scope / follow-ups

- Wiring derived WPs back into `_TAGGER_WPS` / the analysis selection.
- Full top-tag **scale-factor measurement** (semileptonic tag-and-probe + fit) —
  planned as a separate processor + notebook (§6).
- HOTVR or subjet-b-tag–supplemented tagger variants.
- 2023 (`particleNet_XttVsQCD`) — wiki notes its Run-3 top-tag validity is itself
  unconfirmed; separate follow-up.

---

## Note on local testing / NanoAOD files

I have **no local NanoAOD files**, so I can't run or smoke-test the processor
here. To develop and validate it end-to-end on a small scale, it would help to
have **two small NanoAOD v15 2024 files** placed locally:

- one **QCD** file (e.g. one `QCD-4Jets_HT-*` 2024 v15 file), and
- one **TTbar** 2024 v15 file (so gen-matching + signal eff can be exercised).

A few thousand events each is plenty for a smoke test. They must be **v15** so
the `globalParT3_*` branches are present. Put them somewhere like
`data/test_nanoaod/` and tell me the paths. Without them I can still write the
code, but verification will be limited to static review until a full Dask run.
