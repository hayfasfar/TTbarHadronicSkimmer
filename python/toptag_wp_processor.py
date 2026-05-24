"""Coffea processor for deriving GloParTv3 top-tagging working points.

This processor is intentionally separate from ``TTbarResProcessor``: it does NOT
run the full all-hadronic event selection. Instead it fills per-AK8-jet
histograms of the top-tagger discriminant, which is exactly the population
needed to derive working points (the WP threshold is defined by the QCD
background mis-tag efficiency, and the signal efficiency is read off the same
distribution from gen-matched tops in TTbar).

No ntuples are produced — only histograms.

Year handling is data-driven via ``TAGGER_CONFIG`` so that adding a new IOV
later (once NanoAOD v15 is available for it) is a one-line change. There is no
hard-coded ``'2024'`` in the logic.

Outputs (keys in the accumulator dict):
  - ``score``        : Hist[dataset, jettype, pt, disc] inside the mass window.
                       Used to derive pT-binned WP thresholds and signal eff.
  - ``score_vs_msd`` : Hist[dataset, jettype, msd, disc] WITHOUT the mass window.
                       Used for the mass-decorrelation cross-check.
  - ``sumw``         : defaultdict(float), sum of genWeight per dataset.
  - ``nevents``      : defaultdict(int),   events kept after preprocessing.
  - ``nevents_raw``  : defaultdict(int),   input events before preprocessing.
  - ``qcd_genweight_rejected`` : defaultdict(int), QCD events removed by the
                       large-genWeight filter mirrored from ``TTbarResProcessor``.

``jettype`` is ``"incl"`` (all preselected jets) or ``"matched"`` (AK8 matched
to a gen hadronic top, filled only for signal/TTbar samples). The background
(QCD) mis-tag uses ``"incl"``; the signal efficiency uses ``"matched"``.
"""

import os
import sys

import numpy as np
import awkward as ak
import hist
from coffea import processor

# Mirror the import convention used by ttbarprocessor.py so this works both
# locally and on Dask workers (where the package path is not importable).
sys.path.append(os.getcwd() + '/python/')
from truthstudy import get_hadronic_tops, _ensure_p4  # noqa: E402


# ---------------------------------------------------------------------------
# Per-IOV tagger configuration. To add a year, add an entry here with a
# discriminant function and the NanoAOD FatJet fields it requires.
# ---------------------------------------------------------------------------
def _topvsqcd_globalParT3(fatjets):
    """GloParTv3 mass-decorrelated TopvsQCD discriminant.

    (TopbWqq + TopbWq) / (TopbWqq + TopbWq + QCD)

    Identical to TTbarResProcessor._tscore so derived WPs match the analysis.
    """
    num = fatjets.globalParT3_TopbWqq + fatjets.globalParT3_TopbWq
    den = num + fatjets.globalParT3_QCD
    return num / den


TAGGER_CONFIG = {
    '2024': {
        'score': _topvsqcd_globalParT3,
        'required_fields': [
            'globalParT3_TopbWqq',
            'globalParT3_TopbWq',
            'globalParT3_QCD',
        ],
        'label': 'globalParT3_TopvsQCD',
    },
    # Add future IOVs once NanoAOD v15 exists for them, e.g.:
    # '2025': {'score': _topvsqcd_globalParT3,
    #          'required_fields': [...], 'label': 'globalParT3_TopvsQCD'},
}

# pT bin edges used for the pT-dependent WP derivation (GeV). The last edge is a
# finite stand-in for "infinity".
PT_BIN_EDGES = [400.0, 500.0, 600.0, 800.0, 1200.0, 3000.0]

# Number of fine discriminant bins. A fine binning is required so the cumulative
# tail used to invert for a target mis-tag rate is smooth.
N_DISC_BINS = 1000

# Mass window (GeV) applied when deriving WPs (kept, per analysis convention).
MSD_MIN = 105.0
MSD_MAX = 210.0

# Jet pre-selection.
JET_PT_MIN = 400.0
JET_ETA_MAX = 2.5

# gen-top match radius.
DR_MATCH = 0.8


def _qcd_genweight_mask(gen_weight, nsigma=2.0):
    """Mirror TTbarResProcessor's QCD large-genWeight event rejection."""
    vals = ak.to_numpy(gen_weight)
    if len(vals) == 0:
        return ak.ones_like(gen_weight, dtype=bool)

    average = np.average(vals)
    stddev = np.std(vals)
    if stddev == 0 or not np.isfinite(stddev):
        return ak.ones_like(gen_weight, dtype=bool)

    return np.abs((gen_weight - average) / stddev) < nsigma


def _make_score_hist():
    return hist.Hist(
        hist.axis.StrCategory([], name="dataset", growth=True),
        hist.axis.StrCategory([], name="jettype", growth=True),
        hist.axis.Variable(PT_BIN_EDGES, name="pt", label=r"AK8 $p_T$ [GeV]"),
        hist.axis.Regular(N_DISC_BINS, 0.0, 1.0, name="disc", label="TopvsQCD"),
        storage="weight",
        name="Counts",
    )


def _make_score_vs_msd_hist():
    return hist.Hist(
        hist.axis.StrCategory([], name="dataset", growth=True),
        hist.axis.StrCategory([], name="jettype", growth=True),
        hist.axis.Regular(60, 0.0, 300.0, name="msd", label=r"$m_{SD}$ [GeV]"),
        hist.axis.Regular(200, 0.0, 1.0, name="disc", label="TopvsQCD"),
        storage="weight",
        name="Counts",
    )


class TopTagWPProcessor(processor.ProcessorABC):
    """Fill per-AK8-jet tagger-discriminant histograms for WP derivation.

    Parameters
    ----------
    iov : str
        IOV key into ``TAGGER_CONFIG`` (e.g. ``'2024'``). No year is hard-coded
        in the logic; an unknown IOV raises immediately.
    match_gen_top : bool or None
        If True, also fill the ``"matched"`` jettype (AK8 matched to a gen
        hadronic top). If None (default), auto-detect from sample metadata /
        dataset name (samples starting with ``TT``).
    sample_metadata : dict or None
        Optional per-sample metadata (sample, subsample, year, is_mc, xsec_pb).
    """

    def __init__(self, iov='2024', match_gen_top=None, sample_metadata=None):
        if iov not in TAGGER_CONFIG:
            raise KeyError(
                f"IOV '{iov}' not in TAGGER_CONFIG (known: {list(TAGGER_CONFIG)}). "
                "Add an entry with a discriminant + required NanoAOD fields."
            )
        self.iov = iov
        self._cfg = TAGGER_CONFIG[iov]
        self.match_gen_top = match_gen_top
        self.sample_metadata = dict(sample_metadata or {})

    # -- helpers ------------------------------------------------------------
    def _should_match(self, dataset):
        if self.match_gen_top is not None:
            return bool(self.match_gen_top)
        sample = str(self.sample_metadata.get('sample', '')).upper()
        return sample.startswith('TT') or 'TT' in str(dataset).upper()

    @staticmethod
    def _matched_to_gen_top(events):
        """Return a per-jet boolean: AK8 matched to a gen hadronic top (dR<0.8)."""
        tops = _ensure_p4(get_hadronic_tops(events.GenPart))
        fatjets = _ensure_p4(events.FatJet)
        pairs = ak.cartesian({"fat": fatjets, "top": tops}, axis=1, nested=True)
        dr = pairs["fat"].p4.delta_r(pairs["top"].p4)
        min_dr = ak.fill_none(ak.min(dr, axis=2), 999.0)
        return min_dr < DR_MATCH

    # -- coffea API ---------------------------------------------------------
    def process(self, events):
        dataset = events.metadata['dataset']
        cfg = self._cfg

        output = {
            'score': _make_score_hist(),
            'score_vs_msd': _make_score_vs_msd_hist(),
            'sumw': processor.defaultdict_accumulator(float),
            'nevents': processor.defaultdict_accumulator(int),
            'nevents_raw': processor.defaultdict_accumulator(int),
            'qcd_genweight_rejected': processor.defaultdict_accumulator(int),
        }

        n_raw = len(events)
        n_rejected = 0
        if 'QCD' in dataset and 'genWeight' in events.fields:
            genweight_mask = _qcd_genweight_mask(events.genWeight)
            n_rejected = n_raw - int(ak.sum(genweight_mask))
            events = events[genweight_mask]

        fj = events.FatJet
        missing = [f for f in cfg['required_fields'] if f not in fj.fields]
        if missing:
            raise RuntimeError(
                f"FatJet missing {missing} for IOV {self.iov} (dataset {dataset}). "
                "Confirm the sample is NanoAOD v15 with GloParTv3 branches."
            )

        # Per-event weight (genWeight for MC, 1.0 for data / when absent).
        if 'genWeight' in events.fields:
            w_evt = events.genWeight
        else:
            w_evt = ak.ones_like(ak.num(fj, axis=1), dtype=np.float64) * 1.0

        output['sumw'][dataset] += float(ak.sum(w_evt))
        output['nevents'][dataset] += int(len(events))
        output['nevents_raw'][dataset] += int(n_raw)
        output['qcd_genweight_rejected'][dataset] += int(n_rejected)

        disc = cfg['score'](fj)
        pt = fj.pt
        eta = fj.eta
        msd = fj.msoftdrop

        presel = (pt > JET_PT_MIN) & (abs(eta) < JET_ETA_MAX)
        window = presel & (msd > MSD_MIN) & (msd < MSD_MAX)

        # broadcast event weight to per-jet
        w_jet = ak.broadcast_arrays(w_evt, pt)[0]

        def _fill(h, jettype, mask, extra_axes):
            sel_disc = ak.to_numpy(ak.flatten(disc[mask]))
            sel_w = ak.to_numpy(ak.flatten(w_jet[mask]))
            kw = {name: ak.to_numpy(ak.flatten(arr[mask])) for name, arr in extra_axes.items()}
            h.fill(dataset=dataset, jettype=jettype, disc=sel_disc, weight=sel_w, **kw)

        # inclusive (background mis-tag denominator/numerator)
        _fill(output['score'], 'incl', window, {'pt': pt})
        _fill(output['score_vs_msd'], 'incl', presel, {'msd': msd})

        # gen-matched tops (signal efficiency) — only for signal/TTbar samples
        if self._should_match(dataset) and 'GenPart' in events.fields:
            is_matched = self._matched_to_gen_top(events)
            _fill(output['score'], 'matched', window & is_matched, {'pt': pt})
            _fill(output['score_vs_msd'], 'matched', presel & is_matched, {'msd': msd})

        return output

    def postprocess(self, accumulator):
        return accumulator
