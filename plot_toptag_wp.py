#!/usr/bin/env python
"""Derive GloParTv3 top-tagging working points and make plots.

Reads the histogram accumulator written by ``run_toptag_wp.py`` and:

  1. For each AK8 pT bin, inverts the QCD background mis-tag efficiency to find
     the TopvsQCD threshold for each target mis-tag rate (the "working point"),
     using the same convention as the official top-tag table:

         very_tight 0.1% | tight 0.5% | medium 1.0% | loose 2.5% | very_loose 5.0%

  2. Reads the signal (gen-matched top) efficiency at each threshold.
  3. Writes the pT-binned thresholds + efficiencies to a JSON file.
  4. Produces validation plots (score distributions, ROC, threshold vs pT,
     signal-eff vs pT, mis-tag closure vs pT, and mis-tag vs mSD decorrelation).

The WP threshold is defined purely by the QCD background; the signal sample only
supplies the efficiency reported alongside each WP.

Usage::

    coffea-dask/bin/python plot_toptag_wp.py outputs/toptag_wp_2024.coffea \
        --iov 2024 --json data/toptag/toptag_wp_2024.json \
        --plotdir plots/images/toptag_wp/2024
"""

import os
import json
import argparse
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import mplhep as hep
    plt.style.use(hep.style.CMS)
except Exception:
    hep = None

from coffea import util

# Target mis-tag rates, named to mirror the official top-tag table.
TARGETS = {
    'very_tight': 0.001,
    'tight':      0.005,
    'medium':     0.010,
    'loose':      0.025,
    'very_loose': 0.050,
}
TARGET_ORDER = ['very_tight', 'tight', 'medium', 'loose', 'very_loose']

# CMS style uses large axis labels by default, so keep canvases roomy enough
# that labels, ticks, legends, and the CMS header do not dominate the plot.
SINGLE_PANEL_FIGSIZE = (10, 8)
SCORE_PANEL_SIZE = (6.7, 5.2)
DATAMC_PANEL_SIZE = (6.7, 6.4)
DEFAULT_SCORE_REBIN = 10
QCD_MIN_SUBSAMPLE_PT = 300.0
QCD_SUBSAMPLE_RE = re.compile(r'QCD_(?:Bin-)?PT-?(\d+(?:\.\d+)?)to')


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _is_data_metadata(m):
    return (not m.get('is_mc', True)) or str(m.get('sample', '')).lower() in {'data', 'jetmet'}


def qcd_subsample_min_pt(name):
    match = QCD_SUBSAMPLE_RE.search(str(name))
    return float(match.group(1)) if match else None


def keep_qcd_subsample(ds, metadata, min_pt=QCD_MIN_SUBSAMPLE_PT):
    """Keep QCD generated-pT bins starting at ``min_pt`` GeV."""
    sample = str(metadata.get('sample', ds)).upper()
    if not sample.startswith('QCD') and 'QCD' not in str(ds).upper():
        return True

    subsample = metadata.get('subsample', ds)
    low = qcd_subsample_min_pt(subsample)
    return low is None or low >= min_pt


def classify_datasets(meta):
    """Return (signal_datasets, background_datasets) from datasets_metadata."""
    sig, bkg = [], []
    for ds, m in meta.items():
        if _is_data_metadata(m):
            continue
        if not keep_qcd_subsample(ds, m):
            continue
        if str(m.get('sample', ds)).upper().startswith('TT'):
            sig.append(ds)
        else:
            bkg.append(ds)
    return sig, bkg


def classify_data_datasets(meta):
    """Return datasets marked as collision data."""
    return [ds for ds, m in meta.items() if _is_data_metadata(m)]


def classify_mc_datasets(meta):
    """Return MC datasets used for inclusive Data/MC score comparisons."""
    return [
        ds for ds, m in meta.items()
        if not _is_data_metadata(m) and keep_qcd_subsample(ds, m)
    ]


def metadata_with_sumw(output):
    meta = dict(output.get('datasets_metadata', {}))
    for ds in meta:
        meta[ds]['_sumw'] = float(output['sumw'].get(ds, 0.0))
    return meta


def dataset_scale(ds, meta, lumi_pb):
    """lumi*xsec/sumw if available, else 1.0 (irrelevant for a single sample)."""
    m = meta.get(ds, {})
    xsec = m.get('xsec_pb')
    sumw = m.get('_sumw')
    if xsec and sumw and lumi_pb:
        return lumi_pb * xsec / sumw
    return 1.0


def combine_disc(hscore, datasets, jettype, pt_index, scales):
    """Sum (value, variance) over datasets for one jettype & pt bin -> disc dist."""
    val = var = None
    for ds in datasets:
        view = hscore[{'dataset': ds, 'jettype': jettype, 'pt': pt_index}].view(flow=False)
        s = scales.get(ds, 1.0)
        v = view['value'] * s
        e = view['variance'] * (s ** 2)
        val = v if val is None else val + v
        var = e if var is None else var + e
    return val, var


def combine_1d(histo, datasets, jettype, scales):
    """Sum a 1D weighted hist over datasets for one jettype."""
    val = var = None
    for ds in datasets:
        view = histo[{'dataset': ds, 'jettype': jettype}].view(flow=False)
        s = scales.get(ds, 1.0)
        v = view['value'] * s
        e = view['variance'] * (s ** 2)
        val = v if val is None else val + v
        var = e if var is None else var + e
    return val, var


def tail_fraction(counts, edges):
    """Return (frac_at_edges, edges): fraction of weight with disc >= each edge."""
    total = counts.sum()
    if total <= 0:
        return np.zeros(len(edges)), edges
    tail = np.cumsum(counts[::-1])[::-1]            # tail[k] = sum(counts[k:])
    frac_left = tail / total                         # at edges[:-1]
    frac = np.append(frac_left, 0.0)                 # add eff=0 at the last edge
    return frac, total


def invert_threshold(frac_at_edges, edges, target):
    """Find threshold t where fraction(disc >= t) == target (monotone interp)."""
    # frac decreasing with edge; reverse so x is increasing for np.interp
    return float(np.interp(target, frac_at_edges[::-1], edges[::-1]))


def eff_at_threshold(counts, edges, threshold):
    """Fraction of weight with disc >= threshold."""
    total = counts.sum()
    if total <= 0:
        return 0.0
    tail = np.cumsum(counts[::-1])[::-1]
    frac = np.append(tail / total, 0.0)
    return float(np.interp(threshold, edges, frac))


def effective_entries(counts, variance):
    """Neff = (sum w)^2 / sum(w^2) — for binomial efficiency uncertainty."""
    sw = counts.sum()
    sw2 = variance.sum()
    return (sw ** 2 / sw2) if sw2 > 0 else 0.0


def rebin_counts(counts, edges, factor):
    """Merge neighboring 1D bins for display without changing total weight."""
    if factor <= 1:
        return counts, edges

    rebinned = []
    rebinned_edges = [edges[0]]
    for start in range(0, len(counts), factor):
        stop = min(start + factor, len(counts))
        rebinned.append(counts[start:stop].sum())
        rebinned_edges.append(edges[stop])
    return np.asarray(rebinned), np.asarray(rebinned_edges)


def mc_shape_scale_to_data(data_total, mc_total):
    """Per-panel scale for Data/MC shape comparisons."""
    if mc_total <= 0:
        return 0.0
    return data_total / mc_total


def data_mc_ratio(data_counts, data_variance, mc_counts):
    """Return Data/MC ratio and data statistical uncertainty."""
    ratio = np.full_like(data_counts, np.nan, dtype=float)
    ratio_err = np.full_like(data_counts, np.nan, dtype=float)
    mask = mc_counts > 0
    ratio[mask] = data_counts[mask] / mc_counts[mask]
    ratio_err[mask] = np.sqrt(data_variance[mask]) / mc_counts[mask]
    return ratio, ratio_err


# ---------------------------------------------------------------------------
# derivation
# ---------------------------------------------------------------------------
def derive(output, iov, lumi_pb):
    hscore = output['score']
    meta = metadata_with_sumw(output)
    sig_ds, bkg_ds = classify_datasets(meta)
    scales = {ds: dataset_scale(ds, meta, lumi_pb) for ds in meta}

    pt_axis = hscore.axes['pt']
    disc_edges = hscore.axes['disc'].edges
    pt_edges = pt_axis.edges

    result = {
        'iov': iov,
        'discriminant': output.get('run_info', {}).get('tagger_label', 'TopvsQCD'),
        'mass_window': [105.0, 210.0],
        'eta_max': 2.5,
        'pt_min': 400.0,
        'signal_datasets': sig_ds,
        'background_datasets': bkg_ds,
        'pt_bin_edges': pt_edges.tolist(),
        'working_points': {},
    }

    # per-pt-bin distributions cached for plotting
    cache = {'pt_edges': pt_edges, 'disc_edges': disc_edges,
             'bkg': [], 'sig': [], 'bkg_var': [], 'sig_var': []}

    for wp in TARGET_ORDER:
        result['working_points'][wp] = {
            'target_mistag': TARGETS[wp],
            'pt_bins': [], 'threshold': [], 'threshold_unc': [],
            'mistag_achieved': [], 'signal_eff': [],
        }

    for i in range(pt_axis.size):
        bkg_c, bkg_v = combine_disc(hscore, bkg_ds, 'incl', i, scales)
        sig_c, sig_v = combine_disc(hscore, sig_ds, 'matched', i, scales)
        cache['bkg'].append(bkg_c); cache['bkg_var'].append(bkg_v)
        cache['sig'].append(sig_c); cache['sig_var'].append(sig_v)

        frac_bkg, bkg_total = tail_fraction(bkg_c, disc_edges)
        neff = effective_entries(bkg_c, bkg_v)
        ptlo, pthi = float(pt_edges[i]), float(pt_edges[i + 1])

        for wp in TARGET_ORDER:
            tgt = TARGETS[wp]
            wpd = result['working_points'][wp]
            wpd['pt_bins'].append([ptlo, pthi])
            if bkg_total <= 0 or neff <= 0:
                wpd['threshold'].append(None); wpd['threshold_unc'].append(None)
                wpd['mistag_achieved'].append(None); wpd['signal_eff'].append(None)
                continue
            t_star = invert_threshold(frac_bkg, disc_edges, tgt)
            # binomial uncertainty on the target eff -> threshold half-width
            err = np.sqrt(max(tgt * (1 - tgt) / neff, 0.0))
            t_hi = invert_threshold(frac_bkg, disc_edges, max(tgt - err, 0.0))
            t_lo = invert_threshold(frac_bkg, disc_edges, min(tgt + err, 1.0))
            t_unc = 0.5 * abs(t_hi - t_lo)
            mistag_ach = eff_at_threshold(bkg_c, disc_edges, t_star)
            sig_eff = eff_at_threshold(sig_c, disc_edges, t_star) if sig_c is not None and sig_c.sum() > 0 else None
            wpd['threshold'].append(round(t_star, 5))
            wpd['threshold_unc'].append(round(t_unc, 5))
            wpd['mistag_achieved'].append(round(mistag_ach, 6))
            wpd['signal_eff'].append(round(sig_eff, 5) if sig_eff is not None else None)

    return result, cache, (sig_ds, bkg_ds, scales)


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
def _cms(ax, iov):
    if hep is not None:
        hep.cms.label("Preliminary", data=False, year=iov, ax=ax, fontsize=14)


def plot_score_dists(cache, iov, plotdir, score_rebin=DEFAULT_SCORE_REBIN):
    edges = cache['disc_edges']
    pt_edges = cache['pt_edges']
    n = len(cache['bkg'])
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(SCORE_PANEL_SIZE[0] * ncol, SCORE_PANEL_SIZE[1] * nrow),
        squeeze=False,
    )
    for i in range(nrow * ncol):
        ax = axes[i // ncol][i % ncol]
        if i >= n:
            ax.axis('off'); continue
        bkg = cache['bkg'][i]; sig = cache['sig'][i]
        if bkg.sum() > 0:
            bkg_plot, plot_edges = rebin_counts(bkg, edges, score_rebin)
            ax.stairs(bkg_plot / bkg_plot.sum(), plot_edges, label='QCD (bkg)', color='C3')
        if sig is not None and sig.sum() > 0:
            sig_plot, plot_edges = rebin_counts(sig, edges, score_rebin)
            ax.stairs(sig_plot / sig_plot.sum(), plot_edges, label='TTbar matched (sig)', color='C0')
        ax.set_yscale('log')
        ax.set_xlabel('TopvsQCD'); ax.set_ylabel('a.u.')
        ax.set_title(f'{pt_edges[i]:.0f} < pT < {pt_edges[i+1]:.0f} GeV', fontsize=11)
        ax.legend(fontsize=9)
    fig.tight_layout()
    p = os.path.join(plotdir, 'score_distributions.png')
    fig.savefig(p, dpi=120); plt.close(fig); return p


def _load_data_output(data_infile):
    if not data_infile:
        return None
    return util.load(data_infile)


def data_mc_score_hist(output):
    """Score hist for Data/MC shape plots; prefer full mSD sidebands."""
    return output.get('score_full_msd', output['score'])


def plot_data_mc_score_dists(mc_output, data_output, iov, plotdir,
                             score_rebin=DEFAULT_SCORE_REBIN):
    """TopvsQCD Data/MC shape plot with stacked MC scaled per pT bin."""
    if data_output is None:
        data_output = mc_output

    hmc = data_mc_score_hist(mc_output)
    hdata = data_mc_score_hist(data_output)
    mc_meta = metadata_with_sumw(mc_output)
    data_meta = metadata_with_sumw(data_output)
    ttbar_ds, qcd_ds = classify_datasets(mc_meta)
    data_ds = classify_data_datasets(data_meta)
    if not data_ds or (not qcd_ds and not ttbar_ds):
        return None

    lumi_pb = mc_output.get('run_info', {}).get('lumi_pb')
    mc_scales = {ds: dataset_scale(ds, mc_meta, lumi_pb) for ds in mc_meta}
    data_scales = {ds: 1.0 for ds in data_meta}
    edges = hmc.axes['disc'].edges
    pt_edges = hmc.axes['pt'].edges
    n = hmc.axes['pt'].size
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig = plt.figure(figsize=(DATAMC_PANEL_SIZE[0] * ncol, DATAMC_PANEL_SIZE[1] * nrow))
    outer = fig.add_gridspec(nrow, ncol, hspace=0.40, wspace=0.42)
    for i in range(nrow * ncol):
        inner = outer[i // ncol, i % ncol].subgridspec(
            2,
            1,
            height_ratios=(3.0, 1.0),
            hspace=0.05,
        )
        ax = fig.add_subplot(inner[0])
        rax = fig.add_subplot(inner[1], sharex=ax)
        if i >= n:
            ax.axis('off')
            rax.axis('off')
            continue

        qcd_c, qcd_v = combine_disc(hmc, qcd_ds, 'incl', i, mc_scales)
        ttbar_c, ttbar_v = combine_disc(hmc, ttbar_ds, 'incl', i, mc_scales)
        data_c, data_v = combine_disc(hdata, data_ds, 'incl', i, data_scales)

        if data_c is not None and data_c.sum() > 0:
            data_plot, plot_edges = rebin_counts(data_c, edges, score_rebin)
            data_var, _ = rebin_counts(data_v, edges, score_rebin)
            centers = 0.5 * (plot_edges[:-1] + plot_edges[1:])
            widths = np.diff(plot_edges)

            if qcd_c is not None and qcd_c.sum() > 0:
                qcd_plot, _ = rebin_counts(qcd_c, edges, score_rebin)
                qcd_var, _ = rebin_counts(qcd_v, edges, score_rebin)
            else:
                qcd_plot = np.zeros_like(data_plot)
                qcd_var = np.zeros_like(data_plot)
            if ttbar_c is not None and ttbar_c.sum() > 0:
                ttbar_plot, _ = rebin_counts(ttbar_c, edges, score_rebin)
                ttbar_var, _ = rebin_counts(ttbar_v, edges, score_rebin)
            else:
                ttbar_plot = np.zeros_like(data_plot)
                ttbar_var = np.zeros_like(data_plot)

            mc_shape_scale = mc_shape_scale_to_data(
                data_plot.sum(),
                qcd_plot.sum() + ttbar_plot.sum(),
            )
            qcd_plot = qcd_plot * mc_shape_scale
            qcd_var = qcd_var * (mc_shape_scale ** 2)
            ttbar_plot = ttbar_plot * mc_shape_scale
            ttbar_var = ttbar_var * (mc_shape_scale ** 2)
            bottom = np.zeros_like(qcd_plot)
            ax.bar(
                plot_edges[:-1],
                qcd_plot,
                width=widths,
                align='edge',
                bottom=bottom,
                label='QCD',
                color='C3',
                alpha=0.65,
                linewidth=0,
            )
            bottom = bottom + qcd_plot
            ax.bar(
                plot_edges[:-1],
                ttbar_plot,
                width=widths,
                align='edge',
                bottom=bottom,
                label='TTbar',
                color='C0',
                alpha=0.65,
                linewidth=0,
            )
            ax.plot([], [], ' ', label=f'MC shape scale {mc_shape_scale:.2g}')

            mc_total = qcd_plot + ttbar_plot
            mc_var = qcd_var + ttbar_var
            mc_err = np.sqrt(mc_var)
            band_bottom = np.maximum(mc_total - mc_err, 0.0)
            band_top = mc_total + mc_err
            ax.bar(
                plot_edges[:-1],
                band_top - band_bottom,
                width=widths,
                align='edge',
                bottom=band_bottom,
                label='MC stat. unc.',
                facecolor='none',
                edgecolor='0.35',
                hatch='////',
                linewidth=0,
            )

            yerr = np.sqrt(data_var)
            ax.errorbar(centers, data_plot, yerr=yerr, fmt='o', ms=3, lw=1,
                        label='Data', color='black')

            ratio, ratio_err = data_mc_ratio(data_plot, data_var, mc_total)
            rel_mc = np.full_like(mc_total, np.nan, dtype=float)
            mc_mask = mc_total > 0
            rel_mc[mc_mask] = mc_err[mc_mask] / mc_total[mc_mask]
            ratio_band_bottom = np.maximum(1.0 - rel_mc, 0.0)
            ratio_band_top = 1.0 + rel_mc
            finite_band = np.isfinite(ratio_band_bottom) & np.isfinite(ratio_band_top)
            if np.any(finite_band):
                rax.bar(
                    plot_edges[:-1][finite_band],
                    (ratio_band_top - ratio_band_bottom)[finite_band],
                    width=widths[finite_band],
                    align='edge',
                    bottom=ratio_band_bottom[finite_band],
                    facecolor='none',
                    edgecolor='0.35',
                    hatch='////',
                    linewidth=0,
                )
            rax.errorbar(centers, ratio, yerr=ratio_err, fmt='o', ms=3, lw=1,
                         color='black')

        ax.set_yscale('log')
        ax.set_ylabel('Events / bin', fontsize=18, labelpad=4)
        ax.set_title(f'{pt_edges[i]:.0f} < pT < {pt_edges[i+1]:.0f} GeV', fontsize=12)
        ax.tick_params(labelsize=14)
        ax.legend(fontsize=9)
        plt.setp(ax.get_xticklabels(), visible=False)
        rax.axhline(1.0, color='0.35', lw=1, ls='--')
        rax.set_xlabel('TopvsQCD', fontsize=20, labelpad=2)
        rax.set_ylabel('Data/MC', fontsize=16, labelpad=3)
        rax.tick_params(labelsize=14)
        rax.set_ylim(0.0, 2.0)
        rax.grid(axis='y', color='0.85', lw=0.7)

    p = os.path.join(plotdir, 'data_mc_score_distributions.png')
    fig.savefig(p, dpi=120); plt.close(fig); return p


def plot_pt_dists(mc_output, data_output, iov, plotdir):
    """Preselected AK8 pT control plot for QCD stitching / smoothness checks."""
    if 'jet_pt' not in mc_output:
        return None
    if data_output is not None and 'jet_pt' not in data_output:
        data_output = None

    hmc = mc_output['jet_pt']
    mc_meta = metadata_with_sumw(mc_output)
    data_meta = metadata_with_sumw(data_output) if data_output is not None else {}
    ttbar_ds, qcd_ds = classify_datasets(mc_meta)
    data_ds = classify_data_datasets(data_meta)
    if not qcd_ds and not ttbar_ds:
        return None

    lumi_pb = mc_output.get('run_info', {}).get('lumi_pb')
    mc_scales = {ds: dataset_scale(ds, mc_meta, lumi_pb) for ds in mc_meta}
    data_scales = {ds: 1.0 for ds in data_meta}
    edges = hmc.axes['pt'].edges
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    qcd_plot, qcd_var = combine_1d(hmc, qcd_ds, 'incl', mc_scales)
    ttbar_plot, ttbar_var = combine_1d(hmc, ttbar_ds, 'incl', mc_scales)
    if qcd_plot is None:
        qcd_plot = np.zeros(len(edges) - 1)
        qcd_var = np.zeros_like(qcd_plot)
    if ttbar_plot is None:
        ttbar_plot = np.zeros(len(edges) - 1)
        ttbar_var = np.zeros_like(ttbar_plot)

    data_plot = data_var = None
    mc_shape_scale = None
    if data_output is not None and data_ds:
        data_plot, data_var = combine_1d(data_output['jet_pt'], data_ds, 'incl', data_scales)
        if data_plot is not None and data_plot.sum() > 0:
            mc_shape_scale = mc_shape_scale_to_data(
                data_plot.sum(),
                qcd_plot.sum() + ttbar_plot.sum(),
            )
            qcd_plot = qcd_plot * mc_shape_scale
            qcd_var = qcd_var * (mc_shape_scale ** 2)
            ttbar_plot = ttbar_plot * mc_shape_scale
            ttbar_var = ttbar_var * (mc_shape_scale ** 2)
        else:
            data_plot = data_var = None

    has_data = data_plot is not None
    if has_data:
        fig = plt.figure(figsize=(10, 8))
        gs = fig.add_gridspec(2, 1, height_ratios=(3.0, 1.0), hspace=0.05)
        ax = fig.add_subplot(gs[0])
        rax = fig.add_subplot(gs[1], sharex=ax)
    else:
        fig, ax = plt.subplots(figsize=SINGLE_PANEL_FIGSIZE)
        rax = None

    bottom = np.zeros_like(qcd_plot)
    ax.bar(edges[:-1], qcd_plot, width=widths, align='edge', bottom=bottom,
           label='QCD', color='C3', alpha=0.65, linewidth=0)
    bottom = bottom + qcd_plot
    ax.bar(edges[:-1], ttbar_plot, width=widths, align='edge', bottom=bottom,
           label='TTbar', color='C0', alpha=0.65, linewidth=0)

    mc_total = qcd_plot + ttbar_plot
    mc_var = qcd_var + ttbar_var
    mc_err = np.sqrt(mc_var)
    band_bottom = np.maximum(mc_total - mc_err, 0.0)
    band_top = mc_total + mc_err
    ax.bar(edges[:-1], band_top - band_bottom, width=widths, align='edge',
           bottom=band_bottom, label='MC stat. unc.', facecolor='none',
           edgecolor='0.35', hatch='////', linewidth=0)

    if has_data:
        ax.plot([], [], ' ', label=f'MC shape scale {mc_shape_scale:.2g}')
        ax.errorbar(centers, data_plot, yerr=np.sqrt(data_var), fmt='o',
                    ms=3, lw=1, label='Data', color='black')
        ratio, ratio_err = data_mc_ratio(data_plot, data_var, mc_total)
        rel_mc = np.full_like(mc_total, np.nan, dtype=float)
        mc_mask = mc_total > 0
        rel_mc[mc_mask] = mc_err[mc_mask] / mc_total[mc_mask]
        ratio_band_bottom = np.maximum(1.0 - rel_mc, 0.0)
        ratio_band_top = 1.0 + rel_mc
        finite_band = np.isfinite(ratio_band_bottom) & np.isfinite(ratio_band_top)
        if np.any(finite_band):
            rax.bar(
                edges[:-1][finite_band],
                (ratio_band_top - ratio_band_bottom)[finite_band],
                width=widths[finite_band],
                align='edge',
                bottom=ratio_band_bottom[finite_band],
                facecolor='none',
                edgecolor='0.35',
                hatch='////',
                linewidth=0,
            )
        rax.errorbar(centers, ratio, yerr=ratio_err, fmt='o', ms=3, lw=1,
                     color='black')
        rax.axhline(1.0, color='0.35', lw=1, ls='--')
        rax.set_ylabel('Data/MC')
        rax.set_xlabel('AK8 pT [GeV]')
        rax.set_ylim(0.0, 2.0)
        rax.grid(axis='y', color='0.85', lw=0.7)
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        ax.set_xlabel('AK8 pT [GeV]')

    ax.set_yscale('log')
    ax.set_ylabel('Events / bin')
    ax.legend(fontsize=9)
    _cms(ax, iov)
    fig.tight_layout()
    p = os.path.join(plotdir, 'pt_distributions.png')
    fig.savefig(p, dpi=120); plt.close(fig); return p


def plot_roc(cache, iov, plotdir):
    edges = cache['disc_edges']
    pt_edges = cache['pt_edges']
    fig, ax = plt.subplots(figsize=SINGLE_PANEL_FIGSIZE)
    for i in range(len(cache['bkg'])):
        bkg = cache['bkg'][i]; sig = cache['sig'][i]
        if bkg.sum() <= 0 or sig is None or sig.sum() <= 0:
            continue
        eb = np.append(np.cumsum(bkg[::-1])[::-1] / bkg.sum(), 0.0)
        es = np.append(np.cumsum(sig[::-1])[::-1] / sig.sum(), 0.0)
        ax.plot(eb, es, label=f'{pt_edges[i]:.0f}-{pt_edges[i+1]:.0f} GeV')
    for tgt in TARGETS.values():
        ax.axvline(tgt, color='grey', ls=':', lw=0.8)
    ax.set_xscale('log')
    ax.set_xlabel('QCD mis-tag efficiency')
    ax.set_ylabel('Signal (matched top) efficiency')
    ax.set_xlim(1e-4, 1); ax.set_ylim(0, 1)
    ax.legend(fontsize=9, title='AK8 pT')
    _cms(ax, iov)
    fig.tight_layout()
    p = os.path.join(plotdir, 'roc.png')
    fig.savefig(p, dpi=120); plt.close(fig); return p


def plot_vs_pt(result, key, ylabel, fname, iov, plotdir, logy=False, target_lines=False):
    fig, ax = plt.subplots(figsize=SINGLE_PANEL_FIGSIZE)
    for wp in TARGET_ORDER:
        wpd = result['working_points'][wp]
        x = [0.5 * (b[0] + b[1]) for b in wpd['pt_bins']]
        y = wpd[key]
        xs = [xx for xx, yy in zip(x, y) if yy is not None]
        ys = [yy for yy in y if yy is not None]
        line, = ax.plot(xs, ys, marker='o', label=f"{wp} ({TARGETS[wp]*100:.1f}%)")
        if target_lines:
            ax.axhline(TARGETS[wp], color=line.get_color(), ls=':', lw=0.8)
    if logy:
        ax.set_yscale('log')
    ax.set_xlabel('AK8 pT [GeV]'); ax.set_ylabel(ylabel)
    ax.legend(fontsize=9)
    _cms(ax, iov)
    fig.tight_layout()
    p = os.path.join(plotdir, fname)
    fig.savefig(p, dpi=120); plt.close(fig); return p


def plot_mistag_vs_msd(output, result, iov, plotdir, wp='medium'):
    """Decorrelation cross-check: mis-tag vs mSD at a fixed (pT-integrated) WP."""
    hms = output['score_vs_msd']
    meta = metadata_with_sumw(output)
    _, bkg_ds = classify_datasets(meta)
    lumi_pb = output.get('run_info', {}).get('lumi_pb')
    scales = {ds: dataset_scale(ds, meta, lumi_pb) for ds in meta}

    # pT-integrated medium threshold from the windowed score hist
    hscore = output['score']
    disc_edges = hscore.axes['disc'].edges
    bkg_c = None
    for ds in bkg_ds:
        v = hscore[{'dataset': ds, 'jettype': 'incl'}][{'pt': sum}].view(flow=False)
        vals = v['value'] * scales.get(ds, 1.0)
        bkg_c = vals if bkg_c is None else bkg_c + vals
    frac, _ = tail_fraction(bkg_c, disc_edges)
    t_star = invert_threshold(frac, disc_edges, TARGETS[wp])

    # mis-tag vs mSD using score_vs_msd (no mass window)
    msd_edges = hms.axes['msd'].edges
    d_edges = hms.axes['disc'].edges
    d_centers = 0.5 * (d_edges[:-1] + d_edges[1:])
    pass_mask = d_centers >= t_star
    denom = np.zeros(len(msd_edges) - 1)
    numer = np.zeros(len(msd_edges) - 1)
    for ds in bkg_ds:
        h2 = hms[{'dataset': ds, 'jettype': 'incl'}].view(flow=False)['value']  # [msd, disc]
        scale = scales.get(ds, 1.0)
        denom += h2.sum(axis=1) * scale
        numer += h2[:, pass_mask].sum(axis=1) * scale
    with np.errstate(divide='ignore', invalid='ignore'):
        mistag = np.where(denom > 0, numer / denom, np.nan)
    centers = 0.5 * (msd_edges[:-1] + msd_edges[1:])
    fig, ax = plt.subplots(figsize=SINGLE_PANEL_FIGSIZE)
    ax.step(centers, mistag, where='mid', color='C3')
    ax.axhline(TARGETS[wp], color='grey', ls='--', label=f'target {TARGETS[wp]*100:.1f}%')
    ax.axvspan(105, 210, color='C0', alpha=0.1, label='mass window')
    ax.set_xlabel(r'$m_{SD}$ [GeV]'); ax.set_ylabel('QCD mis-tag efficiency')
    ax.text(0.03, 0.05, f'{wp} WP (pT-integrated, thr={t_star:.3f})',
            transform=ax.transAxes, fontsize=10)
    ax.legend(fontsize=9)
    _cms(ax, iov)
    fig.tight_layout()
    p = os.path.join(plotdir, 'mistag_vs_msd.png')
    fig.savefig(p, dpi=120); plt.close(fig); return p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('infile')
    ap.add_argument('--iov', default='2024')
    ap.add_argument('--json', default=None)
    ap.add_argument('--plotdir', default=None)
    ap.add_argument('--data-infile', default=None,
                    help='optional data-only .coffea output for Data/MC TopvsQCD plots')
    ap.add_argument('--score-rebin', type=int, default=DEFAULT_SCORE_REBIN,
                    help='merge this many fine TopvsQCD bins in score_distributions.png only')
    args = ap.parse_args()

    out_json = args.json or f'data/toptag/toptag_wp_{args.iov}.json'
    plotdir = args.plotdir or f'plots/images/toptag_wp/{args.iov}'
    os.makedirs(os.path.dirname(out_json) or '.', exist_ok=True)
    os.makedirs(plotdir, exist_ok=True)

    output = util.load(args.infile)
    data_output = _load_data_output(args.data_infile)
    lumi_pb = output.get('run_info', {}).get('lumi_pb')

    result, cache, _ = derive(output, args.iov, lumi_pb)
    result['provenance'] = {
        'input': os.path.abspath(args.infile),
        'lumi_pb': lumi_pb,
        'note': 'WP threshold defined by QCD mis-tag; signal eff from gen-matched TTbar tops.',
    }

    with open(out_json, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"wrote {out_json}")

    plots = [
        plot_score_dists(cache, args.iov, plotdir, score_rebin=args.score_rebin),
        plot_data_mc_score_dists(
            output,
            data_output,
            args.iov,
            plotdir,
            score_rebin=args.score_rebin,
        ),
        plot_pt_dists(output, data_output, args.iov, plotdir),
        plot_roc(cache, args.iov, plotdir),
        plot_vs_pt(result, 'threshold', 'TopvsQCD threshold', 'wp_threshold_vs_pt.png',
                   args.iov, plotdir),
        plot_vs_pt(result, 'signal_eff', 'Signal efficiency', 'signal_eff_vs_pt.png',
                   args.iov, plotdir),
        plot_vs_pt(result, 'mistag_achieved', 'Achieved QCD mis-tag', 'mistag_closure_vs_pt.png',
                   args.iov, plotdir, logy=True, target_lines=True),
        plot_mistag_vs_msd(output, result, args.iov, plotdir),
    ]
    for p in plots:
        if p is not None:
            print(f"wrote {p}")

    # console summary table
    print("\nDerived working points (threshold | signal eff) per pT bin:")
    pt_edges = result['pt_bin_edges']
    header = "  WP            " + "".join(f"{pt_edges[i]:.0f}-{pt_edges[i+1]:.0f}".rjust(16)
                                          for i in range(len(pt_edges) - 1))
    print(header)
    for wp in TARGET_ORDER:
        wpd = result['working_points'][wp]
        cells = []
        for thr, se in zip(wpd['threshold'], wpd['signal_eff']):
            if thr is None:
                cells.append("—".rjust(16))
            else:
                cells.append(f"{thr:.3f}|{(se if se is not None else float('nan')):.2f}".rjust(16))
        print(f"  {wp:13s}" + "".join(cells))


if __name__ == '__main__':
    main()
