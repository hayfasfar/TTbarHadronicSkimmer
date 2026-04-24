---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: coffea_latest
    language: python
    name: python3
---

```python
%load_ext autoreload
%autoreload 2
```

```python
import os
import re
from pathlib import Path

if os.path.basename(os.getcwd()) == "plots":
    os.chdir("..")

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from coffea.util import load
from plots import hep_plot as hplot
from plots.hep_plot import CMS_COLORS
import hist

hplot.setup(era="2024")
```

```python
# ── Load MC samples ──────────────────────────────────────────────────────────
coffea_dir = Path("./outputs/dy")


def _qcd_pt_sort_key(path):
    match = re.search(r"QCD_PT-(\d+)to(\d+|Inf)", path.name)
    if match is None:
        return (float("inf"), float("inf"), path.name)
    low = int(match.group(1))
    high = float("inf") if match.group(2) == "Inf" else int(match.group(2))
    return (low, high, path.name)


def _qcd_pt_label(path):
    match = re.search(r"QCD_PT-(\d+)to(\d+|Inf)", path.name)
    if match is None:
        return path.stem
    return f"QCD {match.group(1)}-{match.group(2)}"


def _load_components(paths, label_fn):
    paths = sorted(paths, key=_qcd_pt_sort_key)
    if not paths:
        raise FileNotFoundError(
            f"No QCD pT-bin files found in {coffea_dir} matching "
            "'QCD_2024_QCD_PT-*to*_noSyst.coffea'"
        )
    return {label_fn(path): load(path) for path in paths}


qcd_pt_paths = coffea_dir.glob("QCD_2024_QCD_PT-*to*_noSyst.coffea")

samples = {
    "TTbar": {"TTbar": load(coffea_dir / "TTbar_2024inclusive_noSyst.coffea")},
    "QCD": _load_components(qcd_pt_paths, _qcd_pt_label),
    #'ZPrime4000': {"Z' (4 TeV)": load(coffea_dir / 'ZPrime4000_1_2024_noSyst.coffea')},
}

# Convenience: colors and display labels per sample
sample_style = {
    "TTbar": {"color": CMS_COLORS[0], "label": r"$t\bar{t}$", "stack": True},
    "QCD": {"color": CMS_COLORS[1], "label": "QCD", "stack": True, "cmap": "tab20"},
    "ZPrime4000": {"color": CMS_COLORS[2], "label": "Z' (4 TeV)", "stack": False},
}

print("Loaded MC components:")
for sample_name, components in samples.items():
    print(f"  {sample_name}: {', '.join(components)}")
```

```python
# ── Load data (sum over all eras) ─────────────────────────────────────────────
data_eras = ["2024C", "2024D", "2024E", "2024F", "2024G", "2024H"]
data_outputs = [load(f"./outputs/dy/data_{era}_noSyst.coffea") for era in data_eras]
```

```python
# ── Helpers ───────────────────────────────────────────────────────────────────
def _plot_axis_name(h):
    axis_names = [axis.name for axis in h.axes]
    candidates = [name for name in axis_names if name not in {"systematic", "anacat"}]
    if len(candidates) != 1:
        raise ValueError(f"Could not identify one plot axis from {axis_names}")
    return candidates[0]


def get_hist(output, var, anacat_id=None, syst="nominal"):
    """Fetch histogram `var` from a single MC output, project to 1D."""
    h = output[var][syst, ...]
    axis_name = _plot_axis_name(h)
    if anacat_id is not None:
        h = h[anacat_id, :].project(axis_name)
    else:
        h = h.project(axis_name)
    return h


def sum_hists(hists):
    """Return the bin-by-bin sum of a non-empty list of hist histograms."""
    if not hists:
        raise ValueError("Need at least one histogram to sum")
    result = hists[0]
    for h in hists[1:]:
        result = result + h
    return result


def draw_uncertainty_band(ax, h, *, label, density=False, hatch="///", zorder=2):
    """Draw the statistical uncertainty of a 1D hist as a hatched band."""
    variances = h.variances()
    if variances is None:
        return

    values = h.values()
    edges = h.axes[0].edges
    errors = np.sqrt(np.clip(variances, 0.0, None))

    if density:
        area = np.sum(values * np.diff(edges))
        if area > 0:
            values = values / area
            errors = errors / area

    lower = np.clip(values - errors, 0.0, None)
    upper = values + errors
    ax.fill_between(
        edges,
        np.r_[lower, lower[-1]],
        np.r_[upper, upper[-1]],
        step="post",
        facecolor="none",
        edgecolor="black",
        hatch=hatch,
        linewidth=0.0,
        label=label,
        zorder=zorder,
    )


def get_sample_component_hists(sample_name, var, anacat_id=None, syst="nominal"):
    """Fetch 1D histograms for every loaded component of one sample."""
    hists = []
    labels = []
    for component_label, output in samples[sample_name].items():
        hists.append(get_hist(output, var, anacat_id=anacat_id, syst=syst))
        labels.append(component_label)
    return hists, labels


def get_sample_hist(sample_name, var, anacat_id=None, syst="nominal"):
    """Fetch one sample histogram, summing components when needed."""
    hists, _ = get_sample_component_hists(
        sample_name, var, anacat_id=anacat_id, syst=syst
    )
    return sum_hists(hists)


def component_colors(sample_name, n):
    """Return one color per component, using a colormap for stacked samples."""
    style = sample_style[sample_name]
    if n == 1 or "cmap" not in style:
        return [style["color"]] * n
    cmap = plt.get_cmap(style["cmap"])
    if hasattr(cmap, "colors"):
        color_order = list(range(0, cmap.N, 2)) + list(range(1, cmap.N, 2))
        return [cmap(color_order[i % cmap.N]) for i in range(n)]
    return [cmap(0.35 + 0.5 * i / max(n - 1, 1)) for i in range(n)]


def sort_legend_entries(handles, labels):
    """Order legend entries with QCD pT bins from low to high."""
    def key(item):
        _, label = item
        qcd_match = re.match(r"QCD (\d+)-(\d+|Inf)$", label)
        if label in {r"$t\bar{t}$", "TTbar"}:
            return (0, 0)
        if qcd_match:
            return (1, int(qcd_match.group(1)))
        if label == "MC unc.":
            return (2, 0)
        if label == "Data":
            return (3, 0)
        return (4, label)

    entries = sorted(zip(handles, labels), key=key)
    return zip(*entries) if entries else ([], [])


def get_data_hist(var, anacat_id=None, syst="nominal"):
    """Fetch histogram `var` from all data eras and return their sum."""
    hists = []
    for o in data_outputs:
        h = o[var][syst, ...]
        axis_name = _plot_axis_name(h)
        if anacat_id is not None:
            h = h[anacat_id, :].project(axis_name)
        else:
            h = h.project(axis_name)
        hists.append(h)
    return sum_hists(hists)
```

```python
# Inspect available keys and categories (use any sample as reference)
ref = next(iter(next(iter(samples.values())).values()))
print("Histogram keys:", list(ref.keys()))
print("Categories:", ref["analysisCategories"])
```

## Per-sample distributions — central and forward categories

Each row: one variable.  
Left column: central (`|Δy| < 1`), right column: forward (`|Δy| > 1`).  
TTbar and QCD pT bins are drawn as one visible MC stack, with a hatched total-MC uncertainty band.

```python
plot_specs = [
    ("ttbarmass", r"$m_{t\bar{t}}$ [GeV]"),
    ("jetmsd", r"Leading jet $m_{SD}$ [GeV]"),
    ("jetmsd1", r"Subleading jet $m_{SD}$ [GeV]"),
    ("jet0_pt", r"Leading jet $p_T$ [GeV]"),
    ("jet0_eta", r"Leading jet $\eta$"),
    ("jet0_phi", r"Leading jet $\phi$"),
    ("jet0_rapidity", r"Leading jet rapidity"),
    ("jet1_pt", r"Subleading jet $p_T$ [GeV]"),
    ("jet1_eta", r"Subleading jet $\eta$"),
    ("jet1_phi", r"Subleading jet $\phi$"),
    ("jet1_rapidity", r"Subleading jet rapidity"),
    ("jetdy", r"$\Delta y$"),
    ("ht", r"$H_T$ [GeV]"),
]

# anacat IDs: 0=atcen, 1=atfwd, 2=2tcen, 3=2tfwd
cat_pairs = [
    (0, r"$|\Delta y| < 1$  (at least 1 top-tagged)"),
    (1, r"$|\Delta y| > 1$  (at least 1 top-tagged)"),
]

plot_density = False
stack_sample_names = [
    name for name in samples if sample_style[name].get("stack", False)
]
overlay_sample_names = [
    name for name in samples if not sample_style[name].get("stack", False)
]

for var, xlabel in plot_specs:
    fig, axes = plt.subplots(1, len(cat_pairs), figsize=(10 * len(cat_pairs), 8))

    for ax, (cat_id, cat_label) in zip(axes, cat_pairs):
        # ── Stacked MC components, e.g. TTbar and QCD pT bins ──
        stack_hists = []
        stack_labels = []
        stack_colors = []
        for sample_name in stack_sample_names:
            style = sample_style[sample_name]
            try:
                hists, labels = get_sample_component_hists(
                    sample_name, var, anacat_id=cat_id
                )
            except Exception:
                continue
            if len(labels) == 1:
                labels = [style["label"]]
            stack_hists.extend(hists)
            stack_labels.extend(labels)
            stack_colors.extend(component_colors(sample_name, len(hists)))

        if stack_hists:
            h_stack_total = sum_hists(stack_hists)
            hep.histplot(
                stack_hists,
                ax=ax,
                histtype="fill",
                stack=True,
                color=stack_colors,
                edgecolor="black",
                linewidth=0.4,
                label=stack_labels,
                density=plot_density,
                alpha=0.85,
                zorder=1,
            )
            draw_uncertainty_band(
                ax,
                h_stack_total,
                label="MC unc.",
                density=plot_density,
                zorder=2,
            )

        # ── Non-stacked MC overlays, e.g. signal ──
        for sample_name in overlay_sample_names:
            style = sample_style[sample_name]
            try:
                h = get_sample_hist(sample_name, var, anacat_id=cat_id)
            except Exception:
                continue
            hep.histplot(
                h,
                ax=ax,
                histtype="step",
                color=style["color"],
                label=style["label"],
                density=plot_density,
                linewidth=2.0,
                zorder=3,
            )

        # ── Data (sum of all eras) ──
        try:
            h_data = get_data_hist(var, anacat_id=cat_id)
            hep.histplot(
                h_data,
                ax=ax,
                histtype="errorbar",
                color="black",
                label="Data",
                density=plot_density,
                zorder=4,
            )
        except Exception:
            pass

        hplot.quick_label(xlabel=xlabel, data=True, ax=ax)
        ax.text(
            0.97,
            0.97,
            cat_label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=20,
        )
        ax.set_ylabel("A.U." if plot_density else "# Events")
        ax.set_xlabel(xlabel, labelpad=20)
        handles, labels = sort_legend_entries(*ax.get_legend_handles_labels())
        ax.legend(handles, labels, ncol=2, fontsize=14)

    plt.tight_layout()
    plt.show()
```

## Gen-level distributions (TTbar and ZPrime only)

```python
gen_specs = [
    ("gen_mt", r"Gen top mass [GeV]"),
    ("gen_mttbar", r"Gen $m_{t\bar{t}}$ [GeV]"),
    ("jet0_gen_dr", r"Leading jet $\Delta R$ (gen)"),
    ("jet1_gen_dr", r"Subleading jet $\Delta R$ (gen)"),
]

gen_samples = {k: v for k, v in samples.items() if k != "QCD"}

for var, xlabel in gen_specs:
    fig, axes = plt.subplots(1, len(cat_pairs), figsize=(10 * len(cat_pairs), 8))

    for ax, (cat_id, cat_label) in zip(axes, cat_pairs):
        for sample_name in gen_samples:
            style = sample_style[sample_name]
            try:
                h = get_sample_hist(sample_name, var, anacat_id=cat_id)
            except Exception:
                continue
            hep.histplot(
                h, ax=ax, histtype="step", color=style["color"], label=style["label"]
            )

        hplot.quick_label(xlabel=xlabel, data=False, ax=ax)
        ax.text(
            0.97,
            0.97,
            cat_label,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=13,
        )
        ax.set_ylabel("# Events (weighted)")
        ax.set_xlabel(xlabel, labelpad=20)
        ax.legend()

    plt.tight_layout()
    plt.show()
```

## Cutflow

```python
for sample_name, output in samples.items():
    print(f"── {sample_name} ──")
    for component_label, component_output in output.items():
        print(f"  {component_label}")
        print(component_output["cutflow"])
        print()
    print()
```

```python

```
