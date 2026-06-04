from __future__ import annotations

import math
from collections.abc import Mapping


BASE_CUTFLOW_STEPS = [
    {
        "key": "input_events",
        "label": "Input events",
        "group": "Bookkeeping",
        "description": "Raw entries in the input chunk before the processor pre-filter.",
    },
    {
        "key": "analysis_events",
        "label": "Events entering nominal analysis",
        "group": "Bookkeeping",
        "description": "Events after processor pre-filters, lumi mask, and blinding if enabled.",
    },
    {
        "key": "trigger",
        "label": "Trigger",
        "group": "Preselection",
        "description": "Cumulative HLT requirement.",
    },
    {
        "key": "htCut",
        "label": "HT > threshold",
        "group": "Preselection",
        "description": "Cumulative AK4 HT requirement using the configured processor threshold.",
    },
    {
        "key": "metfilter",
        "label": "MET filters",
        "group": "Preselection",
        "description": "Cumulative recommended event filter requirement for the IOV.",
    },
    {
        "key": "jetkincut",
        "label": "At least one AK8 passing kinematics",
        "group": "Preselection",
        "description": "Cumulative AK8 pT and rapidity requirement.",
    },
    {
        "key": "twoFatJets",
        "label": "At least two selected AK8 jets",
        "group": "Preselection",
        "description": "Cumulative requirement after filtering AK8 jets by pT and rapidity.",
    },
    {
        "key": "preselection",
        "label": "All preselection",
        "group": "Preselection",
        "description": "Combined trigger, HT, MET filter, AK8 kinematics, and two-AK8 requirement.",
    },
    {
        "key": "dphi",
        "label": "Delta phi topology",
        "group": "TTbar candidate",
        "denominator_key": "preselection",
        "description": "Cumulative back-to-back requirement on the two leading AK8 candidates.",
    },
    {
        "key": "ttbarcand",
        "label": "TTbar candidate",
        "group": "TTbar candidate",
        "denominator_key": "dphi",
        "description": "Cumulative delta-phi and subjet requirements.",
    },
    {
        "key": "tag_jet0",
        "label": "Leading top-tagged jet",
        "group": "Tagging",
        "denominator_key": "ttbarcand",
        "description": "Events where the higher-score AK8 candidate passes the configured top-tag WP.",
    },
    {
        "key": "tag_2tag",
        "label": "Two top-tagged jets",
        "group": "Tagging",
        "denominator_key": "tag_jet0",
        "description": "Pass-region requirement used by 2DAlphabet categories.",
    },
    {
        "key": "antitag",
        "label": "Antitag region",
        "group": "Tagging",
        "denominator_key": "tag_jet0",
        "description": "Fail-region requirement used by 2DAlphabet categories.",
    },
]

LEGACY_CUTFLOW_ALIASES = {
    "input_events": ("all events 1",),
    "analysis_events": ("all events",),
    "preselection": ("after_eventCut",),
    "ttbarcand": ("after_ttbarcandCuts",),
}

_KNOWN_LEGACY_KEYS = {
    "all events 1",
    "all events",
    "sumw",
    "sumw2",
    "after_lumimask",
    "trigger",
    "htCut",
    "metfilter",
    "jetkincut",
    "twoFatJets",
    "after_eventCut",
    "after_ttbarcandCuts",
}


def cutflow_step_metadata(anacats=None):
    """Return ordered cutflow step metadata, including configured categories."""

    steps = [dict(step) for step in BASE_CUTFLOW_STEPS]
    for label in anacats or []:
        denominator_key = "ttbarcand"
        fallback_denominator_key = None
        if str(label).startswith("at"):
            denominator_key = "antitag"
            fallback_denominator_key = "ttbarcand"
        elif str(label).startswith("2t"):
            denominator_key = "tag_2tag"
            fallback_denominator_key = "ttbarcand"
        steps.append(
            {
                "key": f"category_{label}",
                "label": f"Category {label}",
                "group": "Analysis category",
                "denominator_key": denominator_key,
                "fallback_denominator_key": fallback_denominator_key,
                "description": "Configured analysis category after candidate and top-tag selections.",
            }
        )
    return steps


def _plain_mapping(value):
    if value is None:
        return {}
    if isinstance(value, Mapping):
        return dict(value)
    return value


def _lookup(mapping, key, default=0.0):
    mapping = _plain_mapping(mapping)
    if not mapping:
        return default
    try:
        return mapping.get(key, default)
    except AttributeError:
        return default


def _legacy_category_labels(cutflow):
    cutflow = _plain_mapping(cutflow)
    labels = []
    for key in cutflow:
        if key in _KNOWN_LEGACY_KEYS:
            continue
        if key.startswith("after_"):
            continue
        labels.append(str(key))
    return labels


def _steps_for_output(output):
    if output.get("cutflow_table_steps"):
        return output["cutflow_table_steps"]
    return cutflow_step_metadata(_legacy_category_labels(output.get("cutflow", {})))


def _legacy_cutflow_to_table_counts(cutflow):
    cutflow = _plain_mapping(cutflow)
    counts = {}
    for new_key, old_keys in LEGACY_CUTFLOW_ALIASES.items():
        for old_key in old_keys:
            if old_key in cutflow:
                counts[new_key] = cutflow[old_key]
                break

    for key in ("trigger", "htCut", "metfilter", "jetkincut", "twoFatJets"):
        if key in cutflow:
            counts[key] = cutflow[key]

    for label in _legacy_category_labels(cutflow):
        counts[f"category_{label}"] = cutflow[label]
    return counts


def _cutflow_counts_for_output(output):
    if output.get("cutflow_unweighted"):
        return output["cutflow_unweighted"]
    return _legacy_cutflow_to_table_counts(output.get("cutflow", {}))


def _as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _percent(numerator, denominator):
    numerator = _as_float(numerator)
    denominator = _as_float(denominator)
    if denominator == 0.0:
        return None
    return 100.0 * numerator / denominator


def _format_number(value, precision=3):
    if value is None:
        return "-"
    value = _as_float(value)
    if math.isclose(value, round(value), rel_tol=0.0, abs_tol=1e-9):
        return f"{int(round(value)):,}"
    return f"{value:,.{precision}f}"


def _format_percent(value):
    if value is None:
        return "-"
    return f"{value:.2f}"


def _normalization_text(output):
    normalization = output.get("normalization", {}) if isinstance(output, Mapping) else {}
    if not normalization:
        return None
    if normalization.get("applied"):
        scale = normalization.get("scale_factor", 1.0)
        lumi = normalization.get("lumi_pb")
        xsec = normalization.get("xsec_pb")
        pieces = [f"scaled by {scale:.6g}"]
        if lumi is not None:
            pieces.append(f"lumi = {lumi:g} pb^-1")
        if xsec is not None:
            pieces.append(f"xsec = {xsec:g} pb")
        return ", ".join(pieces)
    reason = normalization.get("reason")
    if reason:
        return f"not scaled ({reason})"
    return "not scaled"


def choose_weight_mapping(output, weight_mode="scaled"):
    """Return the requested weight mapping and the table column label."""

    if weight_mode == "none":
        return None, None
    if weight_mode == "scaled" and "cutflow_weighted_scaled" in output:
        return output["cutflow_weighted_scaled"], "Yield"
    if weight_mode == "scaled" and "cutflow_scaled" in output:
        return _legacy_cutflow_to_table_counts(output["cutflow_scaled"]), "Yield"
    if weight_mode == "scaled":
        return output.get("cutflow_weighted", {}), "Sumw"
    if weight_mode == "raw":
        return output.get("cutflow_weighted", {}), "Sumw"
    raise ValueError("weight_mode must be one of: scaled, raw, none")


def build_cutflow_rows(output, weight_mode="scaled", include_zero=False):
    """Build ordered table rows from a processor output dictionary.

    The processor writes explicit ``cutflow_unweighted`` and ``cutflow_weighted``
    maps. For older outputs, this falls back to the legacy flat ``cutflow`` map.
    """

    output = output or {}
    steps = _steps_for_output(output)
    counts = _cutflow_counts_for_output(output)
    weights, weight_label = choose_weight_mapping(output, weight_mode)

    rows = []
    previous_count = None
    first_count = None
    row_counts = {}
    for step in steps:
        key = step["key"]
        count = _as_float(_lookup(counts, key, 0.0))
        raw_weighted = _lookup(weights, key, None) if weights is not None else None
        weighted = _as_float(raw_weighted) if raw_weighted is not None else None
        if not include_zero and count == 0.0 and (weighted is None or weighted == 0.0):
            continue

        if first_count is None and count != 0.0:
            first_count = count

        denominator = previous_count
        denominator_key = step.get("denominator_key")
        if denominator_key:
            denominator = row_counts.get(denominator_key)
            if not denominator:
                fallback_key = step.get("fallback_denominator_key")
                denominator = row_counts.get(fallback_key) if fallback_key else previous_count

        rows.append(
            {
                "key": key,
                "group": step.get("group", ""),
                "step": step.get("label", key),
                "events": count,
                "weighted": weighted,
                "eff_prev_percent": _percent(count, denominator),
                "eff_total_percent": _percent(count, first_count),
                "weight_label": weight_label,
            }
        )
        row_counts[key] = count
        previous_count = count

    return rows


def format_markdown_table(rows, title=None, normalization=None):
    if not rows:
        return "No cutflow rows found."

    weight_label = rows[0].get("weight_label")
    headers = ["Group", "Step", "Events"]
    if weight_label:
        headers.append(weight_label)
    headers.extend(["Eff. prev [%]", "Eff. total [%]"])

    lines = []
    if title:
        lines.extend([f"## {title}", ""])
    if normalization:
        lines.extend([f"_Normalization: {normalization}_", ""])

    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        cells = [
            row["group"],
            row["step"],
            _format_number(row["events"], precision=0),
        ]
        if weight_label:
            cells.append(_format_number(row["weighted"]))
        cells.extend([
            _format_percent(row["eff_prev_percent"]),
            _format_percent(row["eff_total_percent"]),
        ])
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _latex_escape(value):
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    text = str(value)
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def format_latex_table(rows, caption=None, label=None, normalization=None):
    if not rows:
        return "% No cutflow rows found."

    weight_label = rows[0].get("weight_label")
    headers = ["Group", "Step", "Events"]
    if weight_label:
        headers.append(weight_label)
    headers.extend([r"Eff. prev [\%]", r"Eff. total [\%]"])
    alignment = "llr" + ("r" if weight_label else "") + "rr"

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\small",
    ]
    if caption:
        lines.append(rf"\caption{{{_latex_escape(caption)}}}")
    if label:
        lines.append(rf"\label{{{_latex_escape(label)}}}")
    if normalization:
        lines.append(rf"\par\medskip\footnotesize Normalization: {_latex_escape(normalization)}\par\medskip")
    lines.extend([
        rf"\begin{{tabular}}{{{alignment}}}",
        r"\hline",
        " & ".join(headers) + r" \\",
        r"\hline",
    ])
    for row in rows:
        cells = [
            _latex_escape(row["group"]),
            _latex_escape(row["step"]),
            _format_number(row["events"], precision=0),
        ]
        if weight_label:
            cells.append(_format_number(row["weighted"]))
        cells.extend([
            _format_percent(row["eff_prev_percent"]),
            _format_percent(row["eff_total_percent"]),
        ])
        lines.append(" & ".join(cells) + r" \\")
    lines.extend([
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
    ])
    return "\n".join(lines)


def cutflow_normalization_text(output):
    return _normalization_text(output)
