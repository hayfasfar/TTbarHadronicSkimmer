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
    display_name: coffea-dask (3.14.4)
    language: python
    name: python3
---

# ttbaranalysis

Use the configuration widgets below to choose the datasets, year, taggers, thresholds, systematics, and execution mode. The widget state is saved to `.last_config.json`, then converted into the `args` object used by the same analysis flow as `ttbaranalysis.py`.

```python
from coffea import util
from coffea.nanoevents import NanoAODSchema, BaseSchema
import coffea.processor as processor

import itertools
import time
from datetime import date
import json
import os
from types import SimpleNamespace

from dask.distributed import Client, performance_report

import warnings

warnings.filterwarnings("ignore")
import logging
import dask

dask.config.set({"logging.distributed": "error"})

for name in [
    "distributed",
    "distributed.scheduler",
    "distributed.core",
    "distributed.nanny",
    "distributed.worker",
]:
    logging.getLogger(name).setLevel(logging.CRITICAL)

default_datastets = ["data", "TTbar", "QCD"]
default_signals = ["RSGluon", "ZPrime10", "ZPrime30", "ZPrimeDM", "ZPrime1"]

from ttbarprocessor import TTbarResProcessor
from python.functions import printTime, makeSaveDirectories
```

```python
%load_ext autoreload
%autoreload 2
```

```python
# ── Widgets for interactive configuration ─────────────────────────────────────
import ipywidgets as widgets
from IPython.display import clear_output, display
import json, os

CONFIG_FILE = ".last_config.json"

for _old_widget in list(globals().get("WIDGETS", {}).values()) + [
    globals().get("btn_reset"),
]:
    if _old_widget is not None:
        try:
            _old_widget.close()
        except Exception:
            pass
clear_output(wait=True)

DEFAULTS = dict(
    dataset=["ZPrimeLocal"],
    signals=False,
    iov="2024",
    subsample=[],
    mass="",
    blind=False,
    bkgest=None,
    toptagger="deepak8",
    redirector="rootfiles/",
    ttagWP="medium",
    btagger="deepcsv",
    ht="1500",
    noSyst=False,
    ntuple=False,
    overwrite=False,
    dask=False,
    env="lpc",
    test=False,
    nocluster=False,
)


def load_config():
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE) as f:
                return {**DEFAULTS, **json.load(f)}
        except Exception:
            pass
    return dict(DEFAULTS)


cfg = load_config()

style = {"description_width": "80px"}
layout = widgets.Layout(width="210px")
layout_wide = widgets.Layout(width="260px")

_dataset_opts = [
    "data",
    "QCD",
    "QCD_flat",
    "TTbar",
    "ZPrime1",
    "ZPrime10",
    "ZPrime30",
    "ZPrimeDM",
    "RSGluon",
    "ZPrimeLocal",
]
_redirector_opts = [
    ("Local (rootfiles/)", "rootfiles/"),
    ("FNAL XRootD (root://cmsxrootd.fnal.gov/)", "root://cmsxrootd.fnal.gov/"),
    ("CMS xcache (root://xcache/)", "root://xcache/"),
    ("Winterfell (/mnt/data/cms/)", "/mnt/data/cms/"),
]
_redirector_vals = [v for _, v in _redirector_opts]
_env_opts = ["casa", "lpc", "winterfell", "local"]
_iov_opts = ["2022", "2023", "2024"]
_bkgest_opts = [("None", None), "2dalphabet", "mistag"]
_toptagger_opts = ["deepak8", "cmsv2"]
_ttagWP_opts = ["loose", "medium", "tight"]
_btagger_opts = ["deepcsv", "csvv2"]
_ht_opts = ["1400", "950"]
_manifest_files = {
    "data": "data/nanoAOD/data.json",
    "QCD": "data/nanoAOD/QCD.json",
    "QCD_flat": "data/nanoAOD/QCD_flat.json",
    "TTbar": "data/nanoAOD/TTbar.json",
    "ZPrime1": "data/nanoAOD/ZPrime1.json",
    "ZPrime10": "data/nanoAOD/ZPrime10.json",
    "ZPrime30": "data/nanoAOD/ZPrime30.json",
    "ZPrimeDM": "data/nanoAOD/ZPrimeDM.json",
    "RSGluon": "data/nanoAOD/RSGluon.json",
    "ZPrimeLocal": "data/nanoAOD/local_xsec_test.json",
}


def _option_values(options):
    return [option[1] if isinstance(option, tuple) else option for option in options]


def _valid_choice(value, options, default):
    values = _option_values(options)
    return value if value in values else default


def _valid_multi(values, options, default=()):
    allowed = set(_option_values(options))
    selected = [value for value in (values or []) if value in allowed]
    if selected:
        return tuple(selected)
    return tuple(value for value in default if value in allowed)


def _manifest_subsections(dataset, iov):
    path = _manifest_files.get(dataset)
    if not path or not os.path.exists(path):
        return []

    try:
        with open(path) as f:
            manifest = json.load(f)
    except Exception:
        return []

    entry = manifest.get(iov)
    if isinstance(entry, dict) and "files" not in entry:
        return list(entry.keys())
    return []


def _available_subsamples(datasets, iov):
    subsamples = []
    seen = set()
    for dataset in datasets:
        for subsection in _manifest_subsections(dataset, iov):
            if subsection not in seen:
                seen.add(subsection)
                subsamples.append(subsection)
    return subsamples


_initial_datasets = _valid_multi(cfg["dataset"], _dataset_opts, DEFAULTS["dataset"])
_initial_iov = _valid_choice(cfg["iov"], _iov_opts, DEFAULTS["iov"])
_initial_subsample_datasets = (
    list(default_signals) if cfg["signals"] else list(_initial_datasets)
)
_initial_subsample_options = _available_subsamples(
    _initial_subsample_datasets, _initial_iov
)
_initial_subsamples = _valid_multi(cfg["subsample"], _initial_subsample_options)


# ── Widget definitions ─────────────────────────────────────────────────────────
w_dataset = widgets.SelectMultiple(
    options=_dataset_opts,
    value=_initial_datasets,
    description="Dataset",
    style=style,
    layout=widgets.Layout(width="210px", height="150px"),
)
w_signals = widgets.Checkbox(
    value=cfg["signals"], description="Signals only", style=style, layout=layout
)
w_iov = widgets.Dropdown(
    options=_iov_opts,
    value=_initial_iov,
    description="IOV",
    style=style,
    layout=layout,
)
w_subsample = widgets.SelectMultiple(
    options=_initial_subsample_options,
    value=_initial_subsamples,
    description="Subsample",
    style=style,
    layout=widgets.Layout(width="260px", height="140px"),
)
w_mass = widgets.Text(
    value=cfg["mass"],
    placeholder="e.g. 1000,2000",
    description="Mass pts",
    style=style,
    layout=layout,
)
w_blind = widgets.Checkbox(
    value=cfg["blind"], description="Blind", style=style, layout=layout
)
w_bkgest = widgets.Dropdown(
    options=_bkgest_opts,
    value=_valid_choice(cfg["bkgest"], _bkgest_opts, DEFAULTS["bkgest"]),
    description="Bkg est",
    style=style,
    layout=layout,
)
w_toptagger = widgets.Dropdown(
    options=_toptagger_opts,
    value=_valid_choice(cfg["toptagger"], _toptagger_opts, DEFAULTS["toptagger"]),
    description="Top tagger",
    style=style,
    layout=layout,
)
w_redirector = widgets.Dropdown(
    options=_redirector_opts,
    value=cfg["redirector"] if cfg["redirector"] in _redirector_vals else "rootfiles/",
    description="Redirector",
    style=style,
    layout=layout_wide,
)
w_ttagWP = widgets.Dropdown(
    options=_ttagWP_opts,
    value=_valid_choice(cfg["ttagWP"], _ttagWP_opts, DEFAULTS["ttagWP"]),
    description="ttag WP",
    style=style,
    layout=layout,
)
w_btagger = widgets.Dropdown(
    options=_btagger_opts,
    value=_valid_choice(cfg["btagger"], _btagger_opts, DEFAULTS["btagger"]),
    description="B tagger",
    style=style,
    layout=layout,
)
w_ht = widgets.Dropdown(
    options=_ht_opts,
    value=_valid_choice(cfg["ht"], _ht_opts, DEFAULTS["ht"]),
    description="HT cut",
    style=style,
    layout=layout,
)
w_noSyst = widgets.Checkbox(
    value=cfg["noSyst"], description="No syst", style=style, layout=layout
)
w_ntuple = widgets.Checkbox(
    value=cfg["ntuple"], description="Ntuple", style=style, layout=layout
)
w_overwrite = widgets.Checkbox(
    value=cfg["overwrite"], description="Overwrite", style=style, layout=layout
)
w_dask = widgets.Checkbox(
    value=cfg["dask"], description="Dask", style=style, layout=layout
)
w_env = widgets.Dropdown(
    options=_env_opts,
    value=cfg["env"] if cfg["env"] in _env_opts else "lpc",
    description="Env",
    style=style,
    layout=layout,
)
w_test = widgets.Checkbox(
    value=cfg["test"], description="Test", style=style, layout=layout
)
w_nocluster = widgets.Checkbox(
    value=cfg["nocluster"], description="No cluster", style=style, layout=layout
)

# ── Central widget registry ────────────────────────────────────────────────────
# To add a new config field: add it to DEFAULTS above and WIDGETS below.
# save_config, build_args, and reset_to_defaults all derive from this dict.
WIDGETS = {
    "dataset": w_dataset,
    "signals": w_signals,
    "iov": w_iov,
    "subsample": w_subsample,
    "mass": w_mass,
    "blind": w_blind,
    "bkgest": w_bkgest,
    "toptagger": w_toptagger,
    "redirector": w_redirector,
    "ttagWP": w_ttagWP,
    "btagger": w_btagger,
    "ht": w_ht,
    "noSyst": w_noSyst,
    "ntuple": w_ntuple,
    "overwrite": w_overwrite,
    "dask": w_dask,
    "env": w_env,
    "test": w_test,
    "nocluster": w_nocluster,
}

_MULTI = widgets.SelectMultiple


def _widget_value(w):
    return list(w.value) if isinstance(w, _MULTI) else w.value


def save_config(_=None):
    cfg = {k: _widget_value(w) for k, w in WIDGETS.items()}
    with open(CONFIG_FILE, "w") as f:
        json.dump(cfg, f, indent=2)


def reset_to_defaults(_):
    for key, w in WIDGETS.items():
        default = DEFAULTS[key]
        w.value = tuple(default) if isinstance(w, _MULTI) else default


def refresh_subsample_options(_=None):
    selected_datasets = (
        list(default_signals) if w_signals.value else list(w_dataset.value)
    )
    options = _available_subsamples(selected_datasets, w_iov.value)
    current = [v for v in w_subsample.value if v in options]
    w_subsample.options = options
    w_subsample.value = tuple(current)


for w in WIDGETS.values():
    w.observe(save_config, names="value")

for w in (w_dataset, w_iov, w_signals):
    w.observe(refresh_subsample_options, names="value")

refresh_subsample_options()

btn_reset = widgets.Button(
    description="↺ Reset to Defaults",
    button_style="warning",
    layout=widgets.Layout(width="160px", margin="8px 0 0 0"),
)
btn_reset.on_click(reset_to_defaults)

# ── Display ───────────────────────────────────────────────────────────────────
_loaded = (
    "restored from last session" if os.path.exists(CONFIG_FILE) else "using defaults"
)
print(f"Config {_loaded}; changes save to {CONFIG_FILE}.")
print("Datasets")
display(w_dataset)
display(w_signals)
display(w_iov)
print("Subsections")
display(w_subsample)
display(w_mass)
print("Analysis options")
for _widget in (
    w_blind,
    w_bkgest,
    w_toptagger,
    w_redirector,
    w_ttagWP,
    w_btagger,
    w_ht,
    w_noSyst,
    w_ntuple,
    w_overwrite,
):
    display(_widget)
print("Run options")
for _widget in (w_dask, w_env, w_test, w_nocluster, btn_reset):
    display(_widget)
print("Adjust widgets above, then run the next cell to apply settings.")
```

```python
from types import SimpleNamespace


def build_args():
    selected_datasets = list(w_dataset.value)
    if w_signals.value:
        selected_datasets = list(default_signals)

    raw_mass = w_mass.value.strip()
    mass_list = []
    if raw_mass:
        parts = [m.strip() for m in raw_mass.split(",")]
        invalid = [p for p in parts if not p.isdigit()]
        if invalid:
            print(f"Warning: invalid mass entries ignored: {invalid}")
        mass_list = [p for p in parts if p.isdigit()]

    cfg = {k: _widget_value(w) for k, w in WIDGETS.items()}
    cfg["dataset"] = selected_datasets
    cfg["era"] = []
    cfg["pt"] = []
    cfg["mass"] = mass_list
    return SimpleNamespace(**cfg)


args = build_args()
print("------args------")
for argname, value in vars(args).items():
    print(argname, "=", value)
print("----------------")
```

```python
import subprocess
import traceback


def _build_sample_metadata(sample, subsection, iov, metadata):
    sample_metadata = {
        "sample": sample,
        "subsample": subsection or sample,
        "year": iov,
        "is_mc": not (("data" in sample.lower()) or ("singlemu" in sample.lower())),
    }
    sample_metadata.update(metadata)
    return sample_metadata


def _output_subsection(sample, subsection):
    if sample == "QCD" and subsection and subsection.startswith("QCD_"):
        return subsection.removeprefix("QCD_")
    return subsection


def _parse_manifest_entry(sample, subsection, iov, entry):
    if isinstance(entry, dict) and "files" in entry:
        files = entry["files"]
        metadata = dict(entry.get("metadata", {}))
    else:
        files = entry
        metadata = {}

    return list(files), _build_sample_metadata(sample, subsection, iov, metadata)


def _collect_manifest_sections(sample, iov, manifest, subsections):
    iov_entry = manifest[iov]

    if isinstance(iov_entry, dict) and "files" not in iov_entry:
        requested_sections = subsections if subsections else list(iov_entry.keys())
        entries = []
        for subsection in requested_sections:
            if subsection not in iov_entry:
                print(f"{subsection} not in {sample} {iov}")
                continue
            files, metadata = _parse_manifest_entry(
                sample, subsection, iov, iov_entry[subsection]
            )
            entries.append((subsection, files, metadata))
        return entries

    files, metadata = _parse_manifest_entry(sample, "", iov, iov_entry)
    return [("", files, metadata)]


def _format_section_label(iov, sample, subsection):
    return f"{iov} {sample} {subsection}".strip()


def _print_runner_block(lines, rule_char="-", width=66):
    print(rule_char * width)
    for line in lines:
        print(line)
    print(rule_char * width)


def _archive_existing_output(path, tag="old"):
    if not os.path.exists(path):
        return None

    root, ext = os.path.splitext(path)
    archive_path = f"{root}_{tag}{ext}"

    os.replace(path, archive_path)
    return archive_path


def _close_dask_resources(client, cluster):
    if client is not None:
        client.close()
    if cluster is not None:
        cluster.close()
    return None, None


def _start_dask_resources(args, repo_root, upload_to_dask, dask_memory, nworkers):
    client = None
    cluster = None

    if not args.dask:
        return client, cluster

    if args.env == "lpc":
        if not args.nocluster:
            cluster = LPCCondorCluster(
                memory=dask_memory,
                transfer_input_files=upload_to_dask,
                scheduler_options={"dashboard_address": ":8787"},
            )
            cluster.adapt(minimum=1, maximum=100)
    elif args.env == "casa":
        if not args.nocluster:
            from coffea_casa import CoffeaCasaCluster

            cluster = CoffeaCasaCluster(memory=dask_memory)
            cluster.adapt(minimum=4, maximum=400)
    else:
        cluster = dask.distributed.LocalCluster(
            n_workers=nworkers,
            threads_per_worker=1,
            scheduler_port=0,
            dashboard_address=":8787",
        )

    client = Client(cluster)

    if args.env == "casa" and not args.nocluster:
        from distributed.diagnostics.plugin import UploadDirectory

        client.register_worker_plugin(
            UploadDirectory(
                os.path.join(repo_root, "data"), restart=True, update_path=True
            ),
            nanny=True,
        )
        client.register_worker_plugin(
            UploadDirectory(
                os.path.join(repo_root, "python"), restart=True, update_path=True
            ),
            nanny=True,
        )
        client.upload_file(os.path.join(repo_root, "ttbarprocessor.py"))

    return client, cluster


def run_analysis(args):
    tic = time.time()

    savedir = f"outputs/dy/"

    if args.dask and args.env == "lpc":
        from lpcjobqueue import LPCCondorCluster

    samples = args.dataset
    IOV = args.iov
    useDeepAK8 = args.toptagger == "deepak8"
    useDeepCSV = args.btagger == "deepcsv"
    htCut = 1400.0 if args.ht == "1400" else 950.0
    dask_memory = "5GB"
    chunksize_dask = 100000
    chunksize_futures = 200000
    maxchunks = 10 if args.test else None

    systematics = [
        "nominal",
        "jes",
        "jer",
        "pileup",
        "pdf",
        "q2",
        "ttag_pt1",
        #'ttag_pt2',
        #'ttag_pt3'
    ]

    if ("2016" in IOV) or ("2017" in IOV):
        systematics.append("prefiring")

    if args.bkgest == "2dalphabet":
        systematics.append("transferFunction")

    ttagcats = ["at", "2t"]
    ycats = ["cen", "fwd"]

    anacats = [t + y for t, y in itertools.product(ttagcats, ycats)]
    label_map = {i: label for i, label in enumerate(anacats)}

    with open("out.log", "w") as f:
        print("\n" + date.today().isoformat(), file=f)
        print("categories =", label_map, file=f)
        print("\n", file=f)
        if not args.noSyst:
            print("systematics =", systematics, file=f)

    print("\n------args------")
    for argname, value in vars(args).items():
        print(argname, "=", value)
    if not args.noSyst:
        print("systematics =", systematics)
    print("----------------\n")

    redirector = args.redirector

    jsonfiles = {
        "data": "data/nanoAOD/data.json",
        "QCD": "data/nanoAOD/QCD.json",
        "QCD_flat": "data/nanoAOD/QCD_flat.json",
        "TTbar": "data/nanoAOD/TTbar.json",
        "ZPrime1": "data/nanoAOD/ZPrime1.json",
        "ZPrime10": "data/nanoAOD/ZPrime10.json",
        "ZPrime30": "data/nanoAOD/ZPrime30.json",
        "ZPrimeDM": "data/nanoAOD/ZPrimeDM.json",
        "RSGluon": "data/nanoAOD/RSGluon.json",
        "ZPrimeLocal": "data/nanoAOD/local_xsec_test.json",
    }

    repo_root = os.path.abspath(os.getcwd())
    upload_to_dask = ["data", "python", "ttbarprocessor.py"]

    if not os.path.exists(savedir):
        os.makedirs(savedir)
        os.makedirs(savedir + "logs/")
        os.makedirs(savedir + "scale/")
        os.makedirs(savedir + "twodalphabet/")
        subprocess.run(
            [
                "cp",
                "ttbarprocessor.py",
                savedir
                + "logs/ttbarprocessor_"
                + date.today().isoformat().replace("-", "")
                + ".py",
            ],
            check=True,
        )
        subprocess.run(
            f"cat out.log >> {savedir}logs/ttbarprocessor_diff.txt",
            shell=True,
            check=True,
        )
    else:
        for f in os.listdir(savedir + "logs/"):
            if "ttbarprocessor" in f and "py" in f:
                subprocess.run(
                    f"cat out.log >> {savedir}logs/ttbarprocessor_diff.txt",
                    shell=True,
                    check=True,
                )
                diff_result = subprocess.run(
                    ["diff", "ttbarprocessor.py", savedir + "logs/" + f],
                    capture_output=True,
                    text=True,
                )
                with open(savedir + "logs/ttbarprocessor_diff.txt", "a") as df:
                    df.write(diff_result.stdout)

        if not os.path.exists(savedir + "scale/"):
            os.makedirs(savedir + "scale/")
        if not os.path.exists(savedir + "twodalphabet/"):
            os.makedirs(savedir + "twodalphabet/")

    makeSaveDirectories(coffea_dir=savedir)

    output = None
    metrics = None
    savefilenames = []
    skipped_outputs = []
    failures = []
    nworkers = 1 if args.test else 4

    # ── Dask cluster/client: created once and reused across all samples ────────
    client = None
    cluster = None
    client, cluster = _start_dask_resources(
        args=args,
        repo_root=repo_root,
        upload_to_dask=upload_to_dask,
        dask_memory=dask_memory,
        nworkers=nworkers,
    )

    for sample_index, sample in enumerate(samples):
        skipbadfiles = False
        inputfile = jsonfiles[sample]

        with open(inputfile) as json_file:
            subsections = (
                args.era + args.mass + args.pt + getattr(args, "subsample", [])
            )
            manifest = json.load(json_file)
            sections = _collect_manifest_sections(
                sample=sample,
                iov=IOV,
                manifest=manifest,
                subsections=subsections,
            )

            for section_index, (subsection, files, sample_metadata) in enumerate(
                sections
            ):
                files = [redirector + f for f in files]
                if args.test:
                    files = [files[int(len(files) / 2)]]
                    maxchunks = 1

                fileset = {
                    sample: {
                        "files": files,
                        "metadata": sample_metadata,
                    }
                }

                print(files[0])

                output_subsection = _output_subsection(sample, subsection)
                subString = f"_{output_subsection}" if output_subsection else ""
                if args.bkgest:
                    subString += "_bkgest"

                if (args.toptagger == "cmsv2") and (args.btagger == "csvv2"):
                    savedir = "outputs/oldanalysis/"

                savefilename = f"{savedir}{sample}_{IOV}{subString}.coffea"
                if "RSGluon" in sample:
                    subString = subString.replace(output_subsection, "")
                    savefilename = (
                        f"{savedir}{sample}{subsection}_{IOV}{subString}.coffea"
                    )
                elif "ZPrime" in sample:
                    subString = subString.replace(output_subsection, "")
                    savefilename = f'{savedir}ZPrime{subsection}_{sample.replace("ZPrime", "")}_{IOV}{subString}.coffea'
                print(f"running {IOV} {sample} {subsection}")

                if args.toptagger == "cmsv2":
                    savefilename = savefilename.replace(".coffea", "_cmsv2.coffea")
                if args.btagger == "csvv2":
                    savefilename = savefilename.replace(".coffea", "_csvv2.coffea")
                if args.ht == "950":
                    savefilename = savefilename.replace(".coffea", "_ht950.coffea")
                if args.blind:
                    savefilename = savefilename.replace(".coffea", "_blind.coffea")
                if args.noSyst:
                    savefilename = savefilename.replace(".coffea", "_noSyst.coffea")
                if args.test:
                    savefilename = savefilename.replace(".coffea", "_test.coffea")

                section_label = _format_section_label(IOV, sample, subsection)
                if section_index + 1 < len(sections):
                    next_label = _format_section_label(
                        IOV, sample, sections[section_index + 1][0]
                    )
                elif sample_index + 1 < len(samples):
                    next_label = f"next sample {samples[sample_index + 1]}"
                else:
                    next_label = "end of requested run"

                if os.path.exists(savefilename) and not args.overwrite:
                    _print_runner_block(
                        [
                            f"output already present: {savefilename}",
                            f"skipping {section_label}",
                        ]
                    )
                    skipped_outputs.append((savefilename, sample, subsection))
                    try:
                        output = util.load(savefilename)
                    except Exception as load_error:
                        print(
                            f"warning: could not load skipped output {savefilename}: {load_error}"
                        )
                    continue
                elif os.path.exists(savefilename) and args.overwrite:
                    archived_output = _archive_existing_output(savefilename)
                    _print_runner_block(
                        [
                            f"archived existing output: {archived_output}",
                            f"new output will use: {savefilename}",
                        ]
                    )

                try:
                    if not args.dask:
                        runner = processor.Runner(
                            executor=processor.FuturesExecutor(workers=nworkers),
                            schema=NanoAODSchema,
                            chunksize=chunksize_futures,
                            maxchunks=maxchunks,
                            skipbadfiles=skipbadfiles,
                            xrootdtimeout=500,
                            savemetrics=True,
                        )

                        output, metrics = runner(
                            fileset,
                            treename="Events",
                            processor_instance=TTbarResProcessor(
                                iov=IOV,
                                bkgEst=args.bkgest,
                                noSyst=args.noSyst,
                                deepAK8Cut=args.ttagWP,
                                useDeepAK8=useDeepAK8,
                                useDeepCSV=useDeepCSV,
                                htCut=htCut,
                                anacats=anacats,
                                systematics=systematics,
                                blinding=args.blind,
                                debug=True,
                                produce_ntuple=args.ntuple,
                                sample_metadata=sample_metadata,
                            ),
                        )
                    else:
                        run_instance = processor.Runner(
                            metadata_cache={},
                            executor=processor.DaskExecutor(client=client, retries=2, treereduction=20),
                            schema=NanoAODSchema,
                            savemetrics=True,
                            skipbadfiles=skipbadfiles,
                            chunksize=chunksize_dask,
                            maxchunks=maxchunks,
                        )

                        output, metrics = run_instance(
                            fileset,
                            treename="Events",
                            processor_instance=TTbarResProcessor(
                                iov=IOV,
                                bkgEst=args.bkgest,
                                noSyst=args.noSyst,
                                deepAK8Cut=args.ttagWP,
                                useDeepAK8=useDeepAK8,
                                useDeepCSV=useDeepCSV,
                                htCut=htCut,
                                anacats=anacats,
                                systematics=systematics,
                                blinding=args.blind,
                                produce_ntuple=args.ntuple,
                                sample_metadata=sample_metadata,
                            ),
                        )

                    output["analysisCategories"] = label_map
                    util.save(output, savefilename)
                    print("saving", savefilename)
                    savefilenames.append((savefilename, sample))
                except Exception as exc:
                    failures.append(
                        {
                            "sample": sample,
                            "subsection": subsection,
                            "savefilename": savefilename,
                            "error": repr(exc),
                        }
                    )
                    _print_runner_block(
                        [
                            f"crashed during {section_label}",
                            f"next queued section: {next_label}",
                            "",
                            traceback.format_exc().rstrip(),
                        ]
                    )
                    if args.dask:
                        _print_runner_block(
                            [
                                "restarting Dask client after section failure",
                                f"will retry scheduling from {next_label}",
                            ]
                        )
                        client, cluster = _close_dask_resources(client, cluster)
                        client, cluster = _start_dask_resources(
                            args=args,
                            repo_root=repo_root,
                            upload_to_dask=upload_to_dask,
                            dask_memory=dask_memory,
                            nworkers=nworkers,
                        )
                    continue

    elapsed = time.time() - tic
    printTime(elapsed)
    if metrics is not None:
        print(f"Events/s: {metrics['entries'] / elapsed:.0f}")

    print(
        "run summary:",
        f"saved={len(savefilenames)}",
        f"skipped={len(skipped_outputs)}",
        f"failed={len(failures)}",
    )

    client, cluster = _close_dask_resources(client, cluster)

    return {
        "elapsed": elapsed,
        "metrics": metrics,
        "output": output,
        "savefilenames": savefilenames,
        "skipped_outputs": skipped_outputs,
        "failures": failures,
    }
```

```python
# ---- Build args ------ #
args = build_args()
```

```python
# ---- Run the process ---- #

run_summary = run_analysis(args)
```

```python
import subprocess, os

if args.ntuple:
    for coffea_file, sample in run_summary["savefilenames"]:
        ntuple_dir = os.path.join(os.path.dirname(coffea_file), "ntuples")
        os.makedirs(ntuple_dir, exist_ok=True)
        root_file = os.path.join(
            ntuple_dir, os.path.basename(coffea_file).replace(".coffea", "_ntuple.root")
        )
        print(f"writing ntuple: {coffea_file} -> {root_file}")
        subprocess.run(
            ["python", "write_ntuple.py", coffea_file, root_file, sample], check=True
        )
```

```python
output = run_summary["output"]
print(output["cutflow"])
for key in output:
    print(key)
```

```python
output["normalization"]["applied"]
```

```python
output["ttbarmass"].project("ttbarmass").plot()
```

```python
import matplotlib.pyplot as plt

output["gen_jetmsd_reco_jetmsd"].project("genjetmass").plot(
    label=r"Matched Gen Jet AK8 $m_{SD}$", density=True
)
output["gen_jetmsd_reco_jetmsd"].project("jetmsd").plot(
    label=r"Reco Jet $m_{SD}$", density=True
)
plt.xlabel("Mass (GeV)")
plt.xlim(0, 500)
plt.legend()
```

```python
import mplhep as hep

output["gen_jetmsd_reco_jetmsd"]["nominal", 0, ...].project("genjetmass").plot(
    label=r"Matched Gen Jet AK8 $m_{SD}$", density=True
)
output["gen_jetmsd_reco_jetmsd"]["nominal", 0, ...].project("jetmsd").plot(
    label=r"Reco Jet $m_{SD}$", density=True
)
plt.legend(title=r"$|\eta| < 2.1$")
hep.cms.label(rlabel="Z' 4000 GeV")
```

```python
output["gen_jetmsd_reco_jetmsd"]["nominal", 1, ...].project("genjetmass").plot(
    label=r"Matched Gen Jet AK8 $m_{SD}$", density=True
)
output["gen_jetmsd_reco_jetmsd"]["nominal", 1, ...].project("jetmsd").plot(
    label=r"Reco Jet $m_{SD}$", density=True
)
plt.legend(title=r"$|\eta| > 2.1$")
hep.cms.label(rlabel="Z' 4000 GeV")
```

```python
output["gen_jetmsd_reco_jetmsd"]["nominal", ...].project("jet_nearby").plot()
# plt.legend(title = r"$|\eta| > 2.1$")
hep.cms.label(rlabel="Z' 4000 GeV")
```

```python
output["gen_jetmsd_reco_jetmsd"]["nominal", :, :, "ak4_nearby", ...].project(
    "genjetmass"
).plot(label=r"Matched Gen Jet AK8 $m_{SD}$", density=True)
output["gen_jetmsd_reco_jetmsd"]["nominal", :, :, "ak4_nearby", ...].project(
    "jetmsd"
).plot(label=r"Reco Jet $m_{SD}$", density=True)
plt.legend(title=r"AK4 nearby")
hep.cms.label(rlabel="Z' 4000 GeV")
```

```python
output["gen_jetmsd_reco_jetmsd"]["nominal", :, :, "no_jet_nearby", ...].project(
    "genjetmass"
).plot(label=r"Matched Gen Jet AK8 $m_{SD}$", density=True)
output["gen_jetmsd_reco_jetmsd"]["nominal", :, :, "no_jet_nearby", ...].project(
    "jetmsd"
).plot(label=r"Reco Jet $m_{SD}$", density=True)
plt.legend(title=r"No AK4 nearby")
hep.cms.label(rlabel="Z' 4000 GeV")
```

```python
import coffea
from coffea.dataset_tools import dataset_query
```

```python
from coffea.dataset_tools.rucio_utils import get_dataset_files_replicas
```

```python
outfiles, outsites, sites_counts = get_dataset_files_replicas(
    dataset="/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM",
)
```

```python

```
