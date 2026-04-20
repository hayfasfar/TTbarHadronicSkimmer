import hist
import numpy as np
from coffea import processor


manual_bins = [400, 500, 600, 800, 1000, 1500, 2000, 3000, 7000, 10000]


def _flatten_hist_tree(tree):
    flat = {}
    for value in tree.values():
        for key, histo in value.items():
            flat[key] = histo
    return flat


def build_output_histograms(anacats, systematics, no_syst, produce_ntuple=False):
    syst_category_strings = ["nominal"]
    if not no_syst:
        for s in systematics:
            if s == "nominal":
                continue
            if "hem" in s:
                syst_category_strings.append(s)
            else:
                syst_category_strings.extend([s + "Down", s + "Up"])

    syst_axis = hist.axis.StrCategory(syst_category_strings, name="systematic")
    cats_axis = hist.axis.IntCategory(range(len(anacats)), name="anacat", label="Analysis Category")
    ttbarmass2D_axis = hist.axis.Regular(92, 800, 10000, name="ttbarmass", label=r"$m_{t\bar{t}}$ [GeV]")
    jetmass2D_axis = hist.axis.Regular(100, 0, 500, name="jetmass", label=r"Jet $m_{SD}$ [GeV]")
    jetmsd_axis = hist.axis.Regular(100, 0, 500, name="jetmsd", label=r"Jet $m_{SD}$ [GeV]")
    ht_axis = hist.axis.Regular(40, 400, 4400, name="ht", label=r"$H_T$ [GeV]")
    manual_axis = hist.axis.Variable(manual_bins, name="jetp", label=r"Jet Momentum [GeV]")
    jetdy_axis = hist.axis.Regular(50, -3, 3, name="jetdy", label=r"$\Delta y$")
    chi_axis   = hist.axis.Regular(50, 1, 30, name="chi",   label=r"$\chi_\mathrm{dijet}$")
    jetpt_axis = hist.axis.Regular(50, 300, 2000, name="jetpt", label=r"Jet $p_T$ [GeV]")
    jeteta_axis = hist.axis.Regular(50, -3, 3, name="jeteta", label=r"Jet $\eta$")
    jetphi_axis = hist.axis.Regular(50, -3.14159, 3.14159, name="jetphi", label=r"Jet $\phi$")
    jety_axis = hist.axis.Regular(50, -3, 3, name="jety", label=r"Jet rapidity")
    jetdr_axis = hist.axis.Regular(50, 0, 5, name="dr", label=r"$\Delta R$")
    gentopmass_axis = hist.axis.Regular(100, 0, 500, name="gentopmass", label=r"Gen top mass [GeV]")
    genjetmass_axis = hist.axis.Regular(100, 0, 500, name="genjetmass", label=r"Gen jet mass [GeV]")
    massres_axis = hist.axis.Regular(80, -1.5, 1.5, name="massres", label=r"$(m_{reco}-m_{gen})/m_{gen}$")
    abs_eta_axis = hist.axis.Variable([0.0, 2.1, 5.0], name="abs_eta", label=r"$|\eta|$")
    jet_nearby_axis = hist.axis.StrCategory(["no_jet_nearby", "ak4_nearby", "ak8_nearby"], name="jet_nearby")

    hist_tree = {
        "mass": {
            "ttbarmass": hist.Hist(syst_axis, cats_axis, ttbarmass2D_axis, storage="weight", name="Counts"),
            "mtt_unwgt": hist.Hist(syst_axis, cats_axis, ttbarmass2D_axis, storage="weight", name="Counts"),
            "mtt_vs_mt": hist.Hist(
                syst_axis, cats_axis, jetmass2D_axis, ttbarmass2D_axis, storage="weight", name="Counts"
            ),
        },
        "mistag": {
            "numerator": hist.Hist(cats_axis, manual_axis, storage="weight", name="Counts"),
            "denominator": hist.Hist(cats_axis, manual_axis, storage="weight", name="Counts"),
        },
        "jets": {
            "jetmass": hist.Hist(syst_axis, cats_axis, jetmass2D_axis, storage="weight", name="Counts"), # Selected jet mass
            "jetmsd": hist.Hist(syst_axis, cats_axis, jetmsd_axis, storage="weight", name="Counts"), # Selected jet softdrop mass
            "jetdy": hist.Hist(syst_axis, cats_axis, jetdy_axis, storage="weight", name="Counts"),
            "chi":   hist.Hist(syst_axis, cats_axis, chi_axis,   storage="weight", name="Counts"),
            "jetmass1": hist.Hist(syst_axis, cats_axis, jetmass2D_axis, storage="weight", name="Counts"),
            "jetmsd1": hist.Hist(syst_axis, cats_axis, jetmsd_axis, storage="weight", name="Counts"),
            "dR_min_jet2": hist.Hist(syst_axis, cats_axis, jetdr_axis, ttbarmass2D_axis, storage="weight", name="Counts"),
            "jet0_pt": hist.Hist(syst_axis, cats_axis, jetpt_axis, storage="weight", name="Counts"),
            "jet0_eta": hist.Hist(syst_axis, cats_axis, jeteta_axis, storage="weight", name="Counts"),
            "jet0_phi": hist.Hist(syst_axis, cats_axis, jetphi_axis, storage="weight", name="Counts"),
            "jet0_rapidity": hist.Hist(syst_axis, cats_axis, jety_axis, storage="weight", name="Counts"),
            "jet1_pt": hist.Hist(syst_axis, cats_axis, jetpt_axis, storage="weight", name="Counts"),
            "jet1_eta": hist.Hist(syst_axis, cats_axis, jeteta_axis, storage="weight", name="Counts"),
            "jet1_phi": hist.Hist(syst_axis, cats_axis, jetphi_axis, storage="weight", name="Counts"),
            "jet1_rapidity": hist.Hist(syst_axis, cats_axis, jety_axis, storage="weight", name="Counts"),
        },
        "truth": {
            "gen_mt": hist.Hist(syst_axis, cats_axis, gentopmass_axis, storage="weight", name="Counts"),
            "gen_mttbar": hist.Hist(syst_axis, cats_axis, ttbarmass2D_axis, storage="weight", name="Counts"),
            "jet0_gen_dr": hist.Hist(syst_axis, cats_axis, jetdr_axis, storage="weight", name="Counts"),
            "jet1_gen_dr": hist.Hist(syst_axis, cats_axis, jetdr_axis, storage="weight", name="Counts"),
            "jet_mass_resolution": hist.Hist(
                syst_axis, cats_axis, abs_eta_axis, jet_nearby_axis, massres_axis, storage="weight", name="Counts"
            ),
            "gen_jetmsd_reco_jetmsd": hist.Hist(
                syst_axis,
                cats_axis,
                abs_eta_axis,
                jet_nearby_axis,
                genjetmass_axis,
                jetmsd_axis,
                storage="weight",
                name="Counts",
            ),
        },
        "event": {
            "ht": hist.Hist(syst_axis, cats_axis, ht_axis, storage="weight", name="Counts"),
        },
    }

    output = _flatten_hist_tree(hist_tree)
    output.update(
        {
            "cutflow": processor.defaultdict_accumulator(int),
            "weights": processor.defaultdict_accumulator(float),
            "systematics": processor.defaultdict_accumulator(float),
            "event_list": processor.dict_accumulator(
                {
                    "run": processor.list_accumulator([]),
                    "lumi": processor.list_accumulator([]),
                    "event": processor.list_accumulator([]),
                }
            ),
            "truthstudy": processor.defaultdict_accumulator(int),
        }
    )
    if produce_ntuple:
        output["ntuple"] = build_ntuple_accumulators()
    return output


_NTUPLE_FLOAT_COLS = [
    "jet0_pt", "jet0_eta", "jet0_phi", "jet0_msd", "jet0_tdisc",
    "jet1_pt", "jet1_eta", "jet1_phi", "jet1_msd", "jet1_tdisc",
    "ttbarmass", "ht", "dy", "chi", "weight",
    "jet0_rapidity", "jet1_rapidity",
]
_NTUPLE_INT_COLS = ["anacat", "run", "lumi", "event"]


def build_ntuple_accumulators():
    acc = {col: processor.column_accumulator(np.array([], dtype=np.float32))
           for col in _NTUPLE_FLOAT_COLS}
    acc.update({col: processor.column_accumulator(np.array([], dtype=np.int64))
                for col in _NTUPLE_INT_COLS})
    return processor.dict_accumulator(acc)
