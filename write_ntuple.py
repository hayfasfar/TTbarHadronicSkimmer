#!/usr/bin/env python
"""Convert the ntuple section of one or more coffea output files to a ROOT TTree.

Usage (single file):
    python write_ntuple.py <input.coffea> <output.root> [tree_name]

Usage (merge multiple eras into one file):
    python write_ntuple.py <input1.coffea> <input2.coffea> ... --out <output.root> [--tree <tree_name>]

Usage (merge chunk ROOT files, optionally scaling the weight branch):
    python write_ntuple.py --merge-root <chunk1.root> <chunk2.root> ... --out <output.root> [--weight-scale <scale>]

No ROOT installation required — uproot writes the file in pure Python.
"""
import sys
import argparse
import numpy as np
import uproot


def load_branches(coffea_file):
    from coffea import util

    output = util.load(coffea_file)
    ntuple = output["ntuple"]
    return {col: ntuple[col].value for col in ntuple.keys()}


def merge_branches(branch_list):
    keys = branch_list[0].keys()
    return {col: np.concatenate([b[col] for b in branch_list]) for col in keys}


def write_ntuple(coffea_files, root_file, tree_name="ttbar"):
    all_branches = []
    for f in coffea_files:
        print(f"Loading {f} ...")
        all_branches.append(load_branches(f))

    branches = merge_branches(all_branches) if len(all_branches) > 1 else all_branches[0]

    n_events = len(next(iter(branches.values())))
    print(f"Writing {n_events} events, {len(branches)} branches → {root_file}:{tree_name}")

    with uproot.recreate(root_file) as f:
        f.mktree(tree_name, {col: arr.dtype for col, arr in branches.items()})
        f[tree_name].extend(branches)

    print("Done.")


def merge_root_ntuples(root_files, output_file, tree_name="ttbar", weight_scale=1.0):
    root_files = list(root_files)
    if not root_files:
        print(f"No ntuple chunk files found for {output_file}; skipping merge.", flush=True)
        return 0

    n_events = 0
    tree_created = False
    with uproot.recreate(output_file) as fout:
        for root_file in root_files:
            print(f"Merging {root_file} ...", flush=True)
            with uproot.open(root_file) as fin:
                tree = fin[tree_name]
                arrays = tree.arrays(library="np")

            if weight_scale != 1.0 and "weight" in arrays:
                arrays["weight"] = (arrays["weight"] * weight_scale).astype(arrays["weight"].dtype)

            if not tree_created:
                fout.mktree(tree_name, {col: arr.dtype for col, arr in arrays.items()})
                tree_created = True

            fout[tree_name].extend(arrays)
            n_events += len(next(iter(arrays.values()))) if arrays else 0

    print(f"Wrote {n_events} events from {len(root_files)} chunks → {output_file}:{tree_name}", flush=True)
    return n_events


if __name__ == "__main__":
    if "--merge-root" in sys.argv:
        parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("--merge-root", nargs="+", metavar="chunk.root", help="Chunk ROOT files to merge")
        parser.add_argument("--out", required=True, metavar="output.root", help="Output ROOT file")
        parser.add_argument("--tree", default="ttbar", metavar="tree_name", help="TTree name (default: ttbar)")
        parser.add_argument("--weight-scale", type=float, default=1.0, help="Scale factor applied to the weight branch")
        args = parser.parse_args()
        merge_root_ntuples(args.merge_root, args.out, args.tree, args.weight_scale)
    elif "--out" not in sys.argv:
        # Legacy single-file usage: write_ntuple.py input.coffea output.root [tree_name]
        if len(sys.argv) < 3:
            print(__doc__)
            sys.exit(1)
        tree_name = sys.argv[3] if len(sys.argv) > 3 else "ttbar"
        write_ntuple([sys.argv[1]], sys.argv[2], tree_name)
    else:
        parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("inputs", nargs="+", metavar="input.coffea", help="One or more coffea files to merge")
        parser.add_argument("--out", required=True, metavar="output.root", help="Output ROOT file")
        parser.add_argument("--tree", default="ttbar", metavar="tree_name", help="TTree name (default: ttbar)")
        args = parser.parse_args()
        write_ntuple(args.inputs, args.out, args.tree)
