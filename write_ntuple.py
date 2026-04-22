#!/usr/bin/env python
"""Convert the ntuple section of one or more coffea output files to a ROOT TTree.

Usage (single file):
    python write_ntuple.py <input.coffea> <output.root> [tree_name]

Usage (merge multiple eras into one file):
    python write_ntuple.py <input1.coffea> <input2.coffea> ... --out <output.root> [--tree <tree_name>]

No ROOT installation required — uproot writes the file in pure Python.
"""
import sys
import argparse
import numpy as np
import uproot
from coffea import util


def load_branches(coffea_file):
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


if __name__ == "__main__":
    if "--out" not in sys.argv:
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
