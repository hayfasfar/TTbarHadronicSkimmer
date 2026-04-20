#!/usr/bin/env python
"""Convert the ntuple section of a coffea output file to a ROOT TTree.

Usage:
    python write_ntuple.py <input.coffea> <output.root> [tree_name]

No ROOT installation required — uproot writes the file in pure Python.
"""
import sys
import uproot
from coffea import util


def write_ntuple(coffea_file, root_file, tree_name="ttbar"):
    output = util.load(coffea_file)
    ntuple = output["ntuple"]

    branches = {col: ntuple[col].value for col in ntuple.keys()}

    n_events = len(next(iter(branches.values())))
    print(f"Writing {n_events} events, {len(branches)} branches → {root_file}:{tree_name}")

    with uproot.recreate(root_file) as f:
        f.mktree(tree_name, {col: arr.dtype for col, arr in branches.items()})
        f[tree_name].extend(branches)

    print("Done.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    tree_name = sys.argv[3] if len(sys.argv) > 3 else "ttbar"
    write_ntuple(sys.argv[1], sys.argv[2], tree_name)
