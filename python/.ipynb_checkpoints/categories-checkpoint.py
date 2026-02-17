import itertools
import numpy as np


def build_analysis_categories(antitag, ttag_s0, ttag_s1, rapidity, anacats):
    """Build category masks used by the processor from top-tag and rapidity regions."""

    ttag2 = (ttag_s0 & ttag_s1)
    cen = (np.abs(rapidity) < 1.0)
    fwd = (~cen)

    regs = {"cen": cen, "fwd": fwd}
    ttags = {
        "at": antitag,  # 2Dalphabet fail region
        "2t": ttag2,    # 2Dalphabet pass region
    }

    categories = {
        t[0] + y[0]: (t[1] & y[1])
        for t, y in itertools.product(ttags.items(), regs.items())
    }
    return {label: categories[label] for label in anacats}
