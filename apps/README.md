# PKL/Coffea Explorer

Launch from the repository root:

```bash
streamlit run apps/coffea_explorer.py
```

The app discovers `.coffea`, `.pkl`, and `.pickle` files under `outputs/` and also accepts manual paths. It can inspect top-level keys, plot `hist.Hist` objects, show cutflows and analysis categories, and plot embedded `output["ntuple"]` branches when present.

Use **File slots** to load up to three files. In the histogram tab, use **Traces** to overlay compatible 1D projections, such as `data` and `ttbar` from the same pickle file, or to compare matching histograms from different files. The **Axis ranges** expander slices numeric axes before projection, so you can plot `ttbarmass` after restricting `jetmass` to a signal window.
