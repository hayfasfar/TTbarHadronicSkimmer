
# TTbarHadronicSkimmer

The following setup runs on lpc-el9:

## Login

```bash
ssh -Y -L 8883:127.0.0.1:8883 LPCUSERNAME@cmslpc-el9.fnal.gov
```

## Setup

```bash
cd TTbarHadronicSkimmer
voms-proxy-init --voms cms 
./shell coffeateam/coffea-dask:0.7.22-py3.10-gf48fa  
```

## Running Jobs

```bash
python ttbaranalysis.py --iov 2024 --dataset ZPrime1 
```

> **Note:** You can add `--noSyst` to run without systematics or `--test` to run on 1 chunk.

For now the analysis can run on 2023 and 2024 datasets. JEC and JES has to be added. Please run with `--noSyst`.

## Jupyter Lab 
```bash
jupyter lab --no-browser --ip=127.0.0.1 --port=8883
```

Then copy the provided link to your browser.

## Viewing Histograms

To view basic histograms and systematic variations after running, use [`plots/syst_viewer.ipynb`](plots/syst_viewer.ipynb).