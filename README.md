
# TTbarHadronicSkimmer

This code can be run on either LPC or coffea.casa (hosted at UNL). We highly recommend using coffea.casa because it is upto 5x faster.

## LPC setup

### Login

```bash
ssh -Y -L 8XXX:127.0.0.1:8XXX LPCUSERNAME@cmslpc-el9.fnal.gov
```


### Setup

It is recommended to work on `nobackup` area in LPC.
```bash
cd nobackup
```
Initialise the `voms-proxy`
```bash
voms-proxy-init --rfc --voms cms -valid 192:00
```

Clone coffea-2025 branch of this repository.
```bash
git clone -b coffea-2025 https://github.com/mandalaritra1/TTbarHadronicSkimmer.git
```

```bash
cd TTbarHadronicSkimmer 
```
Setup lpcjobqueue by following instruction from [here](https://github.com/CoffeaTeam/lpcjobqueue). Afterwards, the sigularity container can be run with:

```bash
./shell coffeateam/coffea-dask-almalinux9:2025.12.0-py3.12  
```

### (optional) Jupyter Lab

For interactive jupyter lab environment run the following command inside the singularity container: 
```bash
jupyter lab --no-browser --ip=127.0.0.1 --port=8XXX
```

Then copy the provided link to your browser.

## coffea.casa setup

Go to [coffea.casa](https://coffea.casa) and login using your preferred method.

Select the latest coffea image for 2025 and press `Start`.
![alt text](docs/image.png)

Clone this repository:
```bash
git clone -b coffea-2025 https://github.com/mandalaritra1/TTbarHadronicSkimmer.git
```
Go to the `Dask` tab from left and click the *SHUTDOWN* button.

![alt text](docs/image-1.png)

Open [`ttbaranalysis.ipynb`](ttbaranalysis.ipynb) and run with the CASA configuration.

## Running Jobs

```bash
python ttbaranalysis.py --iov 2024 --dataset ZPrime1 
```

> **Note:** You can add `--noSyst` to run without systematics or `--test` to run on 1 chunk.

> **coffea.casa:** Use the CASA configuration in [`ttbaranalysis.ipynb`](ttbaranalysis.ipynb), or set `--env casa` when running from the command line.

For now the analysis can run on 2023 and 2024 datasets. 

## Viewing Histograms

To view basic histograms and systematic variations after running, use [`plots/syst_viewer.ipynb`](plots/syst_viewer.ipynb).