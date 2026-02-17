
The following setup runs on lpc-el9:


to login:
 
```
ssh -Y -L 8883:127.0.0.1:8883 LPCUSERNAME@cmslpc-el9.fnal.gov
```

To start the setup : 
```
cd TTbarHadronicSkimmer
voms-proxy-init --voms cms 
./shell coffeateam/coffea-dask:0.7.22-py3.10-gf48fa  

```
To run jobs: 

```
python ttbaranalysis.py --iov 2024 --dataset ZPrime1 

```

You can cadd --noSyst to run without systematic or --test to run on 1 chunk

for now the analysis can run on 2023 and 2024 dataset. JEC and JES has to be added. Please run with --noSyst.

If you want to modify text on shell independent way do: 

```
jupyter lab --no-browser --ip=127.0.0.1 --port=8883
```

Then copy the provided link to your browser.

