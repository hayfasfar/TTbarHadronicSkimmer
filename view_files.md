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
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# add any other imports you need here
# Set the file to inspect — swap in any path from qcd.json
filename = "/store/mc/Run3Winter24NanoAOD/QCD_PT-15to7000_TuneCP5_Flat2022_13p6TeV_pythia8/NANOAODSIM/JMENanoV14_133X_mcRun3_2024_realistic_v10-v1/2530000/01eebb7d-8f73-4e7a-baae-83c9ff84af44.root"
filename = "/store/mc/RunIII2024Summer24NanoAODv15/QCD_Bin-PT-15to7000_Par-PT-Flat_TuneCH3_13p6TeV_herwig7/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/2820000/005f8b3c-b044-48c0-bf57-77db0939ade4.root"
redirector = "root://xcache/"

events = NanoEventsFactory.from_root(
    {redirector + filename: "Events"},
    schemaclass=NanoAODSchema,
    entry_stop=1,  # only need 1 event to read schema
).events()
```

```python
print("=== Top-level event fields ===")
print(sorted(events.fields))
```


```python
print("=== FatJet fields ===")
for f in sorted(events.FatJet.fields):
    print(f"  FatJet_{f}")
```


```python
print("=== FatJet tagger-related fields ===")
keywords = ["top", "ttag", "particleNet", "globalParT", "deepTag", "ParT", "Xtt", "Top"]
for f in sorted(events.FatJet.fields):
    if any(k.lower() in f.lower() for k in keywords):
        print(f"  FatJet_{f}")
```


```python

```
