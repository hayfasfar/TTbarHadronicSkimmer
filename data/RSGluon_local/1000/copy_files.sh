#!/bin/bash

# Create an array of file paths
files=(

"/store/mc/RunIISummer20UL17NanoAODv9/RSGluonToTT_M-1000_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/130000/4531062A-4205-FF48-8DFD-4C8F69A38332.root",
"/store/mc/RunIISummer20UL17NanoAODv9/RSGluonToTT_M-1000_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/130000/BF5F806D-C435-B948-8EC8-927061DC3C6C.root",
"/store/mc/RunIISummer20UL17NanoAODv9/RSGluonToTT_M-1000_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/130000/B4C2DEA9-F044-4545-B906-39A52E700575.root",
"/store/mc/RunIISummer20UL17NanoAODv9/RSGluonToTT_M-1000_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/70000/8993AEAD-C88B-C14E-B2CE-E6ADA1B272DE.root",
"/store/mc/RunIISummer20UL17NanoAODv9/RSGluonToTT_M-1000_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/70000/2B8A9044-E287-4E42-B5DA-5049E4E692AB.root",
"/store/mc/RunIISummer20UL17NanoAODv9/RSGluonToTT_M-1000_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/70000/6F7D46DA-05B3-1044-8F28-21B22394DA29.root",
"/store/mc/RunIISummer20UL17NanoAODv9/RSGluonToTT_M-1000_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/70000/A7A0DFD3-0D6B-DC48-B7D7-AC9EA237DAAF.root",
"/store/mc/RunIISummer20UL17NanoAODv9/RSGluonToTT_M-1000_TuneCP5_13TeV-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/280000/B62B9503-7B87-C94E-904A-026F799B8453.root"
)

# Loop over the files and copy each one
for file in "${files[@]}"; do
    xrdcp "root://cmsxrootd.fnal.gov//$file" .
done

