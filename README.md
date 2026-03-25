# MLAnalyzer_run3

## Setup

To use `MLAnalyzer_run3`, you must start a CMSSW installation and add several repos before this one. Follow the lines below to do so.

```
# Set environment on cmslpc
export SCRAM_ARCH=el8_amd64_gcc11
source /cvmfs/cms.cern.ch/cmsset_default.sh

# Install CMSSW and activate cmsenv
cmsrel CMSSW_13_0_13
cd CMSSW_13_0_13/src
eval `scram runtime -sh`
cmsenv

# Clone the necessary repos
git clone https://github.com/svfit/ClassicSVfit TauAnalysis/ClassicSVfit -b fastMTT_21_06_2018
git clone https://github.com/svfit/SVfitTF TauAnalysis/SVfitTF
git clone https://github.com/crovellac/MLAnalyzer_run3.git
scram b -j 16

# Run the analyzer
cd MLAnalyzer_run3
python3 runRHAnalyzer.py # this equivalent to cmsRun command (pythonic version only)
```

## Plotting

There are several scripts for plotting the mass vs pT distribution, which is used in unbiasing.
- `mVSpt_from_ntuple.py`: Takes a list of ntuples output by `SCRegressor`, and makes `hist.root` which has the mass vs pT plot.
- `mVSpt_from_h5.py`: Takes a list of h5 files created by Run_3_convert_root2h5_SC.py, and makes `hist.root` which has the mass vs pT plot.
- `plot_mVSpt_from_root.py`: Opens `hist.root` and creates the graphical plot, and also prints out the bin counts to be used in [unbiasing](https://github.com/crovellac/MCProduction/blob/run3/GeneratorInterface/Pythia8Interface/plugins/Py8PtV3Gun.cc).

## How to use

Samples can be generated using [MCProduction/E2E-AToEleEle](https://github.com/crovellac/MCProduction/tree/run3/E2E-AToEleEle). Following all three steps results in AOD samples with pileup included.

These AOD samples can be used as input for `SCRegressor` to create ntuples.

If the mass vs pT distribution from the ntuples looks biased after plotting it with `mVSpt_from_ntuple.py` and `plot_mVSpt_from_root.py`, one can use the bin counts output by `plot_mVSpt_from_root.py` in the unbiasing part of [Py8PtV3Gun.cc](https://github.com/crovellac/MCProduction/blob/run3/GeneratorInterface/Pythia8Interface/plugins/Py8PtV3Gun.cc#L117-L155).

When the distribution in the ntuples looks unbiased, convert the ntuples to h5 files using `Run_3_convert_root2h5_SC.py`.

## Other repos

- [MCProduction](https://github.com/crovellac/MCProduction/tree/run3): Code for generating Monte Carlo samples, which are used as input to `SCRegressor`.
- [Gen](https://github.com/crovellac/Gen): Analyzer for looking at gen-level samples from `MCProduction.`
