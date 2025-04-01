import os

cfg='SCRegressor/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/home-c/ccrovell/run3_massregression/CMSSW_13_0_13/src/MCProduction/E2E-HToEleEle/miniAOD_HToEleEle.root'
inputFiles_='file:/eos/home-c/ccrovell/run3_massregression/CMSSW_13_0_13/src/MCProduction/E2E-HToEleEle/SIM_HToEleEle_m0p1To6_pT20To150_ctau0To3_eta0To1p4_pythia8.root'#pixel checks

maxEvents_  = -1
skipEvents_ =  0 
outputFile_ =  'HToEleEle_ntuple.root'

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
