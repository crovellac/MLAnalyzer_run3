import os

cfg='SCRegressor/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/home-c/ccrovell/run3_massregression/CMSSW_13_0_13/src/MCProduction/E2E-HToEleEle/miniAOD_HToEleEle.root'
inputFiles_='file:AOD_AToEleEle_1.root'#pixel checks

maxEvents_  = -1
skipEvents_ =  0 
outputFile_ =  'HToEleEle_ntuple.root'

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
