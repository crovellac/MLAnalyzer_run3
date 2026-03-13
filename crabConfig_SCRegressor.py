from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters

CFG = 'AToEleEle_m0p01To1p2_pT25To160_ctau0To3_eta0To1p4_SCRegressor'

# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName># To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>

# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea = 'crab_SCRegressor'
config.General.requestName = CFG
config.General.transferOutputs = True
config.General.transferLogs = False

# CMS cfg file goes here:
config.JobType.pluginName = 'Analysis'
config.JobType.numCores=4
config.JobType.psetName = 'SCRegressor/python/ConfFile_cfg.py' # analyzer cfg file
config.JobType.maxMemoryMB = 4000

# Define input and units per job here:
#config.Data.userInputFiles = open('MLAnalyzer/list_production.txt'%idx).readlines()
#config.Data.userInputFiles = open('MLAnalyzer/list_prod_unbiased.txt').readlines()
config.Data.userInputFiles = open('list_aod_files.txt').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 # units: as defined by config.Data.splitting
config.Data.totalUnits = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission #config.Data.totalUnits = 10 # test production
config.Data.publication = False

# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/user/ccrovell/' # add your username as subdirectory
config.Site.whitelist = ['T3_US_FNALLPC','T2_CH_CERN']
config.Data.outputPrimaryDataset = 'AToEleEle_m0p01To1p2_pT40To160_ctau0p0_eta0To1p4_SCRegressor_biased'
config.Data.outputDatasetTag = config.General.requestName
