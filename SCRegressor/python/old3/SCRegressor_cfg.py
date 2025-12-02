import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('skipEvents', 
    default=0, 
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
options.parseArguments()

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process("FEVTAnalyzer", Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("DQM.Integration.config.FrontierCondition_GT_Offline_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitConverter_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.TrackRefitter.TTRHBuilder = 'WithAngleAndTemplate'
process.MeasurementTrackerEvent.inactivePixelDetectorLabels = cms.VInputTag() # to prebent worning:  fail to get the list of inactive pixel modules, because of 4.2/4.4 event content change.


process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(1)
    input = cms.untracked.int32(options.maxEvents)
    )

print (" >> Loaded",len(options.inputFiles),"input files from list.")
#evtsToProc = open('list_pi0_ptrecoOgen1p2To1p6_eventsToProcess_.txt').read().splitlines()
#evtsToProc = open('pi0_evtsToProc.txt').read().splitlines()
#evtsToProc = open('eta_evtsToProc.txt').read().splitlines()
#print evtsToProc

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      options.inputFiles
      #'file:miniAOD_HToEleEle.root'
      )
    , skipEvents = cms.untracked.uint32(options.skipEvents)
    #, eventsToProcess = cms.untracked.VEventRange('1:6931:1723687928','1:6932:1723895372')
    #, eventsToProcess = cms.untracked.VEventRange(*evtsToProc)
    #, lumisToProcess = cms.untracked.VLuminosityBlockRange('1:2133-1:2133')
    #, lumisToProcess = cms.untracked.VLuminosityBlockRange('1:3393-1:3393')
    )

process.GlobalTag.globaltag = cms.string('130X_mcRun3_2023_realistic_postBPix_v5') # 2024
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.fevt = cms.EDAnalyzer('SCRegressor'
    #, EBRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEB')
    , gsfElectronCollection = cms.InputTag('gedGsfElectrons')
    , photonCollection = cms.InputTag('slimmedPhotons')
    , jetCollection = cms.InputTag("ak4PFJets")
    , muonCollection = cms.InputTag('slimmedMuons')
    #, electronCollection = cms.InputTag('slimmedElectrons')
    , electronCollection = cms.InputTag("gedGsfElectrons")
    , EBRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEB')
    , EERecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEE')
    , ESRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsES')
    , reducedAODEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
    , reducedAODEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
    , reducedAODESRecHitCollection = cms.InputTag('reducedEcalRecHitsES')
    , reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
    , reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
    , reducedESRecHitCollection = cms.InputTag('reducedEcalRecHitsES')
    #, reducedEBRecHitCollection = cms.InputTag('reducedEgamma:reducedEBRecHits')
    #, reducedEERecHitCollection = cms.InputTag('reducedEgamma:reducedEERecHits')
    #, reducedESRecHitCollection = cms.InputTag('reducedEgamma:reducedESRecHits')
    #, reducedHBHERecHitCollection    = cms.InputTag('reducedEgamma:reducedHBHEHits')
    , reducedHBHERecHitCollection    = cms.InputTag('hbhereco')
    #, genParticleCollection = cms.InputTag('prunedGenParticles')
    , genParticleCollection = cms.InputTag('genParticles')
    , genJetCollection = cms.InputTag('ak4GenJets')
    , trackCollection = cms.InputTag("generalTracks")
    , vertexCollection               = cms.InputTag("offlinePrimaryVertices")
    , transTrackBuilder              = cms.ESInputTag("", "TransientTrackBuilder")
    , rhoLabel = cms.InputTag("fixedGridRhoFastjetAll")
    , trgResults = cms.InputTag("TriggerResults","","HLT")
    , generator = cms.InputTag("generator")
    , lhe = cms.InputTag("lhe")
    )

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string('histo.root')
    fileName = cms.string(options.outputFile)
    )

process.hltFilter = cms.EDFilter("HLTHighLevel",
                                          eventSetupPathsKey = cms.string(''),
                                          TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                          HLTPaths = cms.vstring('*'),
                                          andOr = cms.bool(True),
                                          throw = cms.bool(False)
                                          )

process.p = cms.Path(
  process.siStripMatchedRecHits*process.siPixelRecHits*process.MeasurementTrackerEvent*process.TrackRefitter*
  process.hltFilter*
  #process.fullPatMetSequenceModifiedMET*
  #process.egammaPostRecoSeq*
  process.fevt
)
