

import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams

fevt = cms.EDAnalyzer('SCRegressor'
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
    , secVertexCollection            = cms.InputTag("inclusiveCandidateSecondaryVertices")
    , siPixelRecHitCollection        = cms.InputTag("siPixelRecHits")
    , siStripMatchedRecHitCollection = cms.InputTag("siStripMatchedRecHits", "matchedRecHit")
    , siStripRphiRecHits             = cms.InputTag("siStripMatchedRecHits", "rphiRecHit")
    , siStripUnmatchedRphiRecHits    = cms.InputTag("siStripMatchedRecHits", "rphiRecHitUnmatched")
    , siStripStereoRecHits           = cms.InputTag("siStripMatchedRecHits", "stereoRecHit")
    , siStripUnmatchedStereoRecHits  = cms.InputTag("siStripMatchedRecHits", "stereoRecHitUnmatched")
    , pfCollection                   = cms.InputTag("particleFlow")
    , srcPFCandidates                = cms.InputTag("particleFlow")
    , recoJetsForBTagging            = cms.InputTag("ak4PFJetsCHS")
    , jetTagCollection               = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags")

    , transTrackBuilder              = cms.ESInputTag("", "TransientTrackBuilder")
    , srcLeptons                     = cms.VInputTag("gedGsfElectrons","muons","gedPhotons")
    , rhoLabel = cms.InputTag("fixedGridRhoFastjetAll")
    
    # Jet level cfg
    , nJets     = cms.int32(-1)
    , minJetPt  = cms.double(20.)
    , maxJetEta = cms.double(2.4)
    , z0PVCut   = cms.double(0.1)

    , trgResults = cms.InputTag("TriggerResults","","HLT")
    , generator = cms.InputTag("generator")
    , lhe = cms.InputTag("lhe")
    )
