#ifndef SCRegressor_h
#define SCRegressor_h
// -*- C++ -*-
//
// Package:    MLAnalyzer/SCRegressor
// Class:      SCRegressor
//
//
// Original Author:  Michael Andrews
//         Created:  Mon, 17 Jul 2017 15:59:54 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection defined here
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMath.h"
#include "TProfile2D.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DQM/SiPixelMonitorRecHit/interface/SiPixelRecHitModule.h"
#include "DQM/HcalCommon/interface/Constants.h"
//
// class declaration
//
using pat::MuonCollection;
using pat::MuonRef;
using pat::ElectronCollection;
using pat::ElectronRef;
using pat::JetCollection;
using pat::JetRef;
using pat::PhotonCollection;
using pat::PhotonRef;
//using reco::PhotonCollection;
//using reco::PhotonRef;

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SCRegressor : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit SCRegressor(const edm::ParameterSet&);
    ~SCRegressor();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    //edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronCollectionT_;
    //edm::EDGetTokenT<reco::GsfElectronCollection> electronCollectionT_;
    edm::EDGetTokenT<ElectronCollection> electronCollectionT_;
    edm::EDGetTokenT<MuonCollection> muonCollectionT_;
    edm::EDGetTokenT<PhotonCollection> photonCollectionT_;
    //edm::EDGetTokenT<std::vector<reco::Photon>> photonCollectionT_;
    edm::EDGetTokenT<JetCollection> jetCollectionT_;
    //edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    //edm::EDGetTokenT<EcalRecHitCollection> EERecHitCollectionT_;
    //edm::EDGetTokenT<EcalRecHitCollection> ESRecHitCollectionT_;
    edm::EDGetTokenT<edm::SortedCollection<EcalRecHit>> EBRecHitCollectionT_;
    edm::EDGetTokenT<edm::SortedCollection<EcalRecHit>> EERecHitCollectionT_;
    edm::EDGetTokenT<edm::SortedCollection<EcalRecHit>> ESRecHitCollectionT_;
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit>> HBHERecHitCollectionT_;

    edm::EDGetTokenT<EcalRecHitCollection> AODEBRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> AODEERecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> AODESRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> RECOEBRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> RECOEERecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> RECOESRecHitCollectionT_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetCollectionT_;
    //edm::EDGetTokenT<reco::TrackCollection> trackCollectionT_;
    edm::EDGetTokenT<pat::IsolatedTrackCollection> trackCollectionT_;
    edm::EDGetTokenT<double> rhoLabel_;
    edm::EDGetTokenT<edm::TriggerResults> trgResultsT_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoT_;
    edm::EDGetTokenT<LHEEventProduct> lheEventT_;
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

    static const int nEE = 2;
    static const int nElectrons = 4;

    TProfile2D *hEB_energy;
    TProfile2D *hEB_time;
    std::vector<float> vEB_energy_;
    std::vector<float> vEB_time_;
    TProfile2D *hEE_energy[nEE];
    TProfile2D *hEE_time[nEE];
    std::vector<float> vEE_energy_[nEE];
    std::vector<float> vEE_time_[nEE];
    TH2F *hEvt_EE_energy[nEE];
    TProfile2D *hECAL_energy;
    std::vector<float> vECAL_energy_;
    TH2F *hEvt_HBHE_energy;
    TProfile2D *hHBHE_energy_EB;
    TProfile2D *hHBHE_energy;
    std::vector<float> vHBHE_energy_EB_;
    std::vector<float> vHBHE_energy_;

    //TH1D * histo;
    TProfile2D * hSC_energy;
    TProfile2D * hSC_time;
    TH1F * hSC_mass;
    TH1F * hDR;
    TH1F * hdEta;
    TH1F * hdPhi;
    TH3F * hdPhidEtaM;
    TH2F * hdPhidEta;
    TH1F * hSC_pT;

    TTree* RHTree;
    
    unsigned int nPho;
    unsigned int nEle;
    unsigned long long eventId_;
    unsigned int runId_;
    unsigned int lumiId_;

    void branchesSC ( TTree*, edm::Service<TFileService>& );
    void branchesSCaod ( TTree*, edm::Service<TFileService>& );
    void branchesSCreco ( TTree*, edm::Service<TFileService>& );
    void branchesEB ( TTree*, edm::Service<TFileService>& );
    void branchesEE ( TTree*, edm::Service<TFileService>& );
    void branchesECALstitched ( TTree*, edm::Service<TFileService>& );
    void branchesTracksAtEBEE ( TTree*, edm::Service<TFileService>& );
    void branchesHBHE ( TTree*, edm::Service<TFileService>& );
    //void branchesPhoVars ( TTree*, edm::Service<TFileService>& );
    void branchesEvtWgt ( TTree*, edm::Service<TFileService>& );

    void fillSC     ( const edm::Event&, const edm::EventSetup& );
    void fillSCaod  ( const edm::Event&, const edm::EventSetup& );
    void fillSCreco ( const edm::Event&, const edm::EventSetup& );
    void fillEB     ( const edm::Event&, const edm::EventSetup& );
    void fillEE     ( const edm::Event&, const edm::EventSetup& );
    void fillECALstitched     ( const edm::Event&, const edm::EventSetup& );
    void fillTracksAtEBEE ( const edm::Event&, const edm::EventSetup& );
    void fillHBHE   ( const edm::Event&, const edm::EventSetup& );
    //void fillPhoVars ( const edm::Event&, const edm::EventSetup& );
    void fillEvtWgt ( const edm::Event&, const edm::EventSetup& );

    void branchesPiSel       ( TTree*, edm::Service<TFileService>& );
    void branchesPhotonSel   ( TTree*, edm::Service<TFileService>& );
    void branchesDiPhotonSel ( TTree*, edm::Service<TFileService>& );
    void branchesZJetsEleSel ( TTree*, edm::Service<TFileService>& );
    void branchesZJetsMuSel ( TTree*, edm::Service<TFileService>& );
    void branchesNJetsSel ( TTree*, edm::Service<TFileService>& );
    void branchesH2aaSel     ( TTree*, edm::Service<TFileService>& );
    void branchesQCDSel     ( TTree*, edm::Service<TFileService>& );
    bool runPiSel        ( const edm::Event&, const edm::EventSetup& );
    bool runPhotonSel    ( const edm::Event&, const edm::EventSetup& );
    bool runDiPhotonSel  ( const edm::Event&, const edm::EventSetup& );
    bool runZJetsEleSel  ( const edm::Event&, const edm::EventSetup& );
    bool runZJetsMuSel  ( const edm::Event&, const edm::EventSetup& );
    bool runNJetsSel  ( const edm::Event&, const edm::EventSetup& );
    bool runH2aaSel      ( const edm::Event&, const edm::EventSetup& );
    bool runQCDSel      ( const edm::Event&, const edm::EventSetup& );
    void fillPiSel       ( const edm::Event&, const edm::EventSetup& );
    void fillPhotonSel   ( const edm::Event&, const edm::EventSetup& );
    void fillDiPhotonSel ( const edm::Event&, const edm::EventSetup& );
    void fillZJetsEleSel ( const edm::Event&, const edm::EventSetup& );
    void fillZJetsMuSel ( const edm::Event&, const edm::EventSetup& );
    void fillNJetsSel ( const edm::Event&, const edm::EventSetup& );
    void fillH2aaSel     ( const edm::Event&, const edm::EventSetup& );
    void fillQCDSel     ( const edm::Event&, const edm::EventSetup& );
    void fillECAL_with_EEproj ( TH2F*, int, int );

    std::map<unsigned int, std::vector<unsigned int>> mGenPi0_RecoPho;
    std::vector<int> vPreselPhoIdxs_;
    std::vector<int> vRegressPhoIdxs_;
    std::vector<int> vGenAIdxs_;
    std::vector<float> vIphi_Emax_;
    std::vector<float> vIeta_Emax_;

    //std::vector<std::vector<float>> vEB_SCenergy_;
    std::vector<std::vector<float>> vSC_energy_;
    std::vector<std::vector<float>> vSC_energyT_;
    std::vector<std::vector<float>> vSC_energyZ_;
    std::vector<std::vector<float>> vSC_time_;

    TProfile2D * hSCaod_energy;
    TProfile2D * hSCaod_time;
    std::vector<std::vector<float>> vSCaod_energy_;
    std::vector<std::vector<float>> vSCaod_energyT_;
    std::vector<std::vector<float>> vSCaod_energyZ_;
    std::vector<std::vector<float>> vSCaod_time_;
    TProfile2D * hSCreco_energy;
    TProfile2D * hSCreco_time;
    std::vector<std::vector<float>> vSCreco_energy_;
    std::vector<std::vector<float>> vSCreco_energyT_;
    std::vector<std::vector<float>> vSCreco_energyZ_;
    std::vector<std::vector<float>> vSCreco_time_;

    //TH2F *hTracks_EE[nEE];
    TH2F *hTracks_EB;
    //TH2F *hTracksPt_EE[nEE];
    TH2F *hTracksPt_EB;
    //std::vector<float> vTracksPt_EE_[nEE];
    //std::vector<float> vTracksQPt_EE_[nEE];
    //std::vector<float> vTracks_EE_[nEE];
    std::vector<float> vTracksPt_EB_;
    std::vector<float> vTracksQPt_EB_;
    std::vector<float> vTracks_EB_;

    std::vector<float> vPho_pT_;
    std::vector<float> vPho_E_;
    std::vector<float> vPho_eta_;
    std::vector<float> vPho_phi_;
    std::vector<float> vPho_ecalEPostCorr_;

    std::vector<float> vPho_r9_;
    std::vector<float> vPho_sieie_;
    std::vector<float> vPho_phoIso_;
    std::vector<float> vPho_chgIso_;
    std::vector<float> vPho_chgIsoWrongVtx_;
    std::vector<float> vPho_Eraw_;
    std::vector<float> vPho_phiWidth_;
    std::vector<float> vPho_etaWidth_;
    std::vector<float> vPho_scEta_;
    std::vector<float> vPho_sieip_;
    std::vector<float> vPho_s4_;
    std::vector<float> vPho_rho_;

    std::vector<float> vPho_neuIso_;
    std::vector<float> vPho_ecalIso_;
    std::vector<float> vPho_trkIso_;
    std::vector<float> vPho_hasPxlSeed_;
    std::vector<float> vPho_passEleVeto_;
    std::vector<float> vPho_HoE_;
    std::vector<float> vPho_phoIsoCorr_;
    std::vector<float> vPho_ecalIsoCorr_;
    std::vector<float> vPho_neuIsoCorr_;
    std::vector<float> vPho_chgIsoCorr_;
    std::vector<float> vPho_bdt_;

    std::vector<float> vSC_mass_;
    std::vector<float> vSC_DR_;
    std::vector<float> vSC_E_;
    std::vector<float> vSC_pT_;
    std::vector<float> vSC_eta_;
    std::vector<float> vSC_phi_;

    std::vector<float> vA_E_;
    std::vector<float> vA_pT_;
    std::vector<float> vA_eta_;
    std::vector<float> vA_phi_;
    std::vector<float> vA_mass_;
    std::vector<float> vA_DR_;
    std::vector<float> vA_recoIdx_;
    std::vector<float> vA_pdgId_;
    std::vector<float> vA_mothPdgId_;
    std::vector<float> vA_jetM_;
    std::vector<float> vA_status_;
    float mHgen_;

    std::vector<float> vOutPart_pdgId_;

    std::vector<float> vJet_energy_;
    std::vector<float> vJet_pt_;
    std::vector<float> vJet_eta_;
    std::vector<float> vJet_phi_;
    std::vector<float> vJet_tightId_;

    int nTotal, nPreselPassed, nPassed;
    TH1F * hNpassed_kin;
    TH1F * hNpassed_presel;
    TH1F * hNpassed_mGG;
    TH1F * hNpassed_nRecoPho;
    TH1F * hNpassed_hlt;
    TH1F * hNpassed_img;

    //TProfile2D * hnPho;
    TH2F * hnPho;
    TH2F * hnPhoGt2;
    TH1F * hdR_nPhoGt2;
    TH2F * hdPhidEta_nPhoGt2;
    TProfile2D * hdPhidEta_jphoPt_o_iphoPt;
    TH1F * hjphoPt_o_iphoPt;
    TH1F * hMinDRgenRecoPho;
    TH1F * hMinDRrecoPtoGenPt;
    TH1F * hJetNeuM;

    float m0_;
    std::vector<float> vFC_inputs_;
    int hltAccept_;
    unsigned int nRecoPho_;
    std::vector<float> vMinDR_;
    double evtWeight_;

};

//
// constants, enums and typedefs
//
static const float zs = 0.;

static const int crop_size = 32;
static const bool debug = true;
//static const bool debug = false;

static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;
static const int EB_IETA_MIN = EBDetId::MIN_IETA;//1;
static const int EB_IETA_MAX = EBDetId::MAX_IETA;//85;
static const int EE_MIN_IX = EEDetId::IX_MIN;//1;
static const int EE_MIN_IY = EEDetId::IY_MIN;//1;
static const int EE_MAX_IX = EEDetId::IX_MAX;//100;
static const int EE_MAX_IY = EEDetId::IY_MAX;//100;
static const int EE_NC_PER_ZSIDE = EEDetId::IX_MAX*EEDetId::IY_MAX; // 100*100
static const int HBHE_IETA_MAX_FINE = 20;
static const int HBHE_IETA_MAX_HB = hcaldqm::constants::IETA_MAX_HB;//16;
static const int HBHE_IETA_MIN_HB = hcaldqm::constants::IETA_MIN_HB;//1
static const int HBHE_IETA_MAX_HE = hcaldqm::constants::IETA_MAX_HE;//29;
static const int HBHE_IETA_MAX_EB = hcaldqm::constants::IETA_MAX_HB + 1; // 17
static const int HBHE_IPHI_NUM = hcaldqm::constants::IPHI_NUM;//72;
static const int HBHE_IPHI_MIN = hcaldqm::constants::IPHI_MIN;//1;
static const int HBHE_IPHI_MAX = hcaldqm::constants::IPHI_MAX;//72;
static const int ECAL_IETA_MAX_EXT = 140;

// EE-(phi,eta) projection eta edges
// These are generated by requiring 5 fictional crystals
// to uniformly span each HCAL tower in eta (as in EB).
static const double eta_bins_EEm[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                  {-3.    , -2.93  , -2.86  , -2.79  , -2.72  , -2.65  , -2.62  ,
                   -2.59  , -2.56  , -2.53  , -2.5   , -2.4644, -2.4288, -2.3932,
                   -2.3576, -2.322 , -2.292 , -2.262 , -2.232 , -2.202 , -2.172 ,
                   -2.1462, -2.1204, -2.0946, -2.0688, -2.043 , -2.0204, -1.9978,
                   -1.9752, -1.9526, -1.93  , -1.91  , -1.89  , -1.87  , -1.85  ,
                   -1.83  , -1.812 , -1.794 , -1.776 , -1.758 , -1.74  , -1.7226,
                   -1.7052, -1.6878, -1.6704, -1.653 , -1.6356, -1.6182, -1.6008,
                   -1.5834, -1.566 , -1.5486, -1.5312, -1.5138, -1.4964, -1.479 }; // 56
// EE+(phi,eta) projection eta edges
static const double eta_bins_EEp[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                   {1.479 ,  1.4964,  1.5138,  1.5312,  1.5486,  1.566 ,  1.5834,
                    1.6008,  1.6182,  1.6356,  1.653 ,  1.6704,  1.6878,  1.7052,
                    1.7226,  1.74  ,  1.758 ,  1.776 ,  1.794 ,  1.812 ,  1.83  ,
                    1.85  ,  1.87  ,  1.89  ,  1.91  ,  1.93  ,  1.9526,  1.9752,
                    1.9978,  2.0204,  2.043 ,  2.0688,  2.0946,  2.1204,  2.1462,
                    2.172 ,  2.202 ,  2.232 ,  2.262 ,  2.292 ,  2.322 ,  2.3576,
                    2.3932,  2.4288,  2.4644,  2.5   ,  2.53  ,  2.56  ,  2.59  ,
                    2.62  ,  2.65  ,  2.72  ,  2.79  ,  2.86  ,  2.93  ,  3.    }; // 56

// HBHE eta bin edges
static const double eta_bins_HBHE[2*(hcaldqm::constants::IETA_MAX_HE-1)+1] =
                  {-3.000, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
                   -1.218, -1.131, -1.044, -0.957, -0.870, -0.783, -0.695, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000,
                    0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,  0.695,  0.783,  0.870,  0.957,  1.044,  1.131,  1.218,
                    1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,  1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  3.000}; // 57


//
// static data member definitions
//

#endif
