#include "MLAnalyzer_run3/SCRegressor/interface/SCRegressor.h"

struct gen_obj {
  double E;
  double pT;
  double phi;
  double eta;
  double mass;
  double dR;
};

std::vector<gen_obj> vAs;
std::vector<unsigned int> vGenEleIdxs;

// Initialize branches _____________________________________________________//
void SCRegressor::branchesH2aaSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("mHgen",     &mHgen_);
  //tree->Branch("FC_inputs", &vFC_inputs_);
  //tree->Branch("hltAccept", &hltAccept_);

  tree->Branch("A_mass",    &vA_mass_);
  tree->Branch("A_DR",      &vA_DR_);
  tree->Branch("A_E",       &vA_E_);
  tree->Branch("A_pT",      &vA_pT_);
  tree->Branch("A_eta",     &vA_eta_);
  tree->Branch("A_phi",     &vA_phi_);
  tree->Branch("A_recoIdx", &vA_recoIdx_);

  hdPhidEta = fs->make<TH2F>("dPhidEta_GG", "#Delta(#phi,#eta,m);#Delta#phi(#gamma,#gamma);#Delta#eta(#gamma,#gamma)",
              6, 0., 6.*0.0174, 6, 0., 6.*0.0174);
}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runH2aaSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<ElectronCollection> electrons;
  iEvent.getByToken(electronCollectionT_, electrons);  

  vAs.clear();
  vGenEleIdxs.clear();
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vH;
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );
    //We pick out an electron whose mother is the pseudoscalar which has two daughters.
    //If this is the case, we save the index of this electron, and add the mother to our vector.
    if (std::abs(iGen->pdgId()) != 11) continue;
    if (std::abs(iGen->mother()->pdgId()) != 25) continue;
    if (iGen->mother()->numberOfDaughters() != 2) continue;
    if (std::abs(iGen->mother()->daughter(0)->pdgId()) != 11 || std::abs(iGen->mother()->daughter(1)->pdgId()) != 11) continue;
    //Now be sure that we're not double-counting a pseudoscalar
    bool alreadyCounted = false;
    for ( unsigned int iA = 0; iA < vAs.size(); iA++) {
      double sampleEta = vAs[iA].eta;
      double samplePhi = vAs[iA].phi;
      if (iGen->mother()->eta() == sampleEta && iGen->mother()->phi() == samplePhi) alreadyCounted = true;
    }
    //Only add the pseudoscalar if it's not already counted
    if (!alreadyCounted) {
      float dR = reco::deltaR(iGen->mother()->daughter(0)->eta(),iGen->mother()->daughter(0)->phi(), iGen->mother()->daughter(1)->eta(),iGen->mother()->daughter(1)->phi());
      gen_obj MotherA = { iGen->mother()->energy(), iGen->mother()->pt(), iGen->mother()->phi(), iGen->mother()->eta(), iGen->mother()->mass(), dR };
      vAs.push_back( MotherA );
    } 
    vGenEleIdxs.push_back(iG);
  } // gen particles
  if ( vAs.size() != 2) return false;
  
  //   Find the indices of the Reco Electrons which match to the Gen Electrons we've selected.
   
  float minDR = 100.;
  int minDR_idx = -1;
  //Loop over the gen electrons we just picked out.
  for ( auto& iG : vGenEleIdxs ) {
    reco::GenParticleRef iGenEle( genParticles, iG );
    
    //Now loop over all of the reco electrons
    minDR = 100.;
    minDR_idx = -1;
    for ( unsigned int iP = 0; iP < electrons->size(); iP++ ) {
      ElectronRef iRecoEle( electrons, iP );

      double dR = reco::deltaR( iRecoEle->eta(),iRecoEle->phi(), iGenEle->eta(),iGenEle->phi() );
      //std::cout << "dR: " << dR << std::endl;
      if ( dR > minDR ) continue;

      minDR = dR;
      minDR_idx = iP;      

    }
    //We now have the index of the reco electron which has the lowest dR with this gen electron.

    if ( minDR > 0.04 ) continue;

    //Check to see if this index is already present in 
    bool alreadyPresent = false;
    std::cout << "minDR_idx: " << minDR_idx << std::endl;
    for ( unsigned int i = 0; i < vPreselPhoIdxs_.size(); i++ ) {
      if (minDR_idx == vPreselPhoIdxs_[i]) alreadyPresent = true;
    }
    if (!alreadyPresent) {
      vPreselPhoIdxs_.push_back(minDR_idx);
    }    
  } 
    
  return true;
}

// Fill branches ___________________________________________________________________//
void SCRegressor::fillH2aaSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<ElectronCollection> electrons;
  iEvent.getByToken(electronCollectionT_, electrons);


  vA_E_.clear();
  vA_pT_.clear();
  vA_eta_.clear();
  vA_phi_.clear();
  vA_mass_.clear();
  vA_DR_.clear();
  vA_recoIdx_.clear();
  
  for ( unsigned int iG = 0; iG < vAs.size(); iG++ ) {
    gen_obj MotherA = vAs[iG];
    vA_E_.push_back( MotherA.E );
    vA_pT_.push_back( MotherA.pT );
    vA_eta_.push_back( MotherA.eta );
    vA_phi_.push_back( MotherA.phi );
    vA_mass_.push_back( MotherA.mass );
    vA_DR_.push_back( MotherA.dR );
  }

  /*
  float dPhi, dEta, dR, recoDR;
  int recoDR_idx;
  for ( unsigned int iG : vGenAIdxs_ ) {

    reco::GenParticleRef iGen( genParticles, iG );

    vA_E_.push_back( std::abs(iGen->energy()) );
    vA_pT_.push_back( std::abs(iGen->pt()) );
    vA_eta_.push_back( iGen->eta() );
    vA_phi_.push_back( iGen->phi() );
    vA_mass_.push_back( iGen->mass() );
    vA_DR_.push_back( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) );

    dPhi = reco::deltaPhi( iGen->daughter(0)->phi(), iGen->daughter(1)->phi() );
    dEta = std::abs( iGen->daughter(0)->eta() - iGen->daughter(1)->eta() );
    hdPhidEta->Fill( dPhi, dEta );

    // Get index to dR-matched preselected photon
    recoDR = 2*0.04;
    recoDR_idx = -1;
    // Want vA_recoIdx_ to store vector index in vRegressPhoIdxs_
    // i.e., vRegressPhoIdxs_[0]:leading reco pho, vRegressPhoIdxs_[1]:sub-leading reco pho
    // not position in original photon collection
    for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
      //PhotonRef iPho( photons, vRegressPhoIdxs_[iP] );
      ElectronRef iEle ( electrons, vRegressPhoIdxs_[iP] );
      dR = reco::deltaR(iGen->eta(),iGen->phi(), iEle->eta(),iEle->phi());
      if ( dR < recoDR ) {
        recoDR = dR;
        recoDR_idx = iP;
      }
    } // reco pho
    vA_recoIdx_.push_back( recoDR_idx );
    
  } // gen A
  */

}
