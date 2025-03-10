#include "MLAnalyzer_run3/SCRegressor/interface/SCRegressor.h"

// Fill EB rec hits ////////////////////////////////
// Store event rechits in a vector of length equal
// to number of crystals in EB (ieta:170 x iphi:360)

//TProfile2D *hEB_energy;
//TProfile2D *hEB_time;
//std::vector<float> vEB_energy_;
//std::vector<float> vEB_time_;

// Initialize branches _____________________________________________________//
void SCRegressor::branchesEE ( TTree* tree, edm::Service<TFileService> &fs ) {
  char hname[50], htitle[50];
  for ( int iz(0); iz < nEE; iz++ ) {
    // Branches for images
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "EE%s_energy",zside);
    tree->Branch(hname,        &vEE_energy_[iz]);
    sprintf(hname, "EE%s_time",  zside);
    tree->Branch(hname,        &vEE_time_[iz]);

    // Histograms for monitoring
    sprintf(hname, "EE%s_energy",zside);
    sprintf(htitle,"E(ix,iy);ix;iy");
    hEE_energy[iz] = fs->make<TProfile2D>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
    sprintf(hname, "EE%s_time",zside);
    sprintf(htitle,"t(ix,iy);ix;iy");
    hEE_time[iz] = fs->make<TProfile2D>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
  } // iz

} // branchesEE()

// Fill EE rechits _________________________________________________________________//
void SCRegressor::fillEE ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  std::cout << "Starting fillEE" << std::endl;
  int ix_, iy_, iz_, idx_; // rows:ieta, cols:iphi
  float energy_;

  for ( int iz(0); iz < nEE; iz++ ) {
    vEE_energy_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vEE_time_[iz].assign( EE_NC_PER_ZSIDE, 0. );
  }

  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByToken( EERecHitCollectionT_, EERecHitsH_ );
  
  std::cout << "EERecHitsH_->size(): " << EERecHitsH_->size() << std::endl;

  // Fill EE rechits
  for ( EcalRecHitCollection::const_iterator iRHit = EERecHitsH_->begin();
        iRHit != EERecHitsH_->end(); ++iRHit ) {
    energy_ = iRHit->energy();
    if ( energy_ <= zs ) continue;
    // Get detector id and convert to histogram-friendly coordinates
    EEDetId eeId( iRHit->id() );
    ix_ = eeId.ix() - 1;
    iy_ = eeId.iy() - 1;
    iz_ = (eeId.zside() > 0) ? 1 : 0;
    // Fill histograms for monitoring 
    hEE_energy[iz_]->Fill( ix_, iy_, energy_ );
    hEE_time[iz_]->Fill( ix_, iy_, iRHit->time() );
    // Get Hashed Index: provides convenient index mapping from [iy][ix] -> [idx]
    idx_ = iy_*EE_MAX_IX + ix_;
    // Fill vectors for images
    vEE_energy_[iz_][idx_] = energy_;
    vEE_time_[iz_][idx_] = iRHit->time();


  } // EE rechits

} // fillEE()
