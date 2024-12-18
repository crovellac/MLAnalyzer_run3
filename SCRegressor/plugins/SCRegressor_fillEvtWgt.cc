#include "MLAnalyzer_run3/SCRegressor/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesEvtWgt ( TTree* tree, edm::Service<TFileService> &fs )
{

  tree->Branch("evt_weight", &evtWeight_);

} // branchesEvtWgt()

// Fill EvtWeight _________________________________________________________________//
void SCRegressor::fillEvtWgt ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  // from: https://github.com/cms-analysis/flashgg/blob/master/Taggers/interface/CollectionDumper.h

  std::string LHEWeightName = "";
  int LHEWeightIndex = -1;

  double weight = 1.;
  double lumiWeight_ = 1.;

  if( ! iEvent.isRealData() ) {

      edm::Handle<GenEventInfoProduct> genInfo;
      iEvent.getByToken(genInfoT_, genInfo);

      weight = lumiWeight_;

      if( LHEWeightName != ""){
          edm::Handle<LHEEventProduct> product_lhe;
          iEvent.getByToken(lheEventT_, product_lhe);
          if( LHEWeightIndex < 0 ){
              for(uint wgt_id = 0 ; wgt_id < product_lhe->weights().size() ; wgt_id++){
                  auto wgt = product_lhe->weights()[wgt_id] ;
                  if( wgt.id == LHEWeightName ){
                      LHEWeightIndex = wgt_id ;
                  }
              }
              std::cout << "Lumi Weight : " << lumiWeight_ << "; LHEWeightIndex: " << LHEWeightIndex << "; LHEWeightName: " << LHEWeightName << std::endl;
          }
          if( LHEWeightIndex > -1 )
              weight *= ( product_lhe->weights()[LHEWeightIndex].wgt/product_lhe->originalXWGTUP () );
      }

      if( genInfo.isValid() ) {
          const auto &weights = genInfo->weights();
          //std::cout << "evt weight:" << genInfo->weights()[0] << std::endl;
          // FIXME store alternative/all weight-sets
          if( ! weights.empty() ) {
              weight *= weights[0];
          }
      }

      /*
      if( globalVarsDumper_ && globalVarsDumper_->puReWeight() ) {
          if (globalVarsDumper_->cache().puweight > 999999. || globalVarsDumper_->cache().puweight < -999999.) {
              weight = 0.;
              std::cout << "WARNING we got a puweight of " << globalVarsDumper_->cache().puweight 
                        << " for rho=" << globalVarsDumper_->cache().rho
                        << " nvtx=" << globalVarsDumper_->cache().nvtx
                        << " npu=" << globalVarsDumper_->cache().npu
                        << " so we set weight to 0!!!" << std::endl;

          } else {
              weight *= globalVarsDumper_->cache().puweight;
          }
      }
      */
  }
  evtWeight_ = weight;
  //std::cout << "evt weight:" << evtWeight_ << std::endl;

} // fillEvtWgt()
