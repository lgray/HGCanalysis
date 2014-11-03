#ifndef _HGCSimHitsAnalyzer_h_
#define _HGCSimHitsAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "UserCode/HGCanalysis/interface/HGCSimulationEvent.h"

#include "TH1F.h"
#include "TTree.h"
#include "TString.h"

#include <memory>
#include <string>

/**
   @class HGCSimHitsAnalyzer
   @author P. Silva (CERN)
*/

class HGCSimHitsAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCSimHitsAnalyzer( const edm::ParameterSet& );
  ~HGCSimHitsAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:
  
  //
  inline void resetCounters()
  {
    for(std::map<TString, std::vector< Float_t *> >::iterator keyIt=edeps_.begin();
	keyIt!=edeps_.end();
	keyIt++)
      {
	TString key(keyIt->first);
	for(size_t ithr=0; ithr<thrList_.size(); ithr++)
	  for(size_t ilay=0; ilay<100; ilay++)
	    {
	      edeps_[key][ithr][ilay]=0;
	      emeanPhi_[key][ithr][ilay]=0;
	      emeanEta_[key][ithr][ilay]=0;
	      sihih_[key][ithr][ilay]=0;
	      sipip_[key][ithr][ilay]=0;
	      sipih_[key][ithr][ilay]=0;
	      nhits_[key][ithr][ilay]=0;
	    }
      }
  }
  
  Int_t genId_;
  Float_t genEn_,genEta_,genPhi_;
  Int_t nlay_;
  std::map<TString, std::vector< Float_t *> > edeps_,    emeanPhi_,    emeanEta_,    sihih_,     sipip_,     sipih_;
  std::map<TString, std::vector< Int_t *> >   nhits_;

  //tree and summary ntuple
  TTree *t_;

  //mip energy for each sensitive detector and thresholds to apply
  std::vector<double> mipEn_, thrList_;
 
  //gen level
  std::string genSource_;
  
  //hgcal
  std::vector<std::string> hitCollections_, recHitCollections_, pfClustersCollections_, geometrySource_;
};
 

#endif
