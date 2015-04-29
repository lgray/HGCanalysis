#ifndef _HGCJetAnalyzer_h_
#define _HGCJetAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TTree.h"
#include "TH2F.h"
#include "TString.h"

#include <map>

/**
   @class HGCJetAnalyzer
   @author P. Silva (CERN)
*/

class HGCJetAnalyzer : public edm::EDAnalyzer 
{  
 public:
  
  explicit HGCJetAnalyzer( const edm::ParameterSet& );
  ~HGCJetAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:

  virtual void endJob() ;

  std::map<TString, TH2F *> histMap_;
  TTree *jetTree_;
  Float_t jpt_,  jeta_,  jphi_,  jmass_,jnhf_,jnhe_,jnhm_,jgf_,jge_,jgm_,jchf_,jche_,jchm_;
  Float_t gjpt_, gjeta_, gjphi_, gjmass_;
  Float_t ppt_,  peta_,  pphi_,  pmass_, pid_;
  Int_t nTDCHits_;
  std::string genSource_, genJetsSource_, pfJetsSource_, eeRecHitsSource_, hefRecHitsSource_;
};
 

#endif
