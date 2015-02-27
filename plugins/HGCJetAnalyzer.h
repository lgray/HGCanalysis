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
  std::string genJetsSource_, pfJetsSource_;
};
 

#endif
