#ifndef _HGCHitsAnalyzer_h_
#define _HGCHitsAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "SimG4CMS/Calo/interface/HGCNumberingScheme.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"



#include "UserCode/HGCanalysis/interface/HGCSimulationEvent.h"




#include "TH1F.h"
#include "TTree.h"

#include <string>

/**
   @class HGCHitsAnalyzer
   @author P. Silva (CERN)
*/

class HGCHitsAnalyzer : public edm::EDAnalyzer 
{
  
 public:
  
  explicit HGCHitsAnalyzer( const edm::ParameterSet& );
  ~HGCHitsAnalyzer();
  virtual void analyze( const edm::Event&, const edm::EventSetup& );

 private:

  /**
     @short loops over genparticles and saves a summary of stable (status=1) particles incoming to the detector
   */
  void analyzeGenParticles(edm::Handle<edm::View<reco::Candidate> > &genParticles,edm::Handle<reco::GenJetCollection> &genJets);

  /**
     @short accumulate sim hits
   */
  //void analyzeHits(size_t isd,edm::Handle<edm::PCaloHitContainer> &caloHits,const HGCalGeometry *geom);
  
  //tree and summary ntuple
  TTree *t_;
  HGCSimEvent_t simEvt_;
  
  //gen level
  std::string genSource_, genJetsSource_;
  
  //hgcal
  std::vector<std::string> geometrySource_;
  std::vector<std::string> hitCollections_;
  std::vector<double> mipEn_;
};
 

#endif
