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

#include "UserCode/HGCanalysis/interface/SlimmedRecHit.h"
#include "UserCode/HGCanalysis/interface/SlimmedJet.h"
#include "UserCode/HGCanalysis/interface/SlimmedVertex.h"
#include "UserCode/HGCanalysis/interface/SlimmedCluster.h"

#include "TTree.h"
#include "TVector3.h"

#include <unordered_map>

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
  
  void slimRecHits(const edm::Event &iEvent, const edm::EventSetup &iSetup);
  void doMCJetMatching(edm::Handle<std::vector<reco::PFJet> > &pfJets,
		       edm::Handle<reco::GenJetCollection> &genJets,
		       edm::Handle<edm::View<reco::Candidate> > &genParticles,
		       std::unordered_map<uint32_t,uint32_t> &reco2genJet,
		       std::unordered_map<uint32_t,uint32_t> &genJet2Parton,
		       std::unordered_map<uint32_t,uint32_t> &genJet2Stable);

  virtual void endJob() ;

  std::string eeSimHitsSource_, hefSimHitsSource_;
  std::string eeRecHitsSource_, hefRecHitsSource_;

  TTree *tree_;
  Int_t run_,event_,lumi_;
  std::vector<SlimmedRecHit> *slimmedRecHits_;
  std::vector<SlimmedCluster> *slimmedClusters_;
  std::vector<SlimmedJet> *slimmedJets_;
  std::vector<SlimmedVertex> *slimmedVertices_;
  TVector3 *genVertex_;

  std::string g4TracksSource_, g4VerticesSource_;
  std::string recoVertexSource_;
  std::string genSource_, genJetsSource_, pfJetsSource_;
};
 

#endif
