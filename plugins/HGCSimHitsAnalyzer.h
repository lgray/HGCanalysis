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
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "TH1F.h"
#include "TTree.h"
#include "TString.h"

#include <memory>
#include <string>

/**
   @short useful to sort clusters by decreasing energy
*/
bool sortClustersByEnergy(const reco::PFCluster *a, const reco::PFCluster *b) { return (a->energy()>b->energy()); }

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
  math::XYZVectorD getInteractionPosition(const reco::GenParticle & genp, edm::Handle<edm::SimTrackContainer> &SimTk, edm::Handle<edm::SimVertexContainer> &SimVtx, int barcode, Bool_t &hasChargedInteraction);

  //
  inline void resetCounters()
  {
    genId_=0;
    genEn_=0;   genEta_=0;  genPhi_=0;
    genHitX_=0; genHitY_=0; genHitZ_=0;
    hasChargedInteraction_=false;
    hasInteractionBeforeHGC_=false;
    pfMatchId_=0;
    for(std::map<TString, Float_t * >::iterator keyIt=edeps_.begin();
	keyIt!=edeps_.end();
	keyIt++)
      {
	TString key(keyIt->first);
	hitMax_[key]=0;
	hitMaxX_[key]=0; hitMaxY_[key]=0; hitMaxEta_[key]=0; hitMaxPhi_[key]=0;	hitMaxLayer_[key]=0;
	showerMeanX_[key]=0; showerMeanY_[key]=0; showerMeanZ_[key]=0; showerMeanEta_[key]=0; showerMeanPhi_[key]=0;
	nClusters_[key]=0;
	totalE_[key]=0;
	totalX0WgtE_[key]=0;
	totalLambdaWgtE_[key]=0;
	totalLength_[key]=0;
	totalVolume_[key]=0;
	showerStart_[key]=-1;
	for(size_t iclu=0; iclu<5; iclu++)
	  {
	    clusterEn_[key][iclu]=0;
	    clusterZ_[key][iclu]=0;
	    clusterEta_[key][iclu]=0;
	    clusterPhi_[key][iclu]=0;
	  }
	for(size_t ilay=0; ilay<100; ilay++)
	  {
	    nhits_[key][ilay]=0;          nhits5mip_[key][ilay]=0;      nhits10mip_[key][ilay]=0;
	    edeps_[key][ilay]=0;          edeps3x3_[key][ilay]=0;       edeps5x5_[key][ilay]=0;
	    if(ctrledeps_.find(key)!=ctrledeps_.end()) { 
	      ctrlnhits_[key][ilay]=0;
	      ctrledeps_[key][ilay]=0;
	    }
	    emeanPhi_[key][ilay]=0;       emeanEta_[key][ilay]=0;       emeanX_[key][ilay]=0;       emeanY_[key][ilay]=0;
	    edepdR_[key][ilay]=0;         edepArea_[key][ilay]=0;        sihih_[key][ilay]=0;        sipip_[key][ilay]=0;        sipih_[key][ilay]=0;
	  }
      }
  }
  
  inline float getLayerWeight(int layerIdx,bool X0)
  {
    if(layerIdx==0)                  return X0 ? 0.080 : 0.01;
    if(layerIdx>=1 && layerIdx<=10)  return X0 ? 0.620 : 0.036;
    if(layerIdx>=11 && layerIdx<=20) return X0 ? 0.809 : 0.043;
    if(layerIdx>=21 && layerIdx<=29) return X0 ? 1.239 : 0.056;
    if(layerIdx==30)                 return X0 ? 3.580 : 0.338;
    if(layerIdx>=31 && layerIdx<=41) return X0 ? 3.103 : 0.273;
    if(layerIdx>=42 && layerIdx<=53) return X0 ? 5.228 : 0.475;
    return 0.;
  }

  //variables to store in tree
  Int_t genId_;
  Float_t genEn_,genEta_,genPhi_;
  Float_t genHitX_, genHitY_, genHitZ_;
  Int_t pfMatchId_;
  Bool_t hasChargedInteraction_,hasInteractionBeforeHGC_;
  Int_t nlay_;
  std::map<TString, Float_t> showerMeanX_, showerMeanY_, showerMeanZ_, showerMeanEta_, showerMeanPhi_;
  std::map<TString,Int_t> nClusters_;
  std::map<TString, Float_t *> clusterEn_, clusterZ_, clusterEta_, clusterPhi_;
  std::map<TString, Float_t> hitMax_, hitMaxX_, hitMaxY_, hitMaxEta_, hitMaxPhi_;
  std::map<TString, Int_t> hitMaxLayer_, showerStart_;
  std::map<TString, Float_t> totalE_, totalX0WgtE_, totalLambdaWgtE_, totalLength_, totalVolume_;
  std::map<TString, Float_t *> edeps_, weightedEdeps_, edeps3x3_, edeps5x5_, ctrledeps_;
  std::map<TString, Int_t *> nhits_, nhits5mip_, nhits10mip_, ctrlnhits_;
  std::map<TString, Float_t *> emeanX_, emeanY_, emeanPhi_,    emeanEta_;
  std::map<TString, Float_t *> sihih_,        sipip_,        sipih_,        edepdR_,        edepArea_;

  //tree and summary ntuple
  TTree *t_;

  //mip energy for each sensitive detector and thresholds to apply
  std::vector<double> mipEn_;
 
  //gen level
  std::string genSource_;
  
  //hgcal
  std::vector<std::string> hitCollections_, recHitCollections_, geometrySource_;
  std::string pfClustersCollection_, emPFClustersCollection_;

  //Geant4
  std::string g4TracksSource_, g4VerticesSource_;

  //association to gen candidate
  double pfCandAssociationCone_,pfClusterAssociationCone_;
};
 

#endif
