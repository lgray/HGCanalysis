#include "UserCode/HGCanalysis/plugins/HGCSimHitsAnalyzer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "SimG4CMS/Calo/interface/CaloHitID.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <iostream>

using namespace std;

//
HGCSimHitsAnalyzer::HGCSimHitsAnalyzer( const edm::ParameterSet &iConfig )
{
  //configure analyzer
  hitCollections_   = iConfig.getUntrackedParameter< std::vector<std::string> >("hitCollections");
  geometrySource_   = iConfig.getUntrackedParameter< std::vector<std::string> >("geometrySource");

  //init tree
  edm::Service<TFileService> fs;
  t_=fs->make<TTree>("HGC","Event Summary");
  t_->Branch("genId",  &genId_,  "genId/I");
  t_->Branch("genEn",  &genEn_,  "genEn/F");  
  t_->Branch("genEta", &genEta_, "genEta/F");
  t_->Branch("genPhi", &genPhi_, "genPhi/F");
  t_->Branch("nlay",   &nlay_,   "nlay/I");
  t_->Branch("edeps",  edeps_,   "edeps[nlay]/F");
}

//
HGCSimHitsAnalyzer::~HGCSimHitsAnalyzer()
{
}

//
void HGCSimHitsAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  //generator level particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);
  if(genParticles->size()>1) std::cout << "[Warning] found more than 1 gen particle, will save only first one" << std::endl;
  const reco::GenParticle & p = (*genParticles)[0];
  genId_  = p.pdgId();
  genEn_  = p.energy();
  genEta_ = p.eta();
  genPhi_ = p.phi();
    
  //hits + geometry
  std::map<int, std::vector<float> > allEdeps;
  for(size_t i=0; i<hitCollections_.size(); i++)
    {
      edm::Handle<edm::PCaloHitContainer> caloHits;
      iEvent.getByLabel(edm::InputTag("g4SimHits",hitCollections_[i]),caloHits); 
      edm::ESHandle<HGCalGeometry> geom;
      iSetup.get<IdealGeometryRecord>().get(geometrySource_[i],geom);
      const HGCalTopology &topo=geom->topology();
      const HGCalDDDConstants &dddConst=topo.dddConstants();
      std::vector<float> edepsPerLayer(dddConst.layers(true),0.);
      uint32_t mySubDet(ForwardSubdetector::HGCEE);
      if(i==1) mySubDet=ForwardSubdetector::HGCHEF;
      else if(i==2) mySubDet=ForwardSubdetector::HGCHEB;

      for(edm::PCaloHitContainer::const_iterator hit_it = caloHits->begin(); hit_it != caloHits->end(); ++hit_it) 
	{
	  //gang SIM->RECO cells to get final layer assignment
	  HGCalDetId simId(hit_it->id());
	  int layer(simId.layer()),cell(simId.cell());
	  std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,topo.detectorType());
	  cell  = recoLayerCell.first;
	  layer = recoLayerCell.second;
	  if(layer<0) continue;

	  //get global position
	  uint32_t recoDetId( (i==0) ?
			      (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell) :
			      (uint32_t)HGCHEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell)
			      );
	  const GlobalPoint pos( std::move( geom->getPosition(recoDetId) ) );	  
	  float hitEta(pos.eta());
	  if(hitEta*genEta_<0) continue;    
	  
	  //save it if interesting
	  edepsPerLayer[layer-1]+=hit_it->energy()*TMath::TanH(hitEta);
	}
      allEdeps[i]=edepsPerLayer;
    }
  
  nlay_=0;
  for(std::map<int, std::vector<float> >::iterator it = allEdeps.begin(); it != allEdeps.end(); it++)
    {
      for(std::vector<float>::iterator jt=it->second.begin(); jt!=it->second.end(); jt++)
	{
	  edeps_[nlay_]=*jt;
	  nlay_++;
	}
    }
  
  //all done, save information
  t_->Fill();
}


//define this as a plug-in
DEFINE_FWK_MODULE(HGCSimHitsAnalyzer);
